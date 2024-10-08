# Data reading and flat inputs function
function data_read_function(File_name::String,dk,path = "inst/")
    df = readdlm("$path$File_name",',','\n');
    C = df[2,2];
    m = df[3,2];
    N = C + 1;
    Fnode2=[fnode2(i,"name",0,0,0,0,0,0,0,0,0,0) for i=1:N]
    Fnode2[1]=fnode2(1,"name",0,df[6,2],df[6,3],df[6,2]+df[6,3],0,0,0,0,0,0)
    for i in 1:C
        idi = i+1;
        namei = "name";
        pp = df[8+i,2];
        agi = df[8+i,3];
        coti = df[8+i,4]
        bgi = agi+coti;
        t_max_ti = df[8+i,6]
        delta_gi = df[8+i,7]- t_max_ti;        
        temp_ci = df[8+i,5]
        temp_si = temp_ci + 5
        ntri = cld(coti,dk)
        if dk >15
            TMAXi = delta_gi + t_max_ti - dk
        else
            TMAXi = delta_gi + t_max_ti
        end
        Fnode2[i+1] = fnode2(idi,namei,pp,agi,coti,bgi,delta_gi,t_max_ti,temp_ci,temp_si,ntri,TMAXi)
    end
    ttm=zeros(Int,N,N)
    for i=1:N*N
        ii = df[9+N+i,1]+1
        jj = df[9+N+i,2]+1
        ttm[ii,jj] = df[9+N+i,3]
    end
    for i=2:N
     ttm[i,1]=ttm[i,1]+Fnode2[i].temp_c
    end
    for i=1:N
        for j=2:N
            if i!=j
                ttm[i,j]=ttm[i,j]+Fnode2[j].temp_c
            end
        end
    end
    dm = zeros(Float64,N,N)
    for i=1:N*N
        ii = df[11+N+N*N+i,1]+1
        jj = df[11+N+N*N+i,2]+1
        dm[ii,jj] = df[11+N+N*N+i,3]
    end
    Fnode=[fnode(0,"name",0) for i=1:N]
    Fnode[1]=fnode(1,"name",0)
        for i=1:N-1
            Fnode[i+1]=fnode(Fnode2[i+1].id,Fnode2[i+1].name,Fnode2[i+1].ntr)
        end
    K=sum(Fnode[i].ntr for i=2:N)
    Kdata=[]
    lab_initial_time=[]
    lab_final_time=[]
    idk=1
    for i=1:N-1
        op=Fnode2[i+1].ag
        closetime=Fnode2[i+1].bg
        push!(lab_initial_time,op-ttm[1,i+1])
        for j=Fnode[i+1].ntr-1:-1:0
            push!(Kdata,commk(Fnode[i+1],closetime-j*dk,(closetime-j*dk)+Fnode2[i+1].TMAX,idk))
            idk+=1
        end
        push!(lab_final_time,last(Kdata).td)
    end
    initial_time=minimum(lab_initial_time)
    final_time=maximum(lab_final_time)
    return N,Fnode,ttm,Kdata,K,initial_time,final_time,m,Fnode2,dm

end
#------------------------------------------------------------------------------
# dt approximation function
function dt_approximation_function(N,ttm,dt,Kdata,initial_time,final_time)
    ttm_dt=zeros(Int,N,N)
    for i=1:N
        for j=1:N
            ttm_dt[i,j]=Int64(cld(ttm[i,j],dt)*dt)
        end
    end
    Kdata_dt=[]
    for i in Kdata
        push!(Kdata_dt,commk(i.ok,Int64(cld(i.tp,dt)*dt),Int64(fld(i.td,dt)*dt),i.id))
    end
    initial_time_dt=Int64(cld(initial_time,dt)*dt)
    final_time_dt=Int64(fld(final_time,dt)*dt)
    TH=final_time_dt-initial_time_dt
    T=Int64(TH/dt)
    return ttm_dt,Kdata_dt,initial_time_dt,final_time_dt,T
end
#------------------------------------------------------------------------------
# time expanded imputs fuction
function inputs_fuction(N::Int64,Kdata::Array{Any,1},ttm::Array{Int64,2},Fnode2,Fnode,initial_time_dt,initial_time)
	SKdata=[]
    idsk=1
    for i=1:N-1
        idtr=findall(j -> j.ok.id==i+1,Kdata)
        tp=minimum(Kdata[ii].tp for ii in idtr)
        td=maximum(Kdata[ii].td for ii in idtr)
        tpr=maximum(Kdata[ii].tp for ii in idtr)
        tdr=minimum(Kdata[ii].td for ii in idtr)
        push!(SKdata,scommk(Fnode[i+1],idtr,tp,td,tpr,tdr,idsk))
        idsk+=1
    end
    groupinginfact=[]
    for i in SKdata
        if i.tdr-i.tpr<ttm[i.ok.id,1]
            push!(groupinginfact,true)
        else
            push!(groupinginfact,false)
        end
    end
    if sum(groupinginfact)!=0
        firstgroupingcondition=0
    else
        firstgroupingcondition=1
    end
    while firstgroupingcondition!=1
        for i in eachindex(SKdata)
            if groupinginfact[i]
                sklist=SKdata[i].list
                ntrsk=cld(length(sklist),2)
                list=collect(sklist[1]:sklist[ntrsk])
                tp=minimum(Kdata[ii].tp for ii in list)
                td=maximum(Kdata[ii].td for ii in list)
                tpr=maximum(Kdata[ii].tp for ii in list)
                tdr=minimum(Kdata[ii].td for ii in list)
                SKdata[i]=scommk(SKdata[i].ok,list,tp,td,tpr,tdr,SKdata[i].id)
                if SKdata[i].tdr-SKdata[i].tpr>ttm[SKdata[i].ok.id,1]
                    groupinginfact[i]=false
                end
                list=collect(sklist[ntrsk+1]:last(sklist))
                tp=minimum(Kdata[ii].tp for ii in list)
                td=maximum(Kdata[ii].td for ii in list)
                tpr=maximum(Kdata[ii].tp for ii in list)
                tdr=minimum(Kdata[ii].td for ii in list)
                idsk=last(SKdata).id+1
                push!(SKdata,scommk(SKdata[i].ok,list,tp,td,tpr,tdr,idsk))
                if last(SKdata).tdr-last(SKdata).tpr<ttm[last(SKdata).ok.id,1]
                    push!(groupinginfact,true)
                else
                    push!(groupinginfact,false)
                end
            end
        end
        if sum(groupinginfact)==0
            firstgroupingcondition=1
        end
    end
    sK=length(SKdata)
    tlab=[]
    push!(tlab,initial_time_dt)
    for k in SKdata
        push!(tlab,k.tp-ttm[1,k.ok.id])
        push!(tlab,k.td)
    end
    tlab=unique(sort(tlab))
	NE=[node(Fnode[1],i,i) for i in tlab]
    for i in eachindex(NE)
        NE[i].id=i
    end
    id=last(NE).id
    L0=NE[1]
    LT=NE[id]
    L03=L0
    for fn=2:N
        id+=1
        push!(NE,node(Fnode[fn],initial_time,id))
		id+=1
		push!(NE,node(Fnode[fn],Fnode2[fn].bg+30,id))
    end
    for k in SKdata
        id+=1
        push!(NE,node(k.ok,k.tp,id))
    end
    Ok=[node(Fnode[1],k,k) for k=1:sK]
    Dk=[node(Fnode[1],k,k) for k=1:sK]
    cont=1
    for k in SKdata
        idne=findall(i -> i.n==k.ok && i.t==k.tp,NE)
        Ok[cont]=NE[idne[1]]
        idne=findall(i -> i.n.id==1 && i.t==k.td,NE)
        Dk[cont]=NE[idne[1]]
        cont+=1
    end
    return NE,SKdata,L0,LT,Ok,Dk,L03
end
function inputs_fuction_2(N::Int64,Kdata_dt::Array{Any,1},ttm::Array{Int64,2},Fnode2)
	SKdata=[]
    for i in Kdata_dt
        idtr=i.id
        tp=i.tp
        td=i.td
        tpr=i.tp
        tdr=i.td
        push!(SKdata,scommk(i.ok,idtr,tp,td,tpr,tdr,i.id))
    end
    sK=length(SKdata)
    tlab=[]
    push!(tlab,initial_time_dt)
    for k in SKdata
        push!(tlab,k.tp-ttm[1,k.ok.id])
        push!(tlab,k.td)
    end
    tlab=unique(sort(tlab))
	NE=[node(Fnode[1],i,i) for i in tlab]
    for i in eachindex(NE)
        NE[i].id=i
    end
    id=last(NE).id
    L0=NE[1]
    LT=NE[id]
    L03=L0
    for fn=2:N
        id+=1
        push!(NE,node(Fnode[fn],initial_time,id))
		id+=1
		push!(NE,node(Fnode[fn],Fnode2[fn].bg+30,id))
    end
    for k in SKdata
        id+=1
        push!(NE,node(k.ok,k.tp,id))
    end
    Ok=[node(Fnode[1],k,k) for k=1:sK]
    Dk=[node(Fnode[1],k,k) for k=1:sK]
    cont=1
    for k in SKdata
        idne=findall(i -> i.n==k.ok && i.t==k.tp,NE)
        Ok[cont]=NE[idne[1]]
        idne=findall(i -> i.n.id==1 && i.t==k.td,NE)
        Dk[cont]=NE[idne[1]]
        cont+=1
    end
    return NE,SKdata,L0,LT,Ok,Dk,L03
end

#------------------------------------------------------------------------------
#expanded_arcs_function
function expanded_arcs_function(NE::Array{node,1},Fnode::Array{fnode,1},ttm::Array{Int64,2},LT::node,SKdata::Array{Any,1},Fnode2)
	arc_ev=[]
    idarc=1
    for i in NE
        if i!=LT
            for j in Fnode
                if j != i.n
                    sel=findall(ii -> ii.n==j && ii.t<=ttm[i.n.id,j.id]+i.t,NE)
                    if sel != Int64[]
                        tm=maximum(NE[ii].t for ii in sel)
                        jj=NE[findall(ii -> ii.n == j && ii.t == tm,NE)]
                        if jj[1].t-i.t==ttm[i.n.id,jj[1].n.id] && i.t+ttm[i.n.id,jj[1].n.id] <= Fnode2[jj[1].n.id].bg + 30
                            push!(arc_ev,arc(i,jj[1],1,idarc))
                            idarc+=1
                        elseif i.t+ttm[i.n.id,jj[1].n.id] <= Fnode2[jj[1].n.id].bg + 30
                            push!(arc_ev,arc(i,jj[1],0,idarc))
                            idarc+=1
                        end
                    end
                end
            end
        end
        sel=findall(ii -> ii.n==i.n && ii.t>i.t,NE)
        if sel != Int64[]
            tm=minimum(NE[ii].t for ii in sel)
            jj=NE[findall(ii -> ii.n == i.n && ii.t == tm,NE)]
            push!(arc_ev,arc(i,jj[1],1,idarc))
            idarc+=1
        end
    end
    arc_ek=[]
    for k = 1:length(SKdata)
        for ij in arc_ev
            if ij.tail.n.id!=1 && ij.head.n.id!=SKdata[k].ok.id
                push!(arc_ek,arck(ij.tail,ij.head,k,ij.id))
            elseif ij.tail.n.id == 1 && ij.head.n.id == 1
                push!(arc_ek,arck(ij.tail,ij.head,k,ij.id))
            elseif ij.head.n.id==SKdata[k].ok.id && ij.tail.n.id==SKdata[k].ok.id
                push!(arc_ek,arck(ij.tail,ij.head,k,ij.id))
            end
        end
    end
    arcsk=[]
    selidarck = [findall(a -> a.k == k, arc_ek) for k in 1 : length(SKdata)]
    for k=1:length(SKdata)
        push!(arcsk,[arc_ek[i] for i in selidarck[k]])
    end
    return arc_ev,arc_ek,arcsk
end
#-----------------------------------------------------------------------------
#model creation function
function model_function(model_timelimit,arc_ev,arc_ek,ttm,NE,L0,LT,SKdata,arcsk,Ok,Dk)
    model=direct_model(CPLEX.Optimizer());
    model_env=backend(model).env
    set_optimizer_attribute(model, "CPX_PARAM_TILIM", model_timelimit)
    set_optimizer_attribute(model, "CPXPARAM_Emphasis_MIP", 5)
    set_optimizer_attribute(model, "CPXPARAM_MIP_Tolerances_MIPGap", 0.03)
    set_optimizer_attribute(model, "CPXPARAM_Threads", 8)
    @variable(model,y[ij in arc_ev]>=0, Int)
    @variable(model,x[ijk in arc_ek],Bin)
    @objective(model, Min, sum( cost3(ij,ttm) * y[ij] for ij in arc_ev))
    for i in NE
        if i.id!=L0.id && i.id!=LT.id
            @constraint(model,sum(y[ij] for ij in arc_ev if ij.tail.id==i.id)
                            -sum(y[ij] for ij in arc_ev if ij.head.id==i.id)==0)
       end
    end
    for k=1:length(SKdata)
        @constraint(model, sum(x[ijk] for ijk in arcsk[k] if ijk.tail==Ok[k]) ==1)
        @constraint(model, sum(x[ijk] for ijk in arcsk[k] if ijk.head==Dk[k]) ==1)
        for i in NE
            if i.id!=Ok[k].id && i.id!=Dk[k].id
                @constraint(model,sum(x[ijk] for ijk in arcsk[k] if ijk.tail.id==i.id)
                             -sum(x[ijk] for ijk in arcsk[k] if ijk.head.id==i.id)==0)
            end
        end
        for ijk in arcsk[k]
            if ijk.tail.n != ijk.head.n
                @constraint(model,x[ijk] <= y[arc_ev[ijk.id]])
            end
            if ijk.tail.n == ijk.head.n && ijk.tail.n.id!=Ok[k].id
                @constraint(model,x[ijk] <= y[arc_ev[ijk.id]])
            end
        end
    end
    for ij in arc_ev
        if ij.tail.n.id !=1
            @constraint(model,y[ij] <= sum(x[ijk] for ijk in arc_ek if ijk.id == ij.id))
        end
    end
    for ij in arc_ev
        if ij.tail.n.id == 1 && ij.head.n.id !=1
            @constraint(model,y[ij]<=1)
        end
    end
    return model,x,y
end
#-----------------------------------------------------------------------------
# partial expanded network solution verification function
function partial_expanded_verification_function(y,arc_ev::Array{Any,1})
    idarcssolinf=findall(ij -> value(y[ij])!=0 && ij.tt==0,arc_ev)
    if idarcssolinf == Int64[]
        condition=false
    else
        condition=true
    end
    idarcneg = findall(ij -> value(y[ij])!=0 && ij.head.t-ij.tail.t < 0,arc_ev)
    negarcs = true
    if idarcneg == Int64[]
        negarcs = false
    end
    return condition,idarcssolinf,length(idarcssolinf), negarcs
end
#-----------------------------------------------------------------------------
function grouping_verification_function(x,Ok::Array{node,1},Dk::Array{node,1},arcsk::Array{Any,1},SKdata::Array{Any,1})
    kpd=[]
    for k=1:length(SKdata)
        idarcssolp=findall(ijk -> value(x[ijk])>0.5 && ijk.tail.n.id==Ok[k].n.id,arcsk[k])
        idarcssold=findall(ijk -> value(x[ijk])>0.5 && ijk.head.n.id==Dk[k].n.id,arcsk[k])
        tps=maximum(arcsk[k][i].tail.t for i in idarcssolp)
        tds=minimum(arcsk[k][i].head.t for i in idarcssold)
        push!(kpd,(tps,tds))
    end
    inf=[]
    for k=1:length(SKdata)
        dtp=SKdata[k].tpr-kpd[k][1]
        dtd=kpd[k][2]-SKdata[k].tdr
        push!(inf,maximum((dtp,dtd)))
    end
    inf2=[]
    for k=1:length(SKdata)
        dtp=SKdata[k].tpr-kpd[k][1]
        dtd=kpd[k][2]-SKdata[k].tdr
        if dtp>dtd
            infa=1
        else
            infa=2
        end
        push!(inf2,(dtp,dtd,infa))
    end
    idinfa=findall(i -> i>0,inf)
    if idinfa == Int64[]
        condition=false
    else
        condition=true
    end
    infamax=findmax(inf)
    return condition,infamax[1],infamax[2],inf2[infamax[2]][3],kpd[infamax[2]]
end
#-----------------------------------------------------------------------------
# grouping correction
function grouping_correction_function3(GV,SKdata,Kdata_dt,NE,Fnode,Dk,Ok,ttm_dt)
    sklist=SKdata[GV[3]].list
    cond3=0
    inf3=[]
    if GV[4]==1
        for ii in 1:length(sklist)-1
            infa1=Kdata_dt[sklist[ii]].tp-GV[5][1]
            infa2=Kdata_dt[sklist[ii]+1].tp-GV[5][1]
            push!(inf3,(infa1,infa2))
        end
    else
    for ii in 1:length(sklist)-1
        infa1=GV[5][2]-Kdata_dt[sklist[ii]].td
        infa2=GV[5][2]-Kdata_dt[sklist[ii]+1].td
        push!(inf3,(infa1,infa2))
    end
    end
    ntrsk=findmax(inf3)[2]
    list=collect(sklist[1]:sklist[ntrsk])
    lenlist1=length(list)
    tp=minimum(Kdata_dt[ii].tp for ii in list)
    td=maximum(Kdata_dt[ii].td for ii in list)
    tpr=maximum(Kdata_dt[ii].tp for ii in list)
    tdr=minimum(Kdata_dt[ii].td for ii in list)
    SKdata[GV[3]]=scommk(SKdata[GV[3]].ok,list,tp,td,tpr,tdr,SKdata[GV[3]].id)
    list=collect(sklist[ntrsk+1]:last(sklist))
    lenlist2=length(list)
    tp=minimum(Kdata_dt[ii].tp for ii in list)
    td=maximum(Kdata_dt[ii].td for ii in list)
    tpr=maximum(Kdata_dt[ii].tp for ii in list)
    tdr=minimum(Kdata_dt[ii].td for ii in list)
    idsk=last(SKdata).id+1
    push!(SKdata,scommk(SKdata[GV[3]].ok,list,tp,td,tpr,tdr,idsk))
    idne=findall(ne -> ne.n.id==1 && ne.t==SKdata[GV[3]].td,NE)
    if idne == Int64[]
        idnn=last(NE).id+1
        push!(NE,node(Fnode[1],SKdata[GV[3]].td,idnn))
        Dk[GV[3]]=NE[idnn]
    else
        Dk[GV[3]]=NE[idne][1]
    end
    idne=findall(ne -> ne.n.id==1 && ne.t==SKdata[idsk].td,NE)
    if idne != Int64[]
        push!(Dk,NE[idne][1])
    else
        idnn=last(NE).id+1
        push!(NE,node(Fnode[1],SKdata[idsk].td,idnn))
        push!(Dk,NE[idnn])
    end
    idne=findall(ne -> ne.n==SKdata[idsk].ok && ne.t==SKdata[idsk].tp,NE)
    if idne==Int64[]
        idnn=last(NE).id+1
        push!(NE,node(SKdata[idsk].ok,SKdata[idsk].tp,idnn))
        push!(Ok,NE[idnn])
    else
        push!(Ok,NE[idne][1])
    end
    idne=findall(ne->ne.n.id==1 && ne.t==SKdata[idsk].tp-ttm_dt[1,SKdata[idsk].ok.id],NE)
    if idne==Int64[]
        idnn=last(NE).id+1
        push!(NE,node(Fnode[1],SKdata[idsk].tp-ttm_dt[1,SKdata[idsk].ok.id],idnn))
    end
    if lenlist1>=2*lenlist2
        cond3=1
        idk=GV[3]
    elseif lenlist2>2*lenlist1
        cond3=1
        idk=last(SKdata).id
    end
    if cond3==1
        sklist=SKdata[idk].list
        ntrsk=cld(length(sklist),2)
        list=collect(sklist[1]:sklist[ntrsk])
        tp=minimum(Kdata_dt[ii].tp for ii in list)
        td=maximum(Kdata_dt[ii].td for ii in list)
        tpr=maximum(Kdata_dt[ii].tp for ii in list)
        tdr=minimum(Kdata_dt[ii].td for ii in list)
        SKdata[idk]=scommk(SKdata[idk].ok,list,tp,td,tpr,tdr,SKdata[idk].id)
        list=collect(sklist[ntrsk+1]:last(sklist))
        tp=minimum(Kdata_dt[ii].tp for ii in list)
        td=maximum(Kdata_dt[ii].td for ii in list)
        tpr=maximum(Kdata_dt[ii].tp for ii in list)
        tdr=minimum(Kdata_dt[ii].td for ii in list)
        idsk=last(SKdata).id+1
        push!(SKdata,scommk(SKdata[idk].ok,list,tp,td,tpr,tdr,idsk))
        idne=findall(ne -> ne.n.id==1 && ne.t==SKdata[idk].td,NE)
        if idne == Int64[]
            idnn=last(NE).id+1
            push!(NE,node(Fnode[1],SKdata[idk].td,idnn))
            Dk[idk]=NE[idnn]
        else
            Dk[idk]=NE[idne][1]
        end
        idne=findall(ne -> ne.n.id==1 && ne.t==SKdata[idsk].td,NE)
        if idne != Int64[]
            push!(Dk,NE[idne][1])
        else
            idnn=last(NE).id+1
            push!(NE,node(Fnode[1],SKdata[idsk].td,idnn))
            push!(Dk,NE[idnn])
        end
        idne=findall(ne -> ne.n==SKdata[idsk].ok && ne.t==SKdata[idsk].tp,NE)
        if idne==Int64[]
            idnn=last(NE).id+1
            push!(NE,node(SKdata[idsk].ok,SKdata[idsk].tp,idnn))
            push!(Ok,NE[idnn])
        else
            push!(Ok,NE[idne][1])
        end
        idne=findall(ne->ne.n.id==1 && ne.t==SKdata[idsk].tp-ttm_dt[1,SKdata[idsk].ok.id],NE)
        if idne==Int64[]
            idnn=last(NE).id+1
            push!(NE,node(Fnode[1],SKdata[idsk].tp-ttm_dt[1,SKdata[idsk].ok.id],idnn))
        end
    end
end
#-----------------------------------------------------------------------------
# partial expanded correction function
function partial_expanded_correction_function(PEV::Tuple{Bool,Array{Int64,1},Int64, Bool},NE,
    arc_ev,ttm_dt,Fnode2,initial_time)
    id=last(NE).id + 1
    for i in PEV[2]
        idne=findall(ne->ne.n==arc_ev[i].head.n&&ne.t==arc_ev[i].tail.t+ttm_dt[arc_ev[i].tail.n.id,arc_ev[i].head.n.id],NE)
        if idne==Int64[] && arc_ev[i].tail.t+ttm_dt[arc_ev[i].tail.n.id,arc_ev[i].head.n.id] <= Fnode2[arc_ev[i].head.n.id].bg + 30
            push!(NE,node(arc_ev[i].head.n,arc_ev[i].tail.t+ttm_dt[arc_ev[i].tail.n.id,arc_ev[i].head.n.id],id))
            id+=1
        end
        idne=findall(ne->ne.n==arc_ev[i].tail.n&&ne.t==arc_ev[i].head.t-ttm_dt[arc_ev[i].tail.n.id,arc_ev[i].head.n.id],NE)
        if idne==Int64[]&&arc_ev[i].head.t-ttm_dt[arc_ev[i].tail.n.id,arc_ev[i].head.n.id]>=initial_time && arc_ev[i].head.t-ttm_dt[arc_ev[i].tail.n.id,arc_ev[i].head.n.id] <= Fnode2[arc_ev[i].tail.n.id].bg +30
            push!(NE,node(arc_ev[i].tail.n,arc_ev[i].head.t-ttm_dt[arc_ev[i].tail.n.id,arc_ev[i].head.n.id],id))
            id+=1
        end
    end
end
function cost3(ij, ttm)
    if ij.tail.n.id!=ij.head.n.id
        cost=ttm[ij.tail.n.id,ij.head.n.id]
    else
        cost=0
    end
    return cost
end
#-------------------------------------------------------------------------------
function arcsused(y,arc_ev,N)
    idarcssol = findall(ij -> value(y[ij]) >= 0.5, arc_ev)
    cm = zeros(Int, N, N)
    for i in 1:N
        cm[i, i] = 1
    end
    for id in idarcssol
        i = arc_ev[id].tail.n.id
        j = arc_ev[id].head.n.id
        if cm[i, j] == 0
            cm[i, j] = 1
        end
    end
    return cm
end
#GRASPheurisitic
function savings(i,j,ttm)
    s=ttm[i,1]+ttm[1,j]-ttm[i,j]
    return s
end
function solution_cost(routes,ttm)
    cost = 0
    for r in routes
        cost+= ttm[1,first(r[2])]
        for i in 2:length(r[2])
            cost+=ttm[r[2][i-1],r[2][i]]
        end
        cost+=ttm[last(r[2]),1]
    end 
    return cost   
end
function grasp_construction(SKdata, W,Kdata,Fnode2,ttm_dt,ttm,type="SK")
    kdata2 = [k for k in SKdata]
    if W == 0
        for i in eachindex(kdata2)
            if length(kdata2[i].list)==1
                continue
            end
            sklist=kdata2[i].list
            ntrsk=cld(length(sklist),2)
            list=collect(sklist[1]:sklist[ntrsk])
            tp=minimum(Kdata[ii].tp for ii in list)
            td=maximum(Kdata[ii].td for ii in list)
            tpr=maximum(Kdata[ii].tp for ii in list)
            tdr=minimum(Kdata[ii].td for ii in list)
            kdata2[i]=scommk(kdata2[i].ok,list,tp,td,tpr,tdr,kdata2[i].id)
            list=collect(sklist[ntrsk+1]:last(sklist))
            tp=minimum(Kdata[ii].tp for ii in list)
            td=maximum(Kdata[ii].td for ii in list)
            tpr=maximum(Kdata[ii].tp for ii in list)
            tdr=minimum(Kdata[ii].td for ii in list)
            idsk=last(kdata2).id+1
            push!(kdata2,scommk(kdata2[i].ok,list,tp,td,tpr,tdr,idsk))
        end
    end
    routes = []
    if type == "SK"
        for k in kdata2
            tp = k.tpr
            td = tp
            sp = minimum([Fnode2[k.ok.id].bg+30,k.tdr-ttm_dt[k.ok.id,1]])
            sd = k.tdr - (td + ttm_dt[k.ok.id,1])
            slack = minimum([sp,sd])
            push!(routes,[[k.id],[k.ok.id],tp,td,slack])
        end
    else
        for k in Kdata_dt
            tp = k.tp
            td = tp
            sp = minimum([Fnode2[k.ok.id].bg+30,k.td-ttm_dt[k.ok.id,1]])
            sd = k.td - (td + ttm_dt[k.ok.id,1])
            slack = minimum([sp,sd])
            push!(routes,[[k.id],[k.ok.id],tp,td,slack])
        end
    end
    c = 0
    cond = true
    while cond
        if c >= 200
            cond = false
            continue
        end 
        s = []
        for r1 in eachindex(routes)
            for r2 in eachindex(routes)
                if r2 != r1
                    i = last(routes[r1][2])
                    j = first(routes[r2][2])
                    vs = savings(i,j,ttm)
                    if vs < 0
                        continue
                    end
                    tp2 = max(routes[r1][4]+ttm_dt[i,j],routes[r2][3])
                    dtp = tp2 - routes[r2][3]
                    if dtp > routes[r2][5]
                        continue
                    end
                    td2 = routes[r2][4]+dtp
                    ii = last(routes[r2][2])
                    dtd = td2 + ttm_dt[ii,1] - (routes[r1][4] + ttm_dt[i,1])
                    if dtd > routes[r1][5]
                        continue
                    end
                    slack = min(routes[r2][5]-dtp, routes[r1][5]-dtd)
                    push!(s,[r1,r2,vs,slack,td2])
                end
            end
        end
        df_s = DataFrame(r1=[],r2=[],vs=[],s=[],td=[])
        for r in s
            push!(df_s,r)
        end
        sort!(df_s,[:vs,:s],rev=true)
        if size(df_s,1) == 0
            cond = false
            continue
        end
        if size(df_s,1) < 5
            lr = size(df_s,1)
        else
            lr = 5
        end
        row = df_s[rand(1:lr),:]
        r1 = routes[row.r1]
        r2 = routes[row.r2]
        new_route = [vcat(r1[1],r2[1]),vcat(r1[2],r2[2]),r1[3],row.td,row.s]
        routes[row.r1]=new_route
        deleteat!(routes,row.r2)
        c+=1
    end
    value = solution_cost(routes,ttm)
    return routes, value
end