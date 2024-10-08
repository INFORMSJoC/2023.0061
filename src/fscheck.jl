function fscheck()
    idarcssol=findall(ij -> value(y[ij])>0,arc_ev)
    arcsol=arc_ev[idarcssol]
    arcsbrute=arcsol[findall(ij ->  ij.tail.n.id==1 && ij.head.n.id!=1,arcsol)]
    arcs_used = []
    routes = []
    routes_ = []
    for r in eachindex(arcsbrute)
        push!(routes_,arcsbrute[r])
        dec=0
        while dec!=1
            idarcsol_=findall(ij -> ij.tail==last(routes_).head && ij.tail.n.id!=1 && ij ∉ arcs_used ,arcsol)
            if length(idarcsol_)>1
                push!(arcs_used,arcsol[idarcsol_[1]])
            end
            if idarcsol_!=Int64[]
                push!(routes_,arcsol[idarcsol_[1]])
            else
                dec=1
            end
        end
        push!(routes,routes_)
        routes_=[]
    end

    seq = []
    
    for j in eachindex(routes)
        d_arcs=[r for r in routes[j] if r.tail.n != r.head.n]
        for i in eachindex(d_arcs)
            if d_arcs[i].tail.n.id==1
                n = d_arcs[i].tail.n.id
                l = 0
                u = final_time
                push!(seq,[n,l,u,j,1])
            elseif d_arcs[i].head.n.id==1
                n = d_arcs[i].tail.n.id
                k_set=[ij.k for ij in arc_ek[findall(ij -> value(x[ij])>0 && ij.id == d_arcs[i].id, arc_ek)] if SKdata[ij.k].ok == ij.tail.n]
                if length(k_set)>0
                    l = maximum([k.tpr] for k in SKdata[k_set])[1]            
                    f_node = Fnode2[d_arcs[i].tail.n.id]
                    u = f_node.bg+30
                    push!(seq,[n,l,u,j,1])
                end
                n = d_arcs[i].head.n.id
                l = 0
                k_set = [ij.k for ij in arc_ek[findall(ij -> value(x[ij])>0 && ij.id == d_arcs[i].id, arc_ek)]]
                u = minimum([k.tdr] for k in SKdata[k_set])[1]
                push!(seq,[n,l,u,j,1])
            else
                n = d_arcs[i].tail.n.id
                k_set=[ij.k for ij in arc_ek[findall(ij -> value(x[ij])>0 && ij.id == d_arcs[i].id, arc_ek)] if SKdata[ij.k].ok == ij.tail.n]
                if length(k_set)>0
                    l = maximum([k.tpr] for k in SKdata[k_set])[1]
                    f_node = Fnode2[d_arcs[i].tail.n.id]
                    u = f_node.bg+30
                    push!(seq,[n,l,u,j,1])
                end
            end
        end
    end
    p_id=2
    for i in 2:N
        sor_p = sort(seq[findall(j ->j[1]==i,seq)])
        for j in sor_p
            id_seq = findall(s -> s == j,seq)
            seq[id_seq[1]][5]=p_id
            p_id+=1
        end
    end
    P = p_id - 1
    V = [v_i(1,1)]
    for p in 2:P
        info = seq[findall(i ->i[5] == p, seq)[1]]
        push!(V,v_i(info[5],info[1]))
    end
    V
    I=[]
    for c in Fnode2
        push!(I,findall(i -> i.cc == c.id,V))
    end
    p_seq = []
    for i in 1:(length(seq)-1)
        push!(p_seq,[seq[i][5],seq[i+1][5]])
    end
    TMAX = maximum(Fnode2[i].tmax_g+Fnode2[i].delta_g for i in 1:N)
    lp=Model(CPLEX.Optimizer);
    set_optimizer_attribute(model, "CPX_PARAM_TILIM", model_timelimit);
    @variable(lp,u[i in 1:P] >= 0);
    @variable(lp,d[i in 1:P] >=0);
    @variable(lp,f[i in 1:P] >=0);
    @objective(lp,Min,sum(d[i] for i in 2:P));
    @constraints(lp,begin
        [i in p_seq;i[1]!=1 && i[2]!=1],u[i[2]]>=u[i[1]] + ttm[V[i[1]].cc,V[i[2]].cc]
        [c in Fnode2; c.id!=1], u[I[c.id][1]] >= c.ag
        [c in Fnode2, k in I[c.id]; c.id!=1 && length(I[c.id])>2 && k>I[c.id][1]], u[k] >= u[k-1]
        [c in Fnode2; c.id!=1], u[last(I[c.id])] >= c.bg
        [c in Fnode2; c.id!=1], u[last(I[c.id])] <= c.bg + 30
        [i in p_seq;i[1]!=1 && i[2]==1], TMAX - f[i[1]] >= ttm[V[i[1]].cc,1]
        [i in p_seq;i[1]!=1 && i[2]!=1], f[i[2]]-f[i[1]] >= u[i[2]]-u[i[1]]
        [c in Fnode2, i in 2:P; c.id!=1 && i == I[c.id][1] ], f[i] - u[i] + c.ag >= TMAX - (Fnode2[V[i].cc].tmax_g + Fnode2[V[i].cc].delta_g)
        [c in Fnode2, i in I[c.id]; c.id!=1 && length(I[c.id])>=2 && i>I[c.id][1]], f[i] - u[i] + u[i-1] >= TMAX - (Fnode2[V[i].cc].tmax_g + Fnode2[V[i].cc].delta_g)
        [i in p_seq;i[1]==1 && i[2]!=1], d[i[2]] >= TMAX - f[i[2]] + ttm[1,V[i[2]].cc]
    end)
    optimize!(lp)
    if has_values(lp)
        v_times = [value(u[p]) for p in 1:P]
        f_times = [value(f[p]) for p in 1:P]
        for i in eachindex(seq)
            t = v_times[seq[i][5]]
            f = f_times[seq[i][5]]
            append!(seq[i],t,t+f,[],[])
        end
    end
    return has_values(lp),seq
end

function solution(y,arc_ev,x,arc_ek,SKdata,final_time)
    idarcssol=findall(ij -> value(y[ij])>0,arc_ev)
    arcsol=arc_ev[idarcssol]
    arcsbrute=arcsol[findall(ij ->  ij.tail.n.id==1 && ij.head.n.id!=1,arcsol)]
    arcs_used = []
    routes = []
    routes_ = []
    for r in eachindex(arcsbrute)
        push!(routes_,arcsbrute[r])
        dec=0
        while dec!=1
            idarcsol_=findall(ij -> ij.tail==last(routes_).head && ij.tail.n.id!=1 && ij ∉ arcs_used ,arcsol)
            if length(idarcsol_)>1
                push!(arcs_used,arcsol[idarcsol_[1]])
            end
            if idarcsol_!=Int64[]
                push!(routes_,arcsol[idarcsol_[1]])
            else
                dec=1
            end
        end
        push!(routes,routes_)
        routes_=[]
    end

    seq = []
    
    for j in eachindex(routes)
        d_arcs=[r for r in routes[j] if r.tail.n != r.head.n]
        for i in eachindex(d_arcs)
            if d_arcs[i].tail.n.id==1
                n = d_arcs[i].tail.n.id
                t = d_arcs[i].tail.t
                l = 0
                u = final_time
                ks =  []
                k_set =  []
                t_lab =  0
                push!(seq,[n,l,u,j,1,t,t_lab,ks,k_set])
            elseif d_arcs[i].head.n.id==1
                n = d_arcs[i].tail.n.id
                t = d_arcs[i].tail.t
                k_set=[ij.k for ij in arc_ek[findall(ij -> value(x[ij])>0 && ij.id == d_arcs[i].id, arc_ek)] if SKdata[ij.k].ok == ij.tail.n]
                ks = [ij.k for ij in arc_ek[findall(ij -> value(x[ij])>0 && ij.id == d_arcs[i].id, arc_ek)]]
                if length(k_set)>0
                    l = maximum([k.tpr] for k in SKdata[k_set])[1]
                    t_lab = minimum([k.tdr] for k in SKdata[ks])[1]          
                    f_node = Fnode2[d_arcs[i].tail.n.id]
                    u = f_node.bg+30
                    push!(seq,[n,l,u,j,1,t,t_lab,ks,k_set])
                end
                n = d_arcs[i].head.n.id
                t = d_arcs[i].head.t
                l = 0
                k_set = [ij.k for ij in arc_ek[findall(ij -> value(x[ij])>0 && ij.id == d_arcs[i].id, arc_ek)]]
                ks = k_set
                u = minimum([k.tdr] for k in SKdata[k_set])[1]
                t_lab = u
                push!(seq,[n,l,u,j,1,t,t_lab,ks,k_set])
            else
                n = d_arcs[i].tail.n.id
                t = d_arcs[i].tail.t
                k_set=[ij.k for ij in arc_ek[findall(ij -> value(x[ij])>0 && ij.id == d_arcs[i].id, arc_ek)] if SKdata[ij.k].ok == ij.tail.n]
                ks = [ij.k for ij in arc_ek[findall(ij -> value(x[ij])>0 && ij.id == d_arcs[i].id, arc_ek)]]
                if length(k_set)>0
                    l = maximum([k.tpr] for k in SKdata[k_set])[1]
                    t_lab = minimum([k.tdr] for k in SKdata[ks])[1] 
                    f_node = Fnode2[d_arcs[i].tail.n.id]
                    u = f_node.bg+30
                    push!(seq,[n,l,u,j,1,t,t_lab,ks,k_set])
                end
            end
        end
    end
    p_id=2
    for i in 2:N
        sor_p = sort(seq[findall(j ->j[1]==i,seq)])
        for j in sor_p
            id_seq = findall(s -> s == j,seq)
            seq[id_seq[1]][5]=p_id
            p_id+=1
        end
    end
    return seq
end

function solutionup(y,arc_ev,x,arc_ek,SKdata)
    idarcssol=findall(ij -> value(y[ij])>0,arc_ev)
    arcsol=arc_ev[idarcssol]
    arcsbrute=arcsol[findall(ij ->  ij.tail.n.id==1 && ij.head.n.id!=1,arcsol)]
    arcs_used = []
    routes = []
    routes_ = []
    for r in eachindex(arcsbrute)
        push!(routes_,arcsbrute[r])
        dec=0
        while dec!=1
            idarcsol_=findall(ij -> ij.tail==last(routes_).head && ij.tail.n.id!=1 && ij ∉ arcs_used ,arcsol)
            if length(idarcsol_)>1
                push!(arcs_used,arcsol[idarcsol_[1]])
            end
            if idarcsol_!=Int64[]
                push!(routes_,arcsol[idarcsol_[1]])
            else
                dec=1
            end
        end
        push!(routes,routes_)
        routes_=[]
    end

    seq = []
    
    for j in eachindex(routes)
        d_arcs=[r for r in routes[j] if r.tail.n != r.head.n]
        for i in eachindex(d_arcs)
            if d_arcs[i].tail.n.id==1
                n = d_arcs[i].tail.n.id
                t = d_arcs[i].tail.t
                l = 0
                u = final_time
                ks =  []
                k_set =  []
                t_lab =  0
                push!(seq,[n,l,u,j,1,t,t_lab,ks,k_set])
            elseif d_arcs[i].head.n.id==1
                n = d_arcs[i].tail.n.id
                t = d_arcs[i].tail.t
                k_set=[ij.k for ij in arc_ek[findall(ij -> value(x[ij])>0 && ij.id == d_arcs[i].id, arc_ek)] if SKdata[ij.k].ok == ij.tail.n]
                ks = [ij.k for ij in arc_ek[findall(ij -> value(x[ij])>0 && ij.id == d_arcs[i].id, arc_ek)]]
                if length(k_set)>0
                    l = maximum([ k.tp] for k in SKdata[k_set])[1]
                    t_lab = minimum([k.td] for k in SKdata[ks])[1]          
                    f_node = Fnode2[d_arcs[i].tail.n.id]
                    u = f_node.bg+30
                    push!(seq,[n,l,u,j,1,t,t_lab,ks,k_set])
                end
                n = d_arcs[i].head.n.id
                t = d_arcs[i].head.t
                l = 0
                k_set = [ij.k for ij in arc_ek[findall(ij -> value(x[ij])>0 && ij.id == d_arcs[i].id, arc_ek)]]
                ks = k_set
                u = minimum([k.td] for k in SKdata[k_set])[1]
                t_lab = u
                push!(seq,[n,l,u,j,1,t,t_lab,ks,k_set])
            else
                n = d_arcs[i].tail.n.id
                t = d_arcs[i].tail.t
                k_set=[ij.k for ij in arc_ek[findall(ij -> value(x[ij])>0 && ij.id == d_arcs[i].id, arc_ek)] if SKdata[ij.k].ok == ij.tail.n]
                ks = [ij.k for ij in arc_ek[findall(ij -> value(x[ij])>0 && ij.id == d_arcs[i].id, arc_ek)]]
                if length(k_set)>0
                    l = maximum([ k.tp] for k in SKdata[k_set])[1]
                    t_lab = minimum([k.td] for k in SKdata[ks])[1] 
                    f_node = Fnode2[d_arcs[i].tail.n.id]
                    u = f_node.bg+30
                    push!(seq,[n,l,u,j,1,t,t_lab,ks,k_set])
                end
            end
        end
    end
    p_id=2
    for i in 2:N
        sor_p = sort(seq[findall(j ->j[1]==i,seq)])
        for j in sor_p
            id_seq = findall(s -> s == j,seq)
            seq[id_seq[1]][5]=p_id
            p_id+=1
        end
    end
    return seq
end