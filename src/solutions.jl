# solution values function
function solutionvalues(y,x)
    R=cld(sum(value(y[ij]) for ij in arc_ev if ij.tail.n.id==1 && ij.head.n.id!=1 ),1)
    m=sum(value(y[ij]) for ij in arc_ev if ij.head.id==LT.id)
    idarcssol=findall(ij -> value(y[ij])>0.5,arc_ev)
    arcsol=arc_ev[idarcssol]
    arcsbrute=arcsol[findall(ij ->  ij.tail.n.id==1 && ij.head.n.id!=1,arcsol)]
    routes=[]
    routes_=[]
    printcondition=false
    R2=R
    if Int64(R)!=length(arcsbrute)
        printcondition=true
        R=length(arcsbrute)
    end
    for r in eachindex(arcsbrute)
        push!(routes_,arcsbrute[r])
        dec=0
        while dec!=1
            idarcsol_=findall(ij -> ij.tail==last(routes_).head && ij.tail.n.id!=1,arcsol)
            if idarcsol_!=Int64[]
                push!(routes_,arcsol[idarcsol_[1]])
            else
                dec=1
            end
        end
        push!(routes,routes_)
        routes_=[]
    end
    x_routes=[]
    k_solution=[]
    k_condition=[]
    pickuptimes=[]
    for k in eachindex(SKdata)
        pickupid=findall(ijk->value(x[ijk])!=0 && ijk.tail.n==SKdata[k].ok && ijk.head.n!=SKdata[k].ok,arcsk[k])
        arcvpickup=arc_ev[arcsk[k][pickupid[1]].id]
        deliveryid=findall(ijk->value(x[ijk])!=0 && ijk.head.n==Dk[k].n && ijk.tail.n!=Dk[k].n,arcsk[k])
        arcvdelivery=arc_ev[arcsk[k][pickupid[1]].id]
        ktp=arcvpickup.tail.t
        ktd=arcvdelivery.head.t
        push!(pickuptimes,ktp)
        krouteid=[]
        kconditional=true
        for r=1:length(arcsbrute)
            if arcvpickup ∈ routes[r]
                push!(krouteid,r)
                push!(x_routes,[k r])
                kconditional=false
            end
        end
        if kconditional
            printcondition=true
            push!(krouteid,0)
            push!(x_routes,[0 0])
        end
        kpickup_cond=SKdata[k].tpr<=ktp
        kdelivery_cond=SKdata[k].tdr>=ktd
        kcondition=kpickup_cond&&kdelivery_cond
        push!(k_condition,kcondition)
        push!(k_solution,[krouteid[1],kcondition,SKdata[k].tpr,ktp,SKdata[k].tdr,ktd])
    end
        msn="\n"
    if printcondition
        msn="The Routes tracking was not possible\n"
        println("\n")
        println("The Routes tracking was not possible\n")
    end
    pickupsolution=[]
    for i=2:N
        arcspickupsolution=arcsol[findall(ij->ij.tail.n.id!=i&&ij.head.n.id==i,arcsol)]
        timepickup=[ii.head.t for ii in arcspickupsolution]
        push!(pickupsolution,[i,length(arcspickupsolution),timepickup])
    end
    return R,m,routes,k_solution,k_condition,x_routes,pickuptimes,msn,pickupsolution,R2
end
#------------------------------------------------------------------------------
#last solution function
function lastsolution(model,x,y,NE,SKdata,arc_ev,arc_ek,arcsk)
    lastmodel=model
    lastx=x
    lasty=y
    lastne=NE
    lastskdata=SKdata
    lastarc_ev=arc_ev
    lastarc_ek=arc_ek
    lastarcsk=arcsk
    return lastmodel,lastx,lasty,lastne,lastskdata,lastarc_ev,lastarc_ek,lastarcsk
end
#-------------------------------------------------------------------------------
# function solutionvalues2(y,arc_ev,LT,K,x,Kdata,arcsk,Dk)
#     R=cld(sum(value(y[ij]) for ij in arc_ev if ij.tail.n.id==1 && ij.head.n.id!=1 ),1)
#     m=sum(value(y[ij]) for ij in arc_ev if ij.head.id==LT.id)
#     idarcssol=findall(ij -> value(y[ij])>0.5,arc_ev)
#     arcsol=arc_ev[idarcssol]
#     arcsbrute=arcsol[findall(ij ->  ij.tail.n.id==1 && ij.head.n.id!=1,arcsol)]
#     routes=[]
#     routes_=[]
#     printcondition=false
#     for r=1:Int64(R)
#         push!(routes_,arcsbrute[r])
#         dec=0
#         while dec!=1
#             idarcsol_=findall(ij -> ij.tail==last(routes_).head && ij.tail.n.id!=1,arcsol)
#             if idarcsol_!=Int64[]
#                 push!(routes_,arcsol[idarcsol_[1]])
#             else
#                 dec=1
#             end
#         end
#         push!(routes,routes_)
#         routes_=[]
#     end
#     x_routes=[]
#     k_solution=[]
#     k_condition=[]
#     pickuptimes=[]
#     for k=1:K
#         pickupid=findall(ijk->value(x[ijk])!=0 && ijk.tail.n==Kdata[k].ok && ijk.head.n!=Kdata[k].ok,arcsk[k])
#         arcvpickup=arc_ev[arcsk[k][pickupid[1]].id]
#         deliveryid=findall(ijk->value(x[ijk])!=0 && ijk.head.n==Dk[k].n && ijk.tail.n!=Dk[k].n,arcsk[k])
#         arcvdelivery=arc_ev[arcsk[k][pickupid[1]].id]
#         ktp=arcvpickup.tail.t
#         ktd=arcvdelivery.head.t
#         push!(pickuptimes,ktp)
#         krouteid=[]
#         kconditional=true
#         for r=1:Int64(R)
#             if arcvpickup ∈ routes[r]
#                 push!(krouteid,r)
#                 push!(x_routes,[k r])
#                 kconditional=false
#             end
#         end
#         if kconditional
#             printcondition=true
#             push!(krouteid,0)
#             push!(x_routes,[0 0])
#         end
#         kpickup_cond=Kdata[k].tp<=ktp
#         kdelivery_cond=Kdata[k].td>=ktd
#         kcondition=kpickup_cond&&kdelivery_cond
#         push!(k_condition,kcondition)
#         push!(k_solution,[krouteid[1],kcondition,Kdata[k].tp,ktp,Kdata[k].td,ktd])
#     end
#     if printcondition
#         println("\n")
#         println("The Routes tracking was not possible\n")
#     end
#     return R,m,routes,k_solution,k_condition,x_routes,pickuptimes
# end
#-------------------------------------------------------------------------------
function writesol(modelup,yup,arc_evup,LTup,Kup,xup,Kdataup,arcskup,Dkup)
    R=cld(sum(value(yup[ij]) for ij in arc_evup if ij.tail.n.id==1 && ij.head.n.id!=1 ),1)
    m=sum(value(yup[ij]) for ij in arc_evup if ij.head.id==LTup.id)
    idarcssol=findall(ij -> value(yup[ij])>0.5,arc_evup)
    arcsol=arc_evup[idarcssol]
    arcsbrute=arcsol[findall(ij ->  ij.tail.n.id==1 && ij.head.n.id!=1,arcsol)]
    routes=[]
    routes_=[]
    printcondition=false
    R2=R
    if Int64(R)!=length(arcsbrute)
        printcondition=true
        R=length(arcsbrute)
    end
    for r=1:Int64(R)
        push!(routes_,arcsbrute[r])
        dec=0
        while dec!=1
            idarcsol_=findall(ij -> ij.tail==last(routes_).head && ij.tail.n.id!=1,arcsol)
            if idarcsol_!=Int64[]
                push!(routes_,arcsol[idarcsol_[1]])
            else
                dec=1
            end
        end
        push!(routes,routes_)
        routes_=[]
    end
    x_routes=[]
    k_solution=[]
    k_condition=[]
    pickuptimes=[]
    for k=1:Kup
        pickupid=findall(ijk->value(xup[ijk])!=0 && ijk.tail.n==Kdataup[k].ok && ijk.head.n!=Kdataup[k].ok,arcskup[k])
        arcvpickup=arc_evup[arcskup[k][pickupid[1]].id]
        deliveryid=findall(ijk->value(xup[ijk])!=0 && ijk.head.n==Dkup[k].n && ijk.tail.n!=Dkup[k].n,arcskup[k])
        arcvdelivery=arc_evup[arcskup[k][deliveryid[1]].id]
        ktp=arcvpickup.tail.t
        ktd=arcvdelivery.head.t
        push!(pickuptimes,ktp)
        krouteid=[]
        kconditional=true
        for r=1:Int64(R)
            if arcvpickup ∈ routes[r]
                push!(krouteid,r)
                push!(x_routes,[k r])
                kconditional=false
            end
        end
        if kconditional
            printcondition=true
            push!(krouteid,0)
            push!(x_routes,[0 0])
        end
        kpickup_cond=Kdataup[k].tp<=ktp
        kdelivery_cond=Kdataup[k].td>=ktd
        kcondition=kpickup_cond&&kdelivery_cond
        push!(k_condition,kcondition)
        push!(k_solution,[krouteid[1],kcondition,Kdataup[k].tp,ktp,Kdataup[k].td,ktd])
    end
        msn="\n"
        if printcondition
            msn="The Routes tracking was not possible\n"
            println("\n")
            println("The Routes tracking was not possible\n")
        end
    pickupsolution=[]
    for i=2:N
        arcspickupsolution=arcsol[findall(ij->ij.tail.n.id!=i&&ij.head.n.id==i,arcsol)]
        timepickup=[ii.head.t for ii in arcspickupsolution]
        push!(pickupsolution,[i,length(arcspickupsolution),timepickup])
    end
    open("$dir/results/$instance/iterations/dt$(dt)_dk$(dk)_$itnumber-solup.csv", "w") do f
        write(f,"Upperbound Solution for instance $instance iteration $itnumber\n")
        write(f,"SCC;$(N-1);Number of TR;$(length(Kdataup))\n")
        write(f,"Maximum commodity time; $commodity_time_limit minutes ;dk;$dk;;MIP solution time limit;$model_timelimit s\n")
        write(f,"SCC;Number of pickups;pickups time\n")
        for i in pickupsolution
            write(f,"$(i[1]);$(i[2]);$(i[3])\n")
        end
        write(f,"Routes=;$R2\n")
        write(f,"$msn\n")
        for r=1:Int64(R)
            write(f,"route $r\n")
            nodeswrite=transpose(findall(i ->i[2]==r,x_routes))
            write(f,"Commodities=,$nodeswrite\n")
            write(f,"node,time\n")
            write(f,"$(routes[r][1].tail.n.id),$(routes[r][1].tail.t)\n")
            for arcs_ in routes[r]
                writek=findall(i->i==arcs_.head.t,pickuptimes)
                write(f,"$(arcs_.head.n.id),$(arcs_.head.t),$writek\n")
            end
            write(f,"duration=,$(-(routes[r][1].tail.t)+(last(routes[r]).head.t))\n")
            write(f,"\n")
        end
    end
    return pickupsolution
end
