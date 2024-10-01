function upperboundfunction(N::Int64,SKdata::Array{Any,1},ttm::Array{Int64,2},Fnode2,initial_time,final_time,cm,dt,Fnode,
    K,model_timelimit)
    ttm_dt=zeros(Int,N,N)
    for i=1:N
        for j=1:N
            ttm_dt[i,j]=Int64(cld(ttm[i,j],dt)*dt)
        end
    end
    initial_time_dt=Int64(cld(initial_time,dt)*dt)
    final_time_dt=Int64(fld(final_time,dt)*dt)
    TH=final_time_dt-initial_time_dt
    T=Int64(TH/dt)
    Kdata_dt2 = []
    for i in SKdata
        push!(Kdata_dt2,commk(i.ok,Int64(cld(i.tpr, dt) * dt),Int64(fld(i.tdr, dt) * dt),i.id,))
    end
    NE2 =
        [node(Fnode[1], initial_time_dt + Int64((i - 1) * dt), i) for i = 1:T+1]
    id = last(NE2).id
    L02 = NE2[1]
    LT2 = NE2[id]
    for j = 2:N
        for t = 1:T+1
            if (initial_time_dt + dt * t) >=
               (initial_time_dt + ttm_dt[L02.n.id, Fnode[j].id]) &&
               (final_time_dt - dt * t) >=
               (initial_time_dt + ttm_dt[Fnode[j].id, LT2.n.id]) &&
               initial_time_dt + Int64((t - 1) * dt) <= (Fnode2[j].bg + 30)
                id += 1
                push!(
                    NE2,
                    node(Fnode[j], initial_time_dt + Int64((t - 1) * dt), id),
                )
            end
        end
    end
    arc_ev2 = [arc(NE2[i], NE2[j], ttm_dt[NE2[i].n.id, NE2[j].n.id], 1) for  i  in eachindex(NE2) for j in eachindex(NE2) if NE2[i].n != NE2[j].n &&
        NE2[j].t == NE2[i].t + ttm_dt[NE2[i].n.id, NE2[j].n.id] && cm[NE2[i].n.id, NE2[j].n.id] == 1]
    for i  in eachindex(arc_ev2)
        arc_ev2[i].id = i
    end
    id = last(arc_ev2).id
    for i in eachindex(NE2)
        for j  in eachindex(NE2)
            if NE2[i].n == NE2[j].n && NE2[j].t == NE2[i].t + dt
                id += 1
                push!(arc_ev2,arc(NE2[i], NE2[j], ttm_dt[NE2[i].n.id, NE2[j].n.id], id))
            end
        end
    end
    id = last(NE2).id
    for k in Kdata_dt2
        idne = findall(i -> (i.n == k.ok && i.t == k.tp), NE2)
        if idne == Int64[]
            id += 1
            push!(NE2, node(k.ok, k.tp, id))
        end
    end
    Ok2 = [node(Fnode[1], k, k) for k in eachindex(Kdata_dt2)]
    Dk2 = [node(Fnode[1], k, k) for k in eachindex(Kdata_dt2)]
    cont = 1
    for k in Kdata_dt2
        idne = findall(i -> i.n == k.ok && i.t == k.tp, NE2)
        Ok2[cont] = NE2[idne[1]]
        idne = findall(i -> i.n.id == 1 && i.t == k.td, NE2)
        Dk2[cont] = NE2[idne[1]]
        cont += 1
    end
    id = last(arc_ev2).id
    arc_ek2 = []
    for k in eachindex(Kdata_dt2)
        idarc = findall(ij -> ij.tail == Ok2[k], arc_ev2)
        if idarc == Int64[]
            sel = findall(i -> i.n.id == Ok2[k].n.id && i.id != Ok2[k].id, NE2)
            tm = minimum(NE2[i].t for i in sel)
            nh = findall(i -> i.n.id == Ok2[k].n.id && i.t == tm, NE2)
            id += 1
            push!(arc_ek2, arck(Ok2[k], NE2[nh[1]], k, id))
        else
            for i in idarc
                push!(
                    arc_ek2,
                    arck(arc_ev2[i].tail, arc_ev2[i].head, k, arc_ev2[i].id),
                )
            end
        end
    end
    selidarck = [findall(a -> a.k == k, arc_ek2) for k in eachindex(Kdata_dt2)]
    mod = false
    for k in eachindex(Kdata_dt2)
        dec = 0
        flagged = falses(last(arc_ev2).id + K)
        arc_ev_selid = findall(
            ij -> ij.tail.t >= Ok2[k].t && ij.head.t <= Dk2[k].t,
            arc_ev2,
        )
        arc_ev_sel = arc_ev2[arc_ev_selid]
        while dec != 1
            arcsk = [arc_ek2[i] for i in selidarck[k]]
            idarcsk = [arc_ek2[i].id for i in selidarck[k]]
            mod = false
            for ij in arcsk
                if flagged[ij.id]
                    continue
                end
                flagged[ij.id] = true
                if ij.head.n.id != 1
                    selidarcv = findall(j -> j.tail == ij.head, arc_ev_sel)
                    arcsv = arc_ev_sel[selidarcv]
                    for kl in arcsv
                        if kl.id ∉ idarcsk
                            if Dk2[k].t - kl.head.t >=
                               ttm_dt[kl.head.n.id, Dk2[k].n.id]
                                newarc = arck(kl.tail, kl.head, k, kl.id)
                                push!(arc_ek2, newarc)
                                push!(selidarck[k], length(arc_ek2))
                                push!(idarcsk, newarc.id)
                                mod = true
                            end
                        end
                    end
                else
                    selidarcv = findall(
                        j -> j.tail == ij.head && j.head.n.id == 1,
                        arc_ev_sel,
                    )
                    arcsv = arc_ev_sel[selidarcv]
                    for kl in arcsv
                        if kl.id ∉ idarcsk
                            if Dk2[k].t >= kl.head.t
                                newarc = arck(kl.tail, kl.head, k, kl.id)
                                push!(arc_ek2, newarc)
                                push!(selidarck[k], length(arc_ek2))
                                push!(idarcsk, newarc.id)
                                mod = true
                            end
                        end
                    end
                end
            end
            if !mod
                dec = 1
            end
        end
    end
    idarcsktail =
        [unique([ijk.tail.id for ijk in arc_ek2 if ijk.k == k]) for k in eachindex(Kdata_dt2)]
    idarcskhead =
        [unique([ijk.head.id for ijk in arc_ek2 if ijk.k == k]) for k in eachindex(Kdata_dt2)]
    idnek = [sort(unique(union(idarcsktail[k], idarcskhead[k]))) for k in eachindex(Kdata_dt2)]
    arcsk2 = []
    NEK2 = []
    selidarck = [findall(a -> a.k == k, arc_ek2) for k in eachindex(Kdata_dt2)]
    for k in eachindex(Kdata_dt2)
        push!(arcsk2, [arc_ek2[i] for i in selidarck[k]])
        push!(NEK2, [NE2[i] for i in idnek[k]])
    end
    idarcsvtail = unique([ij.tail.id for ij in arc_ev2])
    idarcsvhead = unique([ij.head.id for ij in arc_ev2])
    idnev = sort(unique(union(idarcsvtail, idarcsvhead)))
    NEV2 = [NE2[i] for i in idnev]
    # open("$dir/results/$instance/solution_log_dt$(dt)_dk$(dk)_$instance.txt","a") do f
	# write(f,"Upper bound network created it $itnumber $(Dates.value(Dates.now()-cpu_initial_time)/1000)\n")
    # end
    model2 = Model(optimizer_with_attributes(CPLEX.Optimizer))
    set_optimizer_attribute(model2, "CPX_PARAM_TILIM", model_timelimit)
    set_optimizer_attribute(model2, "CPXPARAM_Emphasis_MIP", 5)
    @variable(model2, y2[ij in arc_ev2] >= 0, Int)
    @variable(model2, x2[ijk in arc_ek2],Bin)
    @objective(model2, Min, sum(cost3(ij,ttm) * y2[ij] for ij in arc_ev2))
    # @constraint(
    #     model2,
    #     sum(y2[ij] for ij in arc_ev2 if ij.tail.id == L02.id) <= m
    # )
    for i in NEV2
        if i.id != L02.id && i.id != LT2.id
            @constraint(
                model2,
                sum(y2[ij] for ij in arc_ev2 if ij.tail.id == i.id) -
                sum(y2[ij] for ij in arc_ev2 if ij.head.id == i.id) == 0
            )
        end
    end
    for k in eachindex(Kdata_dt2)
        @constraint(
            model2,
            sum(x2[ijk] for ijk in arcsk2[k] if ijk.tail == Ok2[k]) == 1
        )
        @constraint(
            model2,
            sum(x2[ijk] for ijk in arcsk2[k] if ijk.head == Dk2[k]) == 1
        )
        for i in NEK2[k]
            if i.id != Ok2[k].id && i.id != Dk2[k].id
                @constraint(
                    model2,
                    sum(x2[ijk] for ijk in arcsk2[k] if ijk.tail.id == i.id) -
                    sum(x2[ijk] for ijk in arcsk2[k] if ijk.head.id == i.id) ==
                    0
                )
            end
        end
        for ijk in arcsk2[k]
            if ijk.tail.n != ijk.head.n
                @constraint(model2, x2[ijk] <= y2[arc_ev2[ijk.id]])
            end
        end
    end
    # open("$dir/results/$instance/solution_log_dt$(dt)_dk$(dk)_$instance.txt","a") do f
	# write(f,"Upper bound model created it $itnumber $(Dates.value(Dates.now()-cpu_initial_time)/1000)\n")
    # end
    optimize!(model2)
    # open("$dir/results/$instance/solution_log_dt$(dt)_dk$(dk)_$instance.txt","a") do f
	# write(f,"Upper bound model solved it $itnumber $(Dates.value(Dates.now()-cpu_initial_time)/1000)\n")
    # end
    if has_values(model2)
        uppb = objective_value(model2)
        cond = true
    else
        uppb = Inf
        cond = false
    end
    return uppb,cm,cond,model2,y2,arc_ev2,LT2,length(Kdata_dt2),x2,Kdata_dt2,arcsk2,Dk2,arc_ek2
end