using CSV, DataFrames, JuMP, Dates, DelimitedFiles,CPLEX,Graphs;
include("structures.jl");include("functions.jl");include("solutions.jl");include("upperbound.jl");
include("fscheck.jl")
function ddd_bstp(instance,dk,dt,path="inst/")
    if dt == 1
        dt2 = 5
    elseif dt == 5
        dt2 = 15
    end
    File_Name="$instance.csv";
    maximum_timilimit=21600;
    model_timelimit=3600;
    if dk==180
        commodity_time_limit=120
    elseif dk==60
        commodity_time_limit=240
    else
        commodity_time_limit=300
    end
    epsilon=0.001;
    println("\n")
    println("$instance dt=$dt dk=$dk commodity_time_limit=$(commodity_time_limit) minutes MIP solution time limit=$model_timelimit s\n")
    cpu_initial_time=Dates.now()
    data_read=data_read_function(File_Name,dk,path)
    global N=data_read[1];global Fnode=data_read[2];global ttm=data_read[3];global Kdata=data_read[4];global K=data_read[5];global m=data_read[8];
    global Fnode2=data_read[9];global dm=data_read[10]
    commodity_time_limit = maximum(fn.TMAX for fn in Fnode2);
    initial_time=data_read[6];final_time=data_read[7];
    data_dt=dt_approximation_function(N,ttm,dt,Kdata,initial_time,final_time)
    ttm_dt=data_dt[1];Kdata_dt=data_dt[2];initial_time_dt=data_dt[3];final_time_dt=data_dt[4];T=data_dt[5];
    imputs=inputs_fuction(N,Kdata_dt,ttm_dt,Fnode2,Fnode,initial_time_dt,initial_time)
    global NE=imputs[1];global SKdata=imputs[2];L0=imputs[3];LT=imputs[4];Ok=imputs[5];Dk=imputs[6]
    L03=imputs[7];
    Fnode2[1].bg = LT.t;
    initial_grouping=length(SKdata);
    kdata2 = imputs[2]
    #----------------------------------------------------------------------------
    grouping_log=[]
    partial_expanded_log=[]
    itnumber=0
    cont_k=1
    cont_n=0
    cont_upb=1
    factibility_condition=0
    lowerbound=-Inf
    llb=-Inf
    global upperbound=Inf
    timelimit_reached=false
    modelnovalues=false
    skdata2=[]
    cm=zeros(Int,N,N)
    upp=Inf
    condupp=false
    upp2=[]
    firstup=0
    cmm=[]
    rgrupup=false
    up_cond = 0
    # First upperbound GRASP
    global sol = []
    for _ in 1:50 global upperbound, sol
        tmp_sol=grasp_construction(kdata2,0,Kdata,Fnode2,ttm_dt,ttm)
        if tmp_sol[2]<upperbound
            sol = tmp_sol
        upperbound = tmp_sol[2]
        end
    end
    ub_update0 = false
    ub_update1 = false
    while factibility_condition!=1 ; 
        initial_it_time=Dates.now()
        itnumber+=1
        expanded_arcs=expanded_arcs_function(NE,Fnode,ttm_dt,LT,SKdata,Fnode2)
        global arc_ev=expanded_arcs[1];global arc_ek=expanded_arcs[2];global arcsk=expanded_arcs[3]
        model_creation=model_function(model_timelimit,arc_ev,arc_ek,ttm,NE,L0,LT,SKdata,arcsk,Ok,Dk)
        global model=model_creation[1];global x=model_creation[2];global y=model_creation[3]
        optimize!(model)
        cont_n+=1
        condition=0
        if !has_values(model)
            factibility_condition=1
            modelnovalues=true
            global model=last_solution[1];global x=last_solution[2];global y=last_solution[3];global NE=last_solution[4];
            global SKdata=last_solution[5];global arc_ev=last_solution[6];global arc_ek=last_solution[7];global arcsk=last_solution[8]
            println("LAST MIP MODEL NO SOLVED\n")
            end_status = 4
            continue
        else
            if objective_bound(model)!=lowerbound
                cont_upb=1
            else
                cont_upb+=1
            end
            llb=objective_bound(model)
            if objective_bound(model)>lowerbound
                lowerbound=objective_bound(model)
            end
            cm2=arcsused(y,arc_ev,N)
            if ((cont_upb > 3 && !(cm2 in cmm)) || (
                itnumber>5 && firstup==1 && !(cm2 in cmm)) || (
                    mod(itnumber,30)==0 && !(cm2 in cmm)))
                    println("Obtaining Upper Bound")
                    upp2=upperboundfunction(N,SKdata,ttm,Fnode2,initial_time,final_time,cm2,dt2,Fnode,
                    K,model_timelimit);
                    upp=upp2[1]
                    cm=upp2[2]
                    condupp=upp2[3]
                    if condupp
                        push!(cmm,cm)
                        modelup = upp2[4];yup=upp2[5];arc_evup=upp2[6];
                        LTup = upp2[7]; Kup = upp2[8];xup=upp2[9]
                        Kdataup=upp2[10];arcskup=upp2[11];Dkup=upp2[12];
                        arc_ekup=upp2[13];
                        firstup=0
                    end
                    if upp<upperbound
                        upperbound=upp
                        ub_update0 = true
                    end
                    sol2 = []
                    upt = Inf
                    for _ in 1:20
                        tmp_sol=grasp_construction(kdata2,1,Kdata,Fnode2,ttm_dt,ttm)
                        if tmp_sol[2]<upt
                            sol2 = tmp_sol
                            upt = tmp_sol[2]
                        end
                    end
                    if upt < upperbound
                        upperbound = upt
                        ub_update1 = true
                    end
            end
            global last_solution=lastsolution(model,x,y,NE,SKdata,arc_ev,arc_ek,arcsk)
            grouping_verification=grouping_verification_function(x,Ok,Dk,arcsk,SKdata)
            partial_expanded_verification=partial_expanded_verification_function(y,arc_ev)
            negarcscond = !partial_expanded_verification[4]
            courrent_time=Dates.value(Dates.now()-cpu_initial_time)/1000
            time_it=Dates.value(Dates.now()-initial_it_time)/1000
            if !partial_expanded_verification[1] & !grouping_verification[1]
                factibility_condition=1
                skinitial=length(SKdata)
                NEinitial=length(NE)
                push!(grouping_log,(cont_k,cont_n,llb,lowerbound,upperbound,(upperbound-lowerbound)/lowerbound*100,grouping_verification[2],grouping_verification[4],grouping_verification[3],skinitial,length(SKdata),time_it ))
                end_status = 2
                continue
            end
            if courrent_time>=maximum_timilimit
                factibility_condition=1
                timelimit_reached=true
                println("TIME LIMIT REACHED\n")
                end_status = 3
                continue
            end
            try global solution_check = fscheck(); catch ;  global solution_check=[false;[]] end
            if solution_check[1]
                factibility_condition=1
                println("SOLUTION\n")
                end_status = 1
                continue
            end
            if upperbound-lowerbound < epsilon
                factibility_condition=1
                println("UpperBound Solution\n")
                end_status = 5
                continue
            end
            if grouping_verification[1] && factibility_condition!=1 && negarcscond
                skinitial=length(SKdata)
                grouping_correction_function3(grouping_verification,SKdata,Kdata_dt,NE,Fnode,Dk,Ok,
                ttm_dt)
                push!(grouping_log,(cont_k,cont_n,llb,lowerbound,upperbound,(upperbound-lowerbound)/lowerbound*100,grouping_verification[2],grouping_verification[4],grouping_verification[3],skinitial,length(SKdata),time_it))
                cont_k+=1
            end
            if partial_expanded_verification[1] && factibility_condition!=1
                NEinitial=length(NE)
                partial_expanded_correction_function(partial_expanded_verification,NE,arc_ev,ttm_dt,
                Fnode2,initial_time)
                push!(partial_expanded_log,(itnumber,llb,lowerbound,upperbound,(upperbound-lowerbound)/lowerbound*100,NEinitial,length(NE),partial_expanded_verification[3],time_it,length(SKdata),grouping_verification[2],grouping_verification[3] ))
            elseif !partial_expanded_verification[1] && factibility_condition!=1
                NEinitial=length(NE)
                push!(partial_expanded_log,(itnumber,llb,lowerbound,upperbound,(upperbound-lowerbound)/lowerbound*100,NEinitial,length(NE),partial_expanded_verification[3],time_it,length(SKdata),grouping_verification[2],grouping_verification[3] ))
            end
        end
        println(itnumber)
    end
    total_time=Dates.value(Dates.now()-cpu_initial_time)/1000
    final_grouping=length(SKdata)
    if end_status == 1
        seq = solution_check[2];
        upp = lowerbound;
    elseif end_status == 2
        seq = solution(y,arc_ev,x,arc_ek,SKdata,final_time);
        upp = lowerbound;
    elseif end_status == 5 && ub_update0 
        seq = solutionup(yup,arc_evup,xup,arc_ekup,Kdataup);
        lowerbound = upp
    elseif end_status == 5 && ub_update1
        seq = sol[1]
        end_status = 6
        upp = upperbound
        lowerbound = upp
    elseif end_status == 5 && !ub_update1
        seq = sol[1]
        end_status = 7
        upp = upperbound
        lowerbound = upp
    else
        cm2=arcsused(y,arc_ev);
        upp2=upperboundfunction(N,SKdata,ttm,Fnode2,initial_time,final_time,cm2,dt2);
        upp=upp2[1];
        modelup = upp2[4];yup=upp2[5];arc_evup=upp2[6];
        LTup = upp2[7]; Kup = upp2[8];xup=upp2[9];
        Kdataup=upp2[10];arcskup=upp2[11];Dkup=upp2[12];arc_ekup=upp2[13];
        if has_values(modelup)
            seq = solutionup(yup,arc_evup,xup,arc_ekup,Kdataup);
        else
            seq = []
        end
        if upp > upperbound
            upp = upperbound
        end
        sol2 = []
        upt = Inf
        for _ in 1:20
            tmp_sol=grasp_construction(kdata2,1,Kdata,Fnode2,ttm_dt,ttm)
                if tmp_sol[2]<upt
                    sol2 = tmp_sol
                    upt = tmp_sol[2]
                end
        end
        if upt < upp
        upp = upt
        end
    end
    return [N-1,K,seq,upp,lowerbound,((upp-lowerbound)/upp),total_time,itnumber,end_status]
end
