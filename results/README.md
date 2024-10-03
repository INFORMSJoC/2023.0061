# RESULTS

Results files are divided in two:

- First, the folders named `end_uni_dt_dk` contain solutions for the uniform distribution case (uni) with $\delta_t = dt$ and $\delta_k = dk$. Each folder describes the evolution of the grouping and network discovery procedures.
    - cplexlog
    - CG_log.csv = Grouping procdure log
    - PE_log.csv = Dynamic network discovery log
    - Iterations_log.csv = Algortithm iterations log 
    - te_solution.csv files gives presents:
        - Number of SCCs
        - Number of TR
        - Lowerbound value
        - Uperbound value
        - Optimality gap
        - CPU time (s)
        - Number of iterations
        - Termination criteria
- Second, folders named `tr_g_end_dis_dt_dk` contain files te_solution.csv that include:
    - Whether to allow (tr1) or not (tr2) transfers
    - Whether to use (g1) or not (g0) the grouping procedure
    - The type of temporal distribution of the commodity generation, namely: uniform (uni), distribution 1 (dist1) and distribution 2 (dist2)
    - Time discretization (dt)
    - Comodities generation (dk) 
