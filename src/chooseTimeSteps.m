function S = chooseTimeSteps(S,substrate_time_advance)

if S.method=="global"
    S.pars.max_pde_dt = S.pars.max_global_pde_dt;
else
    S.pars.max_pde_dt = S.pars.max_local_pde_dt;
end

if S.method=="global" && S.pars.global_method_in_one_ode && S.pars.global_solver=="euler_direct"
    dts = [S.pars.max_pk_dt,S.pars.max_pde_dt,S.pars.max_exchange_dt,S.pars.max_internal_agent_dt,S.pars.max_dt];
    dt = min(dts,[],"omitnan");
    nsteps = ceil(substrate_time_advance/dt);
    S.solver.dt = substrate_time_advance/nsteps;
else
    nsteps = ceil(substrate_time_advance/S.pars.max_dt);
    S.solver.dt = substrate_time_advance/nsteps;
end

if S.method~="global" || ~S.pars.global_method_in_one_ode

    solvers_without_min_step = ["matlab_solver","matrix_exponential"];
    % any step using matlab solver or matrix exponential should not set the
    % time step
    if any(S.pars.pk_solver==solvers_without_min_step)
        S.solver.pk.num_steps = 1;
        S.solver.pk.dt = S.solver.dt;
    else
        S.solver.pk.num_steps = ceil(S.solver.dt/S.pars.max_pk_dt);
        S.solver.pk.dt = S.solver.dt/S.solver.pk.num_steps;
    end

    if S.method=="global" && any(S.pars.pde_solver_global_method==solvers_without_min_step)
        S.solver.pde.num_steps = 1;
        S.solver.pde.dt = S.solver.dt;
    else
        S.solver.pde.num_steps = ceil(S.solver.dt/S.pars.max_pde_dt);
        S.solver.pde.dt = S.solver.dt/S.solver.pde.num_steps;
    end

    if any(S.pars.exchange_solver==solvers_without_min_step)
        S.solver.exchange.num_steps = 1;
        S.solver.exchange.dt = S.solver.dt;
    else
        S.solver.exchange.num_steps = ceil(S.solver.dt/S.pars.max_exchange_dt);
        S.solver.exchange.dt = S.solver.dt/S.solver.exchange.num_steps;
    end

    if ~S.pars.use_internal_agent_ode || any(S.pars.internal_agent_solver==solvers_without_min_step)
        S.solver.internal_agent.num_steps = 1;
        S.solver.internal_agent.dt = S.solver.dt;
    else
        S.solver.internal_agent.num_steps = ceil(S.solver.dt/S.pars.max_internal_agent_dt);
        S.solver.internal_agent.dt = S.solver.dt/S.solver.internal_agent.num_steps;
    end

end