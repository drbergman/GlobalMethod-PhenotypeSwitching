function C = solve_pk(C,pk_pars)

switch pk_pars.solver
    case "euler_direct"
        C_BV = C(pk_pars.region_contains_bv);
        for ti = 1:pk_pars.num_steps
            C_BV = C_BV + pk_pars.dt * pk_pars.fluid_exchange_rate_bv_regions .* (pk_pars.circulation_concentration - C_BV); % distribution and redistribution
        end
        C(pk_pars.region_contains_bv) = C_BV;
    case "matrix_exponential"
        C(pk_pars.region_contains_bv) = (C(pk_pars.region_contains_bv)-pk_pars.circulation_concentration).*exp(-pk_pars.fluid_exchange_rate_bv_regions*pk_pars.dt) + pk_pars.circulation_concentration;
    case "matlab_solver"
        sol = pk_pars.matlab_solver(@(t,x) pk_pars.fluid_exchange_rate_bv_regions.*(pk_pars.circulation_concentration-x),...
            [0 pk_pars.dt],C(pk_pars.region_contains_bv));
        C(pk_pars.region_contains_bv) = sol.y(:,end);
    otherwise
        error("have not set up this pk solver yet")

end