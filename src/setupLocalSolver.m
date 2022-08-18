function S = setupLocalSolver(S,TME_size,BV,substrate_time_advance)

S = chooseTimeSteps(S,substrate_time_advance);
S.solver.regions = [];
S.solver.region_volumes = [];
S.solver.ndims = length(TME_size);


if S.pars.is_pk
    if S.pars.contains_dcs
        warning("make sure the Dirichlet conditions and blood vessels do not overlap otherwise the fluid exchange could violate the Dirichlet conditions")
    end
    S.solver.pk.solver = S.pars.pk_solver;
    S.solver.pk.region_contains_bv = BV>0;
    S.solver.pk.fluid_exchange_rate_bv_regions = S.pars.fluid_exchange_rate * BV(S.solver.pk.region_contains_bv);
    S.solver.pk.circulation_concentration = S.pars.circulation_concentration;
    if S.solver.pk.solver == "matlab_solver"
        S.solver.pk.matlab_solver = S.pars.pk_matlab_solver;
    end
end

S.solver.pde = setupPDESolver(TME_size,S.solver.pde,S.pars.diffusion,S.pars.degradation);
S.solver.pde.contains_dcs = S.pars.contains_dcs;
if S.solver.pde.contains_dcs
    S.solver.pde.dc_inds = setupDC(S,TME_size,BV);
    S.solver.pde.dc_val = S.pars.dc_val;
end

S.solver.ntypes = S.pars.ntypes;

if S.pars.ntypes~=0 && (any(S.pars.uptake_rate ~=0) || any(S.pars.secretion_rate ~= 0))
    if S.pars.contains_dcs
        warning("make sure the exchanging agents cannot end up on the Dirichlet conditions or else they may violate those conditions")
    end
    S.pars.is_exchange = true;
    S.solver.exchange.solver = S.pars.exchange_solver;
    for type_ind_substrate = 1:S.solver.ntypes
        S.solver.exchange.M_exchange_agent(:,:,type_ind_substrate) = [S.pars.secretion_rate(type_ind_substrate)*[-1;1],S.pars.uptake_rate(type_ind_substrate)*[1;-1]];
        switch S.solver.exchange.solver
            case "matlab_solver"
                S.solver.exchange.matlab_solver = S.pars.exchange_matlab_solver;
            case "matrix_exponential"
                [S.solver.exchange.M_exchange_agent_eV(:,:,type_ind_substrate),ev_temp] = eig(S.solver.exchange.M_exchange_agent(:,:,type_ind_substrate));
                S.solver.exchange.M_exchange_agent_ev(:,type_ind_substrate) = diag(ev_temp);
                S.solver.exchange.M_exchange_agent_eV_inv(:,:,type_ind_substrate) = inv(S.solver.exchange.M_exchange_agent_eV(:,:,type_ind_substrate));
        end
    end
end

if any(S.pars.export_rate~=0) && S.pars.contains_dcs
    warning("make sure the exporting agents cannot end up on the Dirichlet conditions or else they may violate those conditions")
end

if S.pars.use_internal_agent_ode
    S.solver.internal_agent.solver = S.pars.internal_agent_solver;
    if S.solver.internal_agent.solver=="matlab_solver"
        S.solver.internal_agent.matlab_solver = S.pars.internal_agent_matlab_solver;
    end
end


S.setup_done = true;