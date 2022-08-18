function [S,internalized_substrates] = solve_hybrid(S,...
    internalized_substrates,dt,cells_Lind)

for type_ind = S.pars.ntypes:-1:1
    if isempty(cells_Lind{type_ind})
        total_internalized_substrates(:,type_ind) = zeros(S.solver.nregions,1);
        cells_per_region(:,type_ind) = zeros(S.solver.nregions,1);
    else
        total_internalized_substrates(:,type_ind) = accumarray(S.solver.regions(cells_Lind{type_ind}),internalized_substrates{type_ind},[S.solver.nregions,1]);
        cells_per_region(:,type_ind) = accumarray(S.solver.regions(cells_Lind{type_ind}),1,[S.solver.nregions,1]);
    end
end
average_internalized_substrates = total_internalized_substrates./cells_per_region;
average_internalized_substrates(cells_per_region==0) = 0;


IC = zeros(2,S.pars.ntypes,S.solver.nregions);
IC(1,:,:) = average_internalized_substrates';

nt = ceil(dt/S.solver.dt);

for ri = S.solver.nregions:-1:1
    for type_ind = S.pars.ntypes:-1:1
        agents_in_region_by_type{type_ind,ri} = find(S.solver.regions(cells_Lind{type_ind})==ri);
        C_ind{type_ind,ri} = cells_Lind{type_ind}(agents_in_region_by_type{type_ind,ri});
        type_in_region(type_ind,ri) = ~isempty(agents_in_region_by_type{type_ind,ri});
    end
end

for ti = 1:nt

    %% solve blood exchange equation
    if S.pars.is_pk
        S.concentration = solve_pk(S.concentration,S.solver.pk);
    end

    %% solve diffusion/degradation equation
    S.concentration = solve_pde(S.concentration,S.solver.pde);

    %% collect averages in regions
    for ri = S.solver.nregions:-1:1
        for type_ind = S.pars.ntypes:-1:1
            if type_in_region(type_ind,ri)
                IC(2,type_ind,ri) = mean(S.concentration(C_ind{type_ind,ri}));
            end
        end
    end

    %% solve cellular exchange
    if S.pars.is_exchange
        IC = permute(IC,[1,3,2]); % rearrange so each page is a 2xnregions containing all the info for a particular cell type
%         IC = reshape(IC,2,[]);
%         IC_occupied = IC(:,type_in_region(:));
        switch S.solver.exchange.solver
            case "euler_direct"
                for type_ind = 1:S.pars.ntypes
                    for ex_ti = 1:S.num_exchange_steps
                        IC(:,:,type_ind) = IC(:,:,type_ind) + S.exchange_dt*S.solver.exchange.M_exchange_agent(:,:,type_ind)*IC(:,:,type_ind);
                    end
                end

            case "matlab_solver"
                for type_ind = 1:S.pars.ntypes
                    for ri = 1:S.solver.nregions
                        sol = S.exchange_matlab_solver(@(t,x) S.solver.exchange.M_exchange_agent(:,:,type_ind)*x,[0 S.solver.dt],IC(:,ri,type_ind));
                        IC_occupied(:,ci) = sol.y(:,end);
                    end
                end

            case "matrix_exponential"
                for type_ind = 1:S.pars.ntypes
                    IC(:,:,type_ind) = S.solver.exchange.M_exchange_agent_eV(:,:,type_ind) * diag(exp(S.solver.exchange.M_exchange_agent_ev(:,type_ind)*S.solver.dt)) * S.solver.exchange.M_exchange_agent_eV_inv(:,:,type_ind) * IC(:,:,type_ind);
                end
        end
        IC = permute(IC,[1,3,2]);
%         IC(:,type_in_region(:)) = IC_occupied;
%         IC = reshape(IC,2,S.pars.ntypes,S.solver.nregions);
    end

    for type_ind = S.pars.ntypes:-1:1
        for ri = S.solver.nregions:-1:1
            if type_in_region(type_ind,ri)
                S.concentration(C_ind{type_ind,ri}) = IC(2,type_ind,ri) + S.solver.dt*S.pars.export_rate(type_ind);
            end
        end
    end

    %% internal agent ODEs
    if S.pars.use_internal_agent_ode
        IC = reshape(IC,2,[]);
        IC_occupied_cell = IC(1,type_in_region(:))';
        switch S.solver.internal_agent.solver
            case "euler_direct"
                for int_ti = 1:S.solver.internal_agent.num_steps
                    IC_occupied_cell = IC_occupied_cell + S.solver.internal_agent.dt*S.pars.internal_agent_ode(IC_occupied_cell,1);
                end

            case "matlab_solver"
                for ci = 1:length(IC_occupied_cell)
                    sol_internal = S.internal_agent_matlab_solver(@(t,x) S.internal_agent_ode(x,S.internal_agent_ode_pars),...
                        [0 S.solver.dt],IC_occupied_cell(ci));
                    IC_occupied_cell(ci) = sol_internal.y(end);
                end

            case "matrix_exponential"
                IC_occupied_cell = S.M_internal_agent_eV * diag(exp(S.M_internal_agent_ev*S.solver.dt)) * S.M_internal_agent_eV_inv * IC_occupied_cell;
        end
        IC(1,type_in_region(:)) = IC_occupied_cell;
        IC = reshape(IC,2,S.pars.ntypes,S.solver.nregions);
    end
end

if S.pars.ntypes==1
    internalized_substrates = reshape(IC(1,1,S.solver.regions(cells_Lind{1})),[],1);
else
    for type_ind = 1:S.pars.ntypes
        internalized_substrates{type_ind} = reshape(IC(1,type_ind,S.solver.regions(cells_Lind{type_ind})),[],1);
    end
end

S.t = S.t + nt*S.solver.dt;
