function [S,internalized_substrates] = solve_local(S,...
    internalized_substrates,dt,cells_Lind)

all_cells_Lind = cat(1,cells_Lind{:});

nt = ceil(dt/S.solver.dt);

for ti = 1:nt

    %% solve blood exchange equation
    if S.pars.is_pk
        S.concentration = solve_pk(S.concentration,S.solver.pk);
    end

    %% solve diffusion/degradation equation
    S.concentration = solve_pde(S.concentration,S.solver.pde);

    %% Cell supply/uptake
    if isempty(all_cells_Lind)
        continue;
    end

    if S.pars.is_exchange
        switch S.solver.exchange.solver
            case "euler_direct"
                for ex_ti = 1:S.solver.exchange.num_steps
%                     for type_ind = S.solver.ntypes:-1:1
%                         d_IS{type_ind} = [internalized_substrates{type_ind},S.concentration(cells_Lind{type_ind})] * S.solver.exchange.M_exchange_agent(1,:)';  % the rate of change of the concentration in the microenvironment is just the opposite of this becuase we assume they have the same volume; if they don't, then only need to multiply by the ratio of these volumes (cell volume / lattice site volume)
%                         internalized_substrates{type_ind} = internalized_substrates{type_ind} + S.solver.exchange.dt * d_IS{type_ind};
%                         sec_rate{type_ind} = S.pars.secretion_rate * internalized_substrates{type_ind};
%                         up_rate{type_ind} = S.pars.uptake_rate * S.concentration(cells_Lind{type_ind});
%                         internalized_substrates{type_ind} = internalized_substrates{type_ind} + S.solver.exchange.dt * (up_rate{type_ind} - sec_rate{type_ind});
%                     end
%                     S.concentration = S.concentration + S.solver.exchange.dt * reshape(accumarray(all_cells_Lind,-cat(1,d_IS{:}),[numel(S.concentration),1]),size(S.concentration));
%                     S.concentration = S.concentration + S.solver.exchange.dt * reshape(accumarray(all_cells_Lind,cat(1,sec_rate{:})-cat(1,up_rate{:}),[numel(S.concentration),1]),size(S.concentration));
                    
                    for type_ind_substrate = S.solver.ntypes:-1:1
                        if S.solver.is_exchange_type(type_ind_substrate)
                            d_IS = S.solver.exchange.M_exchange_agent(1,:,type_ind_substrate) * [internalized_substrates{type_ind_substrate},S.concentration(cells_Lind{type_ind_substrate})];  % the rate of change of the concentration in the microenvironment is just the opposite of this becuase we assume they have the same volume; if they don't, then only need to multiply by the ratio of these volumes (cell volume / lattice site volume)
                            internalized_substrates{type_ind_substrate} = internalized_substrates{type_ind_substrate} + S.solver.exchange.dt * d_IS;
                            S.concentration = S.concentration + S.solver.exchange.dt * reshape(accumarray(cells_Lind{type_ind_substrate},-d_IS,[numel(S.concentration),1]),size(S.concentration));
                        end
                    end
                end
            case "matlab_solver"
                for type_ind_substrate = S.solver.ntypes:-1:1
                    for ai = 1:length(cells_Lind{type_ind_substrate})
                        sol = S.solver.exchange.matlab_solver(@(t,x) S.solver.exchange.M_exchange_agent(:,:,type_ind_substrate)*x,...
                            [0,S.solver.dt],[internalized_substrates{type_ind_substrate}(ai);S.concentration(cells_Lind{type_ind_substrate}(ai))]);
                        internalized_substrates{type_ind_substrate}(ai) = sol.y(1,end);
                        S.concentration(cells_Lind{type_ind_substrate}(ai)) = sol.y(2,end);
                    end
                end
            case "matrix_exponential"
                for type_ind_substrate = S.solver.ntypes:-1:1
                    if ~isempty(internalized_substrates{type_ind_substrate})
                        c_temp = S.solver.exchange.M_exchange_agent_eV(:,:,type_ind_substrate) * diag(exp(S.solver.exchange.M_exchange_agent_ev(:,type_ind_substrate)*S.solver.dt)) * S.solver.exchange.M_exchange_agent_eV_inv(:,:,type_ind_substrate) * [internalized_substrates{type_ind_substrate}';S.concentration(cells_Lind{type_ind_substrate})'];
                        internalized_substrates{type_ind_substrate}(:) = c_temp(1,:);
                        S.concentration(cells_Lind{type_ind_substrate}) = c_temp(2,:);
                    end
                end
        end
    end

    if any(S.pars.export_rate~=0)
        for type_ind_substrate = 1:S.pars.ntypes
            S.concentration(cells_Lind{type_ind_substrate}) = S.concentration(cells_Lind{type_ind_substrate}) + S.solver.dt*S.pars.export_rate(type_ind_substrate);
        end
    end

    %% internal agent ODEs
    if S.pars.use_internal_agent_ode
        switch S.solver.internal_agent.solver
            case "euler_direct"
                for type_ind_substrate = S.solver.ntypes:-1:1
                    for int_ti = 1:S.solver.internal_agent.num_steps
                        internalized_substrates{type_ind_substrate} = internalized_substrates{type_ind_substrate} + S.solver.internal_agent.dt * S.pars.internal_agent_ode(internalized_substrates{type_ind_substrate},type_ind_substrate);
                    end
                end
            case "matlab_solver"
                for type_ind_substrate = S.solver.ntypes:-1:1
                    for ai = 1:length(cells_Lind{type_ind_substrate})
                        sol_internal = S.solver.internal_agent.matlab_solver(@(t,x) S.pars.internal_agent_ode(x,type_ind_substrate),...
                            [0,S.solver.dt],internalized_substrates{type_ind_substrate}(ai));
                        internalized_substrates{type_ind_substrate}(ai) = sol_internal.y(end);
                    end
                end
            case "matrix_exponential"
                for type_ind_substrate = S.solver.ntypes:-1:1
                    internalized_substrates{type_ind_substrate} = (S.M_internal_agent_eV * diag(exp(S.M_internal_agent_ev*S.solver.dt)) * S.M_internal_agent_eV_inv * (internalized_substrates{type_ind_substrate}'))';
                end
        end

    end

end

% if length(internalized_substrates)==1
%     internalized_substrates = internalized_substrates{1};
% end

S.t = S.t + nt*S.solver.dt;