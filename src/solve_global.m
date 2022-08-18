function [S,internalized_substrates] = solve_global(S,...
    internalized_substrates,dt,cells_Lind)

if S.pars.is_cellular
    for type_ind_substrate = S.solver.ntypes:-1:1
        if isempty(internalized_substrates{type_ind_substrate}) % then either there are no cells of this type present or their internalization of this substrate is not tracked
            total_internalized_substrates(:,type_ind_substrate) = zeros(S.solver.nregions,1);
        else
            total_internalized_substrates(:,type_ind_substrate) = accumarray(S.solver.regions(cells_Lind{type_ind_substrate}),internalized_substrates{type_ind_substrate},[S.solver.nregions,1]);
        end
        cells_per_region(:,type_ind_substrate) = accumarray(S.solver.regions(cells_Lind{type_ind_substrate}),1,[S.solver.nregions,1]);
    end
    average_internalized_substrates = total_internalized_substrates./cells_per_region;
    average_internalized_substrates(cells_per_region==0) = 0;

    agent_density = cells_per_region./S.solver.region_volumes;
end

if S.pars.contains_dcs
    S.concentration(S.pars.dc_inds) = S.pars.dc_val;
end

if S.pars.global_method_in_one_ode

    IC = [S.concentration,average_internalized_substrates];

    switch S.pars.global_solver
        case "euler_direct"
            nt = ceil(dt/S.solver.dt);
            dt_temp = dt/nt;
            EC = IC(:);
            for ti = 1:nt
                EC = EC + dt_temp*globalODEFunction(EC,S.pars,...
                    S.solver.M_diffusion_and_efflux,agent_density,S.solver.sz);
            end
            EC = reshape(EC,S.solver.sz); % end conditions
        case "matlab_solver"
            sol = S.pars.global_matlab_solver(@(t,x) globalODEFunction(x,...
                S.pars,S.solver.M_diffusion_and_efflux,agent_density,S.solver.sz),[0 dt],IC(:));
            EC = reshape(sol.y(:,end),S.solver.sz); % end conditions
        case "matrix_exponential"
            if S.pars.use_internal_agent_ode || any(S.pars.export_rate~=0)
                error("not allowing for the matrix exponential if also using an internal agent ode and solving all equations at once; or if cells export this")
            end
            M = S.solver.M_global_template;
            for type_ind_substrate = S.solver.ntypes:-1:1
                M((1:S.solver.nregions),(1:S.solver.nregions)+S.solver.nregions*type_ind_substrate) = ...
                    diag(S.pars.secretion_rate(type_ind_substrate)*agent_density(:,type_ind_substrate));
                M((1:S.solver.nregions)+S.solver.nregions*type_ind_substrate,(1:S.solver.nregions)+S.solver.nregions*type_ind_substrate) = ...
                    diag(-S.pars.secretion_rate(type_ind_substrate)*ones(S.solver.nregions,1));
                M(1:S.solver.nregions,1:S.solver.nregions) = M(1:S.solver.nregions,1:S.solver.nregions) - ...
                    diag(S.pars.uptake_rate(type_ind_substrate)*agent_density(:,type_ind_substrate));
                M((1:S.solver.nregions)+S.solver.nregions*type_ind_substrate,(1:S.solver.nregions)) = ...
                    diag(S.pars.uptake_rate(type_ind_substrate)*ones(S.solver.nregions,1));
            end
            if S.pars.contains_dcs
                M(S.pars.dc_inds,:) = 0;
            end
            if S.pars.is_pk
                Minv_b = S.fluid_exchange_rate*S.circulation_concentration*(M\[S.solver.regional_BV_prop;zeros(S.solver.ntypes*S.solver.nregions,1)]);
            else
                Minv_b = zeros(S.solver.nregions,1);
            end
            c = IC(:) + Minv_b;
            EC = expm(M*dt)*c - Minv_b;
            EC = reshape(EC,S.solver.sz); % end conditions
            % If I want to, I could include the internal agent ode here in
            % the case where it is linear
            %             if C.use_internal_agent_ode
            %                 sol_internal = C.matlab_solver(@(t,x) C.internal_agent_ode(x,C.internal_agent_ode_pars),...
            %                     [0 dt],reshape(EC(:,2:end),[],1));
            %                 EC(:,2:end) = reshape(sol_internal.y(:,end),C.temp_sz);
            %                 % temp = sol_internal.y(:,end); % EC(:,2:end) = temp(:); % I think something like this may be necessary instead of the following line
            %                 %         EC(:,2:end) = sol_internal.y(:,end);
            %             end
    end

    S.concentration = EC(:,1);
    for type_ind_substrate = 1:S.solver.ntypes
        internalized_substrates{type_ind_substrate} = EC(S.solver.regions(cells_Lind{type_ind_substrate}),1+type_ind_substrate);
    end

    S.t = S.t + dt;

else % solve global method by sequentially solving the differential equations

    nt = ceil(dt/S.solver.dt);

    for ti = 1:nt
        %% solve blood exchange equation
        if S.pars.is_pk
            S.concentration = solve_pk(S.concentration,S.solver.pk);
        end

        %% solve diffusion/degradation equation
        switch S.pars.pde_solver_global_method
            case "euler_direct"
                for pde_ti = 1:S.solver.pde.num_steps
                    S.concentration = S.concentration + S.solver.pde.dt*S.solver.M_diffusion*S.concentration;
                end

            case "matlab_solver"
                sol = S.solver.pde.matlab_solver(@(t,x) S.solver.M_diffusion*x,[0 S.solver.dt],S.concentration);
                S.concentration = sol.y(:,end);

            case "matrix_exponential"
                S.concentration = S.solver.M_diffusion_exp * S.concentration;

        end

        %% solve cellular exchange
        if S.pars.is_exchange
            switch S.solver.exchange.solver
                % not using occupied regions to reduce computation
                case "euler_direct"
                    % the euler method here differs from how the below solvers
                    % work. Here, throughout the integration, the effects of
                    % all types on the region concentration is tracked,
                    % requiring the use of the agent density at each step. In other words,
                    % the contributions are averaged as if mimicking the diffusion process.
                    % This makes it a bit unfit for this because diffusion
                    % is handled elsewhere. Additionally, it is unlikely this 
                    % is ever used because tracking the
                    % concentration in a type-specific manner is much faster
                    % and the solution of this ODE can be analytically
                    % determined. Hence, I will usually solve the exchange
                    % dynamics using matrix_exponential.
                    for ex_ti = 1:S.solver.exchange.num_steps
                        sec_rate = S.pars.secretion_rate .* average_internalized_substrates;
                        up_rate = S.pars.uptake_rate .* S.concentration;
                        average_internalized_substrates = average_internalized_substrates + S.solver.exchange.dt*(up_rate-sec_rate);
                        S.concentration = S.concentration + S.solver.exchange.dt*sum((sec_rate-up_rate).*agent_density,2);
                    end

                    % If I ever want to use type-specific concentrations
                    % for exchange which is not done above, this should do it. one benchmarking
                    % exercise showed this to be actually marginally slower
                    % than the above. likely just from time spent
                    % allocating memory and the big sum after the loop
%                     conc_temp = repmat(S.concentration,1,S.solver.ntypes);
%                     for ex_ti = 1:S.num_exchange_steps
%                         sec_rate = S.pars.secretion_rate * average_internalized_substrates;
%                         up_rate = S.pars.uptake_rate * conc_temp;
%                         a_i_s_change = S.solver.exchange.dt*(up_rate-sec_rate);
%                         average_internalized_substrates = average_internalized_substrates + a_i_s_change;
%                         conc_temp = conc_temp - a_i_s_change;
%                     end
%                     S.concentration = S.concentration + sum(agent_density.*(conc_temp-S.concentration),2);

                case "matlab_solver"
                    sol = S.solver.exchange.matlab_solver(@(t,x) globalExchangeFunction(x,S.pars.secretion_rate,S.pars.uptake_rate,[size(agent_density),2]),[0 S.solver.dt],reshape([repmat(S.concentration,S.solver.ntypes,1),average_internalized_substrates(:)],[],1));
                    EC = sol.y(:,end);
                    EC = reshape(EC,[size(agent_density),2]);
                    a_i_s_change = EC(:,:,2) - average_internalized_substrates;
                    average_internalized_substrates = EC(:,:,2);
                    S.concentration = S.concentration - sum(agent_density.*a_i_s_change,2);

                case "matrix_exponential"
                    D = S.pars.uptake_rate.*(average_internalized_substrates+S.concentration)./(S.pars.secretion_rate+S.pars.uptake_rate);
                    const_C = average_internalized_substrates-D;
                    a_i_s_change = const_C.*exp(-(S.pars.secretion_rate+S.pars.uptake_rate)*S.solver.dt) + D - average_internalized_substrates;
                    average_internalized_substrates = average_internalized_substrates + a_i_s_change;
                    S.concentration = S.concentration - sum(agent_density.*a_i_s_change,2);
            end
        end

        if any(S.pars.export_rate~=0)
            S.concentration = S.concentration + S.solver.dt*sum(S.pars.export_rate .* agent_density,2); % this is the exact solution to the export equation: c' = export_rate * agent_density
        end

        %% internal agent ODEs
        if S.pars.use_internal_agent_ode
            switch S.solver.internal_agent.solver
                case "euler_direct"
                    for int_ti = 1:S.solver.internal_agent.num_steps
                        average_internalized_substrates = average_internalized_substrates + S.solver.internal_agent.dt*reshape(S.pars.internal_agent_ode(average_internalized_substrates,1:S.solver.ntypes),[],S.solver.ntypes);
                    end

                case "matlab_solver"
                    sol_internal = S.solver.internal_agent.matlab_solver(@(t,x) S.pars.internal_agent_ode(x,1:S.solver.ntypes),...
                        [0 S.solver.dt],average_internalized_substrates);
                    average_internalized_substrates(:) = sol_internal.y(:,end);

                case "matrix_exponential"
                    error("not allowing for the matrix exponential if also using an internal agent ode and solving all equations at once")
                    average_internalized_substrates = S.M_internal_agent_eV * diag(exp(S.M_internal_agent_ev*S.solver.dt)) * S.M_internal_agent_eV_inv * average_internalized_substrates;

            end
        end

    end

    for type_ind_substrate = 1:S.solver.ntypes
        if ~isempty(internalized_substrates{type_ind_substrate})
            internalized_substrates{type_ind_substrate} = average_internalized_substrates(S.solver.regions(cells_Lind{type_ind_substrate}),type_ind_substrate);
        end
    end

    S.t = S.t + nt*S.solver.dt;


end