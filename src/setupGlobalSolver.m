function S = setupGlobalSolver(S,TME_size,bv_location,BV,substrate_time_advance)

S = chooseTimeSteps(S,substrate_time_advance);

if S.pars.contains_dcs
    S.pars.dc_inds = setupDC(S,TME_size,BV); % dc inds here will correspond to location in cellular mesh, to be (possibly) used in setupRegions
end

S.solver.regions = setupRegions(S,TME_size,bv_location,BV);

if S.pars.contains_dcs
    S.pars.dc_inds = unique(S.solver.regions(S.pars.dc_inds)); % now switch dc inds to reference the region ind they are in (consider adding new field dc_region_inds instead)
end

S.solver.region_volumes = accumarray(S.solver.regions(:),1);
S.solver.nregions = length(S.solver.region_volumes);
S.solver.ndims = length(TME_size);

switch S.solver.ndims
    case 2
        % the row index (first column of the first input) is the neighboring
        % region. the column index (second column of the first input) is the
        % center region.
        M = accumarray([reshape(S.solver.regions(1:end-1,:),[],1),reshape(S.solver.regions(2:end,:),[],1);... % count neighbors to left
            reshape(S.solver.regions(2:end,:),[],1),reshape(S.solver.regions(1:end-1,:),[],1);... % count neighbors to right
            reshape(S.solver.regions(:,1:end-1),[],1),reshape(S.solver.regions(:,2:end),[],1);... % count neighbors below
            reshape(S.solver.regions(:,2:end),[],1),reshape(S.solver.regions(:,1:end-1),[],1);... % count neighbors above
            S.solver.regions(1,:)',S.solver.regions(1,:)';... % count neighbors left of left edge as same region as center
            S.solver.regions(end,:)',S.solver.regions(end,:)';... % count neighbors right of right edge as same region as center
            S.solver.regions(:,1),S.solver.regions(:,1);... % count neighbors below bottom edge as same region as center
            S.solver.regions(:,end),S.solver.regions(:,end);... % count neighbors above top edge as same region as center
            ],1);

        % Since M is symmetric, it is irrevelant (in these loops) whether you
        % consider rows as the region of the center point or of the neighbor
        % point.
        % Given how it is used in globalODEFunction.m (the global method S.solver), i
        % corresponds to the neighboring region and j all the possible
        % neighbors in that region. Then, the four possible center points which
        % have said spot as a neighbor are checked for their region.

        % this was the loops that do the above perhaps more clearly, but slower
        %     uniq_regions = unique(C.regions);
        %     n_regions = length(uniq_regions);
        %     M = zeros(n_regions);

        %     for i = 1:n_regions
        %
        %         [x,y] = find(C.regions==i);
        %
        %         for j = 1:length(x)
        %
        %             if x(j)>1
        %                 M(i,C.regions(x(j)-1,y(j))) = M(i,C.regions(x(j)-1,y(j))) + 1;
        %             else % Neumann condition can be interpreted as neighbors off the lattice are the same region as the center point (so there is no derivative there)
        %                 M(i,i) = M(i,i)+1;
        %             end
        %             if x(j)<size(C.regions,1)
        %                 M(i,C.regions(x(j)+1,y(j))) = M(i,C.regions(x(j)+1,y(j))) + 1;
        %             else % Neumann condition can be interpreted as neighbors off the lattice are the same region as the center point (so there is no derivative there)
        %                 M(i,i) = M(i,i)+1;
        %             end
        %
        %             if y(j)>1
        %                 M(i,C.regions(x(j),y(j)-1)) = M(i,C.regions(x(j),y(j)-1)) + 1;
        %             else % Neumann condition can be interpreted as neighbors off the lattice are the same region as the center point (so there is no derivative there)
        %                 M(i,i) = M(i,i)+1;
        %             end
        %             if y(j)<size(C.regions,2)
        %                 M(i,C.regions(x(j),y(j)+1)) = M(i,C.regions(x(j),y(j)+1)) + 1;
        %             else % Neumann condition can be interpreted as neighbors off the lattice are the same region as the center point (so there is no derivative there)
        %                 M(i,i) = M(i,i)+1;
        %             end
        %         end
        %     end

    case 3

        M = accumarray([reshape(S.solver.regions(1:end-1,:,:),[],1),reshape(S.solver.regions(2:end,:,:),[],1);... % count neighbors to left
            reshape(S.solver.regions(2:end,:,:),[],1),reshape(S.solver.regions(1:end-1,:,:),[],1);... % count neighbors to right
            reshape(S.solver.regions(:,1:end-1,:),[],1),reshape(S.solver.regions(:,2:end,:),[],1);... % count neighbors in front
            reshape(S.solver.regions(:,2:end,:),[],1),reshape(S.solver.regions(:,1:end-1,:),[],1);... % count neighbors behind
            reshape(S.solver.regions(:,:,1:end-1),[],1),reshape(S.solver.regions(:,:,2:end),[],1);... % count neighbors below
            reshape(S.solver.regions(:,:,2:end),[],1),reshape(S.solver.regions(:,:,1:end-1),[],1);... % count neighbors above
            reshape(S.solver.regions(1,:,:),[],1),reshape(S.solver.regions(1,:,:),[],1);... % count neighbors left of left edge as same region as center
            reshape(S.solver.regions(end,:,:),[],1),reshape(S.solver.regions(end,:,:),[],1);... % count neighbors right of right edge as same region as center
            reshape(S.solver.regions(:,1,:),[],1),reshape(S.solver.regions(:,1,:),[],1);... % count neighbors in front of front edge as same region as center
            reshape(S.solver.regions(:,end,:),[],1),reshape(S.solver.regions(:,end,:),[],1);... % count neighbors behind back edge as same region as center
            reshape(S.solver.regions(:,:,1),[],1),reshape(S.solver.regions(:,:,1),[],1);... % count neighbors below bottom edge as same region as center
            reshape(S.solver.regions(:,:,end),[],1),reshape(S.solver.regions(:,:,end),[],1);... % count neighbors above top edge as same region as center
            ],1);

end

assert(issymmetric(M))

S.solver.M = M./sum(M,1);

if S.pars.is_pk
    if S.pars.contains_dcs
        warning("make sure the Dirichlet conditions and blood vessels do not overlap otherwise the fluid exchange could violate the Dirichlet conditions")
    end
    S.solver.regional_BV_prop = accumarray(S.solver.regions(:),BV(:))./S.solver.region_volumes;
    S.solver.pk.region_contains_bv = S.solver.regional_BV_prop>0;
    if S.pars.global_method_in_one_ode
        S.pars.fluid_exchange_rate_all_regions = S.pars.fluid_exchange_rate * S.solver.regional_BV_prop;
    else
        S.pars.fluid_exchange_rate_bv_regions = S.pars.fluid_exchange_rate * S.solver.regional_BV_prop(S.solver.pk.region_contains_bv);
    end
else
    S.solver.regional_BV_prop = zeros(S.solver.nregions,1);
    S.pars.fluid_exchange_rate_all_regions = zeros(S.solver.nregions,1);
end

if S.pars.global_method_in_one_ode
    if S.pars.global_solver=="matrix_exponential"
        S.solver.M_global_template = zeros(S.solver.nregions*(S.solver.ntypes+1));
        S.solver.M_global_template(1:S.solver.nregions,1:S.solver.nregions) = coefficientMatrixRegions(S,S.solver.ndims);
    else
        S.solver.M_diffusion_and_efflux = coefficientMatrixRegions(S,S.solver.ndims);
    end
else
    if S.pars.is_pk
        S.solver.pk.solver = S.pars.pk_solver;
        S.solver.pk.fluid_exchange_rate_bv_regions = S.pars.fluid_exchange_rate_bv_regions;
        S.solver.pk.circulation_concentration = S.pars.circulation_concentration;
        if S.solver.pk.solver == "matlab_solver"
            S.solver.pk.matlab_solver = S.pars.pk_matlab_solver;
        end
    end
    
    S.pars.fluid_exchange_rate_all_regions = zeros(S.solver.nregions,1); % don't do pk dynamics when solving the diffusion equation
    S.solver.M_diffusion = coefficientMatrixRegions(S,S.solver.ndims);
    switch S.pars.pde_solver_global_method
        case "matlab_solver"
            S.solver.pde.matlab_solver = S.pars.pde_matlab_solver;
        case "matrix_exponential"
            S.solver.M_diffusion_exp = expm(S.solver.M_diffusion*S.solver.dt);
    end

    if S.pars.ntypes~=0 && (any(S.pars.uptake_rate ~=0) || any(S.pars.secretion_rate ~= 0))
        S.pars.is_exchange = true;
        S.solver.exchange.solver = S.pars.exchange_solver;
        switch S.solver.exchange.solver
            case "matlab_solver"
                S.solver.exchange.matlab_solver = S.pars.exchange_matlab_solver;
        end
    end

    if S.pars.use_internal_agent_ode
        S.solver.internal_agent.solver = S.pars.internal_agent_solver;
        if S.solver.internal_agent.solver=="matlab_solver"
            S.solver.internal_agent.matlab_solver = S.pars.internal_agent_matlab_solver;
        end
    end

end

S.pars.is_cellular = S.pars.ntypes~=0 && (any(S.pars.uptake_rate~=0) || any(S.pars.secretion_rate~=0) || any(S.pars.export_rate~=0) || S.pars.use_internal_agent_ode);

if S.pars.is_cellular && S.pars.contains_dcs
    warning("make sure the exchanging/exporting agents cannot end up on the Dirichlet conditions or else they may violate those conditions")
end

if any(S.solver.regional_BV_prop>0 & S.solver.regional_BV_prop<1)
%     error("currently only have global method set up for these bv proportions to be 0 or 1")
end

S.solver.ntypes = S.pars.ntypes;

S.solver.sz = [S.solver.nregions,1+S.solver.ntypes];

if S.pars.use_internal_agent_ode
    S.solver.temp_sz = [S.solver.nregions,S.solver.ntypes];
end

S.setup_done = true;