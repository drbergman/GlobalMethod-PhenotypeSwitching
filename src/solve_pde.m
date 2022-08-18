function C = solve_pde(C,pde_pars)

switch ndims(C)
    case 2
        for ti = 1:pde_pars.num_steps

            % set DCs if exist
            if pde_pars.contains_dcs
                C(pde_pars.dc_inds) = pde_pars.dc_val;
            end

            % x-diffusion
            for j = 1:pde_pars.ny
                C(:,j) = pde_pars.Mx\C(:,j);
            end

            % set DCs if exist
            if pde_pars.contains_dcs
                C(pde_pars.dc_inds) = pde_pars.dc_val;
            end

            % y-diffusion
            C = C';
            for i = 1:pde_pars.nx
                C(:,i) = pde_pars.My\C(:,i);
            end
            C = C';

            % set DCs if exist
            if pde_pars.contains_dcs
                C(pde_pars.dc_inds) = pde_pars.dc_val;
            end
        end

    case 3
        for ti = 1:pde_pars.num_steps
            % x-diffusion
            for j = 1:pde_pars.ny
                for k = 1:pde_pars.nz
                    C(:,j,k) = pde_pars.Mx\C(:,j,k);
                end
            end

            % set DCs if exist
            if pde_pars.contains_dcs
                C(pde_pars.dc_inds) = pde_pars.dc_val;
            end

            % y-diffusion
            C = permute(C,[2,3,1]); % set up in order [y,z,x]
            for i = 1:pde_pars.nx
                for k = 1:pde_pars.nz
                    C(:,k,i) = pde_pars.My\C(:,k,i);
                end
            end

            % set DCs if exist
            if pde_pars.contains_dcs
                C = permute(C,[3,1,2]); % return C to original dimension ordering
                C(pde_pars.dc_inds) = pde_pars.dc_val;
                C = permute(C,[3,1,2]); % set up in order [z,x,y]
            else
                C = permute(C,[2,3,1]); % set up in order [z,x,y]
            end

            % z-diffusion
            for i = 1:pde_pars.nx
                for j = 1:pde_pars.ny
                    C(:,i,j) = pde_pars.Mz\C(:,i,j);
                end
            end
            C = permute(C,[2,3,1]); % return C to original dimension ordering

            % set DCs if exist
            if pde_pars.contains_dcs
                C(pde_pars.dc_inds) = pde_pars.dc_val;
            end
        end

end