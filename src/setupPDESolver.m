function pde_pars = setupPDESolver(TME_size,pde_pars,diffusion,degradation)

pde_pars.nx = TME_size(1);
pde_pars.ny = TME_size(2);
if length(TME_size)==3
    pde_pars.nz = TME_size(3);
    ndims = 3;
else
    ndims = 2;
end

% note that diffusion has already been normalized by the lattice size

%% set up in x direction
r_ind_ux = 1:pde_pars.nx-1;
c_ind_ux = 2:pde_pars.nx;
v_ux = -pde_pars.dt*diffusion * ones(1,pde_pars.nx-1);

r_ind_dx = 1:pde_pars.nx;
c_ind_dx = 1:pde_pars.nx;
v_dx(1,2:pde_pars.nx-1) = 1 + pde_pars.dt*degradation/ndims + 2*pde_pars.dt*diffusion;
v_dx([1,pde_pars.nx]) =   1 + pde_pars.dt*degradation/ndims +   pde_pars.dt*diffusion;

r_ind_lx = 2:pde_pars.nx;
c_ind_lx = 1:pde_pars.nx-1;
v_lx = -pde_pars.dt*diffusion * ones(1,pde_pars.nx-1);

pde_pars.Mx = sparse([r_ind_ux,r_ind_dx,r_ind_lx],...
    [c_ind_ux,c_ind_dx,c_ind_lx],...
    [v_ux    ,v_dx   ,v_lx     ],pde_pars.nx,pde_pars.nx);

%% set up in y direction
r_ind_uy = 1:pde_pars.ny-1;
c_ind_uy = 2:pde_pars.ny;
v_uy = -pde_pars.dt*diffusion * ones(1,pde_pars.ny-1);

r_ind_dy = 1:pde_pars.ny;
c_ind_dy = 1:pde_pars.ny;
v_dy(1,2:pde_pars.ny-1) = 1 + pde_pars.dt*degradation/ndims + 2*pde_pars.dt*diffusion;
v_dy([1,pde_pars.ny]) = 1 + pde_pars.dt*degradation/ndims + pde_pars.dt*diffusion;

r_ind_ly = 2:pde_pars.ny;
c_ind_ly = 1:pde_pars.ny-1;
v_ly = -pde_pars.dt*diffusion * ones(1,pde_pars.ny-1);

pde_pars.My = sparse([r_ind_uy,r_ind_dy,r_ind_ly],...
    [c_ind_uy,c_ind_dy,c_ind_ly],...
    [v_uy    ,v_dy   ,v_ly    ]);

if ndims==3
    %% set up in z direction
    r_ind_uz = 1:pde_pars.nz-1;
    c_ind_uz = 2:pde_pars.nz;
    v_uz = -pde_pars.dt*diffusion * ones(1,pde_pars.nz-1);

    r_ind_dz = 1:pde_pars.nz;
    c_ind_dz = 1:pde_pars.nz;
    v_dz(1,2:pde_pars.nz-1) = 1 + pde_pars.dt*degradation/ndims + 2*pde_pars.dt*diffusion;
    v_dz([1,pde_pars.nz]) = 1 + pde_pars.dt*degradation/ndims + pde_pars.dt*diffusion;

    r_ind_lz = 2:pde_pars.nz;
    c_ind_lz = 1:pde_pars.nz-1;
    v_lz = -pde_pars.dt*diffusion * ones(1,pde_pars.nz-1);

    pde_pars.Mz = sparse([r_ind_uz,r_ind_dz,r_ind_lz],...
        [c_ind_uz,c_ind_dz,c_ind_lz],...
        [v_uz    ,v_dz   ,v_lz    ]);
end