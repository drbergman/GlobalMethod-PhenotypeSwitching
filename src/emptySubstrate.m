function S = emptySubstrate()

S.name = "";
S.pars.types = []; % types that will be internalizing/exporting this substrate
S.pars.ntypes = 0; % number of types that will be internalizing/exporting this substrate
S.pars.diffusion_microns = 0;
S.pars.degradation = 0;

S.pars.is_exchange = false;
S.pars.export_rate = 0;
S.pars.secretion_rate = 0;
S.pars.uptake_rate = 0;

S.pars.concentration_initialization = "equal_to_constant";
S.pars.initial_concentration = 0;

S.pars.is_pk = false;
S.pars.circulation_concentration = 0;
S.pars.fluid_exchange_rate = 0;

S.pars.contains_dcs = false;
S.pars.dc_location_type = "none";
S.pars.dc_val = 0;
S.pars.disc_center = 0;
S.pars.disc_radius = NaN;

S.pars.use_regions_for_diffusion = true; % default to global method (rather than hybrid) when using regions
S.pars.global_method_in_one_ode = false;
S.pars.use_internal_agent_ode = false;

S.pars.global_solver = false;
S.pars.pk_solver = "euler_direct";
S.pars.pde_solver_global_method = "euler_direct";
S.pars.exchange_solver = "euler_direct";
S.pars.internal_agent_solver = "euler_direct";

S.pars.global_matlab_solver = @ode45;
S.pars.pk_matlab_solver = @ode45;
S.pars.pde_matlab_solver = @ode45;
S.pars.exchange_matlab_solver = @ode45;
S.pars.internal_agent_matlab_solver = @ode45;

S.pars.regions_type = "rectangular_grid";
S.pars.num_concentric_circle_regions = Inf;
S.pars.rectangle_compartment_size = [1,1];
S.pars.distance_delta = 0;
S.pars.blood_vessels_in_separate_region = false;

S.method = "local";

S.pars.max_dt = NaN;
S.pars.max_pk_dt = NaN;
S.pars.max_local_pde_dt = NaN;
S.pars.max_global_pde_dt = NaN;
S.pars.max_exchange_dt = NaN;
S.pars.max_internal_agent_dt = NaN;

% S.pars.dt = 1;
% 
% S.pars.pk_dt = 1;
% S.pars.pde_dt = 1;
% S.pars.exchange_dt = 1;
% S.pars.interngal_agent_dt = 1;
% 
% S.pars.num_pk_steps = 1;
% S.pars.num_pde_steps = 1;
% S.pars.num_exchange_steps = 1;
% S.pars.num_internal_agent_steps = 1;

% S.pars.M_exchange_agent = zeros(2);
% S.pars.M_exchange_agent_eV = eye(2);
% S.pars.M_exchange_agent_ev = zeros(2,1);
% S.pars.M_exchange_agent_eV_inv = eye(2);

S.t = 0;
S.setup_done = false;

S.pars.internal_agent_ode = @(x,p) zeros(size(x));