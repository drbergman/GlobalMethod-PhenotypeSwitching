function input = baseParsPhenotypeSwitching()

input.type = 'PhenotypeSwitching';
input.TME_size = [20,100];
input.Tend = 2; % censor date in minutes
input.substrate_time_advance = 1e-1; % how far into future to simulate substrate level dynamics (in minutes)
input.allow_necrosis = false; % whether to fix necrotic cells in place permanently after quiescent cell death or not
% input.agent_initialization_locations = "uniform"; % scatter agents uniformly across microenvironment to start simulation
input.track_movement = false; % whether or not to track movement paths for all cells
input.bv_location = "floor";

input.parameter_setup_fn = @finishParameterSetupPhenotypeSwitching;
input.sim_function = @runPhenotypeSwitchingSimulation;

input.type_names = ["proliferating","quiescent"];
input.substrate_names = "oxygen";

input.lattice_h = 20; % distance between lattice points in microns

input.n_saves = 500;
input.continue_without_agents = false;

input.proliferating_p0 = .1; % start with tumors in this proportion of spots
input.proliferating_rp = 1; % proliferation rate in per minute
input.proliferating_rm_microns = 20; % movement rate in microns per minute
input.proliferating_rd = 0.01; % death rate in per minute

input.quiescent_p0 = 0; % start with tumors in this proportion of spots
input.quiescent_rp = 0;
input.quiescent_rm_microns = 0.1;
input.quiescent_rd = 0.01;
input.quiescent_chemotaxis = false; % whether to let quiescent cells chemotax towards high oxygen

input.quiescence_threshold = 5;

input.necrotic_p0 = 0; % start with tumors in this proportion of spots
input.necrotic_rp = 0;
input.necrotic_rm_microns = 0;
input.necrotic_rd = 0.001;

input.proliferating_oxygen_uptake_rate = 10; % in per minute
input.quiescent_oxygen_uptake_rate = 10; % in per minute
input.proliferating_oxygen_secretion_rate = 10; % assume oxygen just does passive diffusion to enter/exit cells (in per minute)
input.quiescent_oxygen_secretion_rate = 10; % assume oxygen just does passive diffusion to enter/exit cells (in per minute)

%% oxygen
input.oxygen_types = ["proliferating","quiescent"];
input.oxygen_regions_type = "rectangular_grid";
input.oxygen_num_concentric_circle_regions = Inf; % if using concentric circles for regions, this runs the local method
input.oxygen_rectangle_compartment_size = [1,1]; % if using rectangular compartments for regions, this runs the local method
input.oxygen_is_pk = true;
input.oxygen_global_method_in_one_ode = true; % whether to solve the global method with a single ODE or by splitting it

input.oxygen_max_dt = 0.01; % max timestep allowed for steeping through the various DEs
input.oxygen_max_pk_dt = 1e-2; % max pk time step (in minutes); NaN means this can be as large as the substrate_time_advance
input.oxygen_max_local_pde_dt = 1e-3; % max pde time step (in minutes)
input.oxygen_max_global_pde_dt = 1e-4; % max pde time step (in minutes)
input.oxygen_max_exchange_dt = 1e-3; % max time step for solving linear exchange equation of substrates between cells and environment (in minutes)
input.oxygen_max_internal_agent_dt = 1e-3; % max time step for solving internal agent ODEs (in minutes)

input.oxygen_degradation = 0.1; % in per minute
input.oxygen_diffusion_microns = 1e5; % in um^2/minute
input.oxygen_fluid_exchange_rate = 20; % rate at which fluid is exchanged between TME and circulation (measured in per minute)
input.oxygen_circulation_concentration = 38; % multiply by fluid_exchange_rate to determine rate of amount of substrate accumulating in TME

input.oxygen_concentration_initialization = "equal_to_constant";
input.oxygen_blood_vessels_in_separate_region = true;

input.oxygen_initial_concentration = 38;

input.oxygen_global_solver = "euler_direct";
input.oxygen_pk_solver = "euler_direct";
input.oxygen_pde_solver_global_method = "euler_direct"; % how to solve the PDE in the global method when the equations are split
input.oxygen_exchange_solver = "euler_direct";
input.oxygen_internal_agent_solver = "euler_direct";

input.oxygen_internal_agent_matlab_solver = @ode45;
input.oxygen_pk_matlab_solver = @ode45;
input.oxygen_pde_matlab_solver = @ode45;
input.oxygen_exchange_matlab_solver = @ode45;

input.oxygen_use_regions_for_diffusion = true; % if true, then solve PDE for global method and do cell reactions separately using global method

%% quiescence as event parameters
input.quiescence_as_event = false; % whether or not to treat quiescence as event or switch type just based on internalized oxygen

input.quiescence_rate_function = "ramp_up"; % this means the function will be a ramp-up function (though plotted it will look like a ramp down as more oxygen means lower quiescence rate)
input.quiescence_saturation = 2.5; % internalized substrate concentration below which the quiescence rate is at the maximal level
% will use the quiescence threshold as the default level at which quiescence rate becomes positive
input.quiescence_rate_max = 0.1; % rate of quiescence (in per minute) when internalized substrate is below oxygen_saturation; at 0.1 this means that it will take 10 minutes, on average

%% use aerobic respiration inside agents to metabolize oxygen
input.oxygen_use_internal_agent_ode = false; % whether or not to use this
% input.k_aer = 0.01; % rate of aerobic respiration reaction; pulled from Toy_Metabolic_Model.xml from PhysiCell's ode-energy intracellular sample project
input.proliferating_k_aer = 1e-5; % rate of aerobic respiration reaction; according to the slide (see ref in glucose par), oxygen levels go very low inside the cell, much lower than the threshold values we are working with
input.quiescent_k_aer = 1e-5; % rate of aerobic respiration reaction; according to the slide (see ref in glucose par), oxygen levels go very low inside the cell, much lower than the threshold values we are working with
input.glucose = 10; % level of glucose in every cell; assumed constant; picked based on the steady-state from Furkan's simulation from 2021 Physicell Hackathon (https://github.com/physicell-training/ws2021/blob/main/pdfs/PhysiCell_ws2021_Session11.pdf)

%% biased movement
input.grad_to_bias_ec50 = 5e-4; % ec50 for the norm of the gradient affecting biased movement probability (in mmHg / micron)

%% quorum factor
input.use_quorum = false; % whether or not to use quorum factor

input.quorum_types = ["proliferating","quiescent"];
input.quorum_num_concentric_circle_regions = Inf; % if using concentric circles for regions, this runs the local method
input.quorum_rectangle_compartment_size = [1,1]; % if using rectangular compartments for regions, this runs the local method
input.quorum_global_method_in_one_ode = true; % whether to solve the global method with a single ODE or by splitting it
input.quorum_regions_type = "rectangular_grid";

input.quorum_max_dt = 0.01; % max timestep allowed for stepping through the various DEs
input.quorum_max_pk_dt = 1e-2; % max pk time step (in minutes); NaN means this can be as large as the substrate_time_advance
input.quorum_max_local_pde_dt = 1e-3; % max pde time step (in minutes)
input.quorum_max_global_pde_dt = 1e-4; % max pde time step (in minutes)
input.quorum_max_exchange_dt = 1e-3; % max time step for solving linear exchange equation of substrates between cells and environment (in minutes)
input.quorum_max_internal_agent_dt = 1e-3; % max time step for solving internal agent ODEs (in minutes)

input.proliferating_quorum_export_rate = 1; % export rate of quorum factor by all proliferating cells
input.quiescent_quorum_export_rate = 1; % export rate of quorum factor by all quiescent cells
input.quorum_degradation = 1; % degradation rate of quorum factor
input.quorum_diffusion_microns = 1e4; % diffusion of quorum factor in um^2/minute
input.grad_to_bias_ec50_quorum = 1e-4; % ec50 for the norm of the quorum gradient affecting biased movement probability



