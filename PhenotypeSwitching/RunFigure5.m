% This script runs both the local and global methods until Day 15. Quiescent cells chemotax
% towards oxygen; proliferating cells chemotax towards quorum factor. Simulation occurs in 3D 
% with blood vessels surrounding the spheroid

clearvars

addpath("custom/")
addpath("../src/")

%% setup run details
summary_flags.make_svgs = false;
summary_flags.show_patch_plots = false;
summary_flags.save_output_by_default = false;
summary_flags.codensity_computation_on = false; % computing codensities takes a long time!
base_name = 'figure5_';

min_parfor_num = 20;
nsamps = 8;

input = baseParsPhenotypeSwitching();

%% tweak input parameters
input.Tend = 15*1440;

input.TME_size = [20,20,20];
input.oxygen_max_local_pde_dt = 1e-2; % physicell uses 0.01 minutes
input.oxygen_max_exchange_dt = 1e-3;
input.oxygen_max_internal_agent_dt = 1e-3;
input.substrate_time_advance = 6e0;

input.proliferating_rp = 2/(24*60);
input.proliferating_rm_microns = 2;
input.proliferating_rd = input.proliferating_rp * 5e-2;

input.quiescent_rm_microns = 0.2;
input.quiescent_rd = 1.0 * input.proliferating_rd;

input.proliferating_initialization_locations = "center";
input.proliferating_p0 = 0.25;

input.quiescent_chemotaxis = true;
input.allow_necrosis = false;

input.oxygen_degradation = 0.1;
input.oxygen_diffusion_microns = 1e5;

input.oxygen_pk_solver = "matrix_exponential";
input.oxygen_exchange_solver = "matrix_exponential";
input.oxygen_pde_solver_global_method = "matrix_exponential";

input.continue_without_agents = true;

input.quiescence_as_event = true;
input.quiescence_threshold = 10;
input.quiescence_rate_max = 0.01;

input.oxygen_use_internal_agent_ode = true;
input.proliferating_k_aer = 1e-6;
input.quiescent_k_aer = 1e-6;
input.proliferating_oxygen_secretion_rate = 1;
input.quiescent_oxygen_secretion_rate = 1;
input.oxygen_circulation_concentration = 38; 
input.oxygen_fluid_exchange_rate = 12.5; % reduced because there are more blood vessels in this simulation due to their arrangement; otherwise the entire TME is flooded with oxygen

input.oxygen_global_method_in_one_ode = false;

input.oxygen_regions_type = "concentric_circles"; % concentric circles match the geometry of the blood vessels and the spheroid
input.bv_location = "weighted_shell";
input.shell_radius = 190;

input.use_quorum = true;
input.quorum_max_local_pde_dt = 1e-2; % physicell uses 0.01 minutes
input.quorum_pde_solver_global_method = "matrix_exponential";
input.quorum_global_method_in_one_ode = false;
input.quorum_regions_type = "concentric_circles";
input.quorum_num_concentric_circle_regions = 40;

%% choose regions to run
M = [Inf;70]; % When using concentric circles for regions, using Inf prompts using the local method; a finite value attempts to use that many regions with equal spacing between

%% prepare runs and do runs
[input,out,use_parpool] = runCohort(input,struct('oxygen_num_concentric_circle_regions',M),nsamps,min_parfor_num);

%% make panels
addpath("./FigureScripts/")

if ~summary_flags.codensity_computation_on
    warning('Codensities will not be computed and their plots not produced.')
end
FigurePanels(out,input,summary_flags.codensity_computation_on);
