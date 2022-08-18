% Same as RunFigure3 (quiescence as event) but also allow quiescent cells to chemotax towards oxygen

% Also creates Figure S8.

clearvars

addpath("custom/")
addpath("../src/")

%% setup run details
summary_flags.make_svgs = false;
summary_flags.show_patch_plots = false;
summary_flags.save_output_by_default = false;
summary_flags.codensity_computation_on = false; % computing codensities takes a long time!
base_name = 'figure4_';

min_parfor_num = 20;
nsamps = 8;

input = baseParsPhenotypeSwitching();

%% tweak input parameters
input.Tend = 1;

input.TME_size = [100,100];
input.oxygen_max_local_pde_dt = 1e-2; % physicell uses 0.01 minutes
input.oxygen_max_exchange_dt = 1e-3;
input.oxygen_max_internal_agent_dt = 1e-3;
input.substrate_time_advance = 6e0;

input.track_movement = true; % will track the movement history of each cell

input.proliferating_rp = 2/(24*60);
input.proliferating_rm_microns = 2;
input.proliferating_rd = input.proliferating_rp * 1e-3;

input.quiescent_rm_microns = 0.2;
input.quiescent_rd = 0.01 * input.proliferating_rd;

input.quiescent_chemotaxis = true;
input.oxygen_diffusion_microns = 1e5;

input.oxygen_pk_solver = "matrix_exponential";
input.oxygen_exchange_solver = "matrix_exponential";
input.oxygen_pde_solver_global_method = "matrix_exponential";

input.continue_without_agents = true;

input.quiescence_as_event = true;
input.quiescence_threshold = 7.5;
input.quiescence_rate_max = 0.01;

input.oxygen_use_internal_agent_ode = true;
input.k_aer = 1e-6;
input.oxygen_secretion_rate = 1;

% input.use_regions_for_diffusion = false; % if true, then solve PDE for global method and do cell reactions separately using global method
input.oxygen_global_method_in_one_ode = false;
%% choose regions to run
M = [ones(1,length(input.TME_size));... % rectangle compartment size of [1,1] or [1,1,1] prompts the use of the local method
    input.TME_size(1:end-1),1];   % any other rectangle compartment size prompts the use of the global method

%% prepare runs and do runs
[input,out,use_parpool] = runCohort(input,struct("oxygen_rectangle_compartment_size",M),nsamps,min_parfor_num);

%% make panels
addpath("./FigureScripts/")

if ~summary_flags.codensity_computation_on
    warning('Codensities will not be computed and their plots not produced.')
end
FigurePanels(out,input,summary_flags.codensity_computation_on);
