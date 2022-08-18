% This script works the same as RunFigure2.m except all internal agent ODEs will be solved with ode45

clearvars

addpath("custom/")
addpath("../src/")

%% setup run details
summary_flags.make_svgs = false;
summary_flags.show_patch_plots = false;
summary_flags.save_output_by_default = false;
base_name = 'figureS6_';

min_parfor_num = 4;
nsamps = 1;

input = baseParsPhenotypeSwitching();

%% tweak input parameters
input.Tend = 2;

input.TME_size = [100,100];
input.oxygen_max_local_pde_dt = 1e-2; % physicell uses 0.01 minutes
input.oxygen_max_exchange_dt = 1e-3;
input.oxygen_max_internal_agent_dt = 1e-3;
input.substrate_time_advance = 6e0;

input.proliferating_rp = 2/(24*60);
input.proliferating_rm_microns = 2;
input.proliferating_rd = input.proliferating_rp * 1e-3;

input.quiescent_rm_microns = 0.2;
input.quiescent_rd = 0.01 * input.proliferating_rd;

input.oxygen_pk_solver = "matrix_exponential";
input.oxygen_pde_solver_global_method = "matrix_exponential";
input.oxygen_exchange_solver = "matrix_exponential";
input.oxygen_internal_agent_solver = "matlab_solver";

input.continue_without_agents = true;

input.oxygen_use_internal_agent_ode = true;
input.k_aer = 1e-6;
input.oxygen_secretion_rate = 1;

input.oxygen_global_method_in_one_ode = false;
%% choose regions to run
M = [ones(1,length(input.TME_size));... % rectangle compartment size of [1,1] or [1,1,1] prompts the use of the local method
    input.TME_size(1:end-1),1];   % any other rectangle compartment size prompts the use of the global method

%% prepare runs and do runs
[input,out,use_parpool] = runCohort(input,struct("oxygen_rectangle_compartment_size",M),nsamps,min_parfor_num);

%% make panels
wall_times = arrayify(out,'wall_time');

figure;
ax = gca;
hold on;
histogram(wall_times(1,:))
histogram(wall_times(2,:))

legend({'Local','Global'})
xlabel('Wall Time (s)')
ax.FontSize = 20;
