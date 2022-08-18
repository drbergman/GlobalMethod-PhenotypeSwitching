% a script to run a cohort of samples. multiple samples can be run for a single setup by simply increasing nsamps.
% To vary a single parameter with nsamps of each, create a structure with the matching parameter name(s)
% as field(s) with the desired values varying along the first dimension. Pass this structure as the second input
% into runCohort

clearvars

addpath("./custom/")
addpath("../src/")

%% setup run details
summary_flags.make_svgs = false;
summary_flags.show_patch_plots = false;
summary_flags.save_output_by_default = false;

time_runs = true;
min_parfor_num = 4;
nsamps = 1;

base_name = 'PhenotypeSwitching';

%%
input = baseParsPhenotypeSwitching();

%% tweak input parameters
input.Tend = 1;
input.n_saves = 10;

% input.TME_size = [100,100];
input.TME_size = [20,20,20];
input.oxygen_max_local_pde_dt = 1e-2; % physicell uses 0.01 minutes
input.oxygen_max_exchange_dt = 1e-3;
input.oxygen_max_internal_agent_dt = 1e-3;
input.substrate_time_advance = 6e0;
% input.oxygen_max_dt = 1e-4;

input.proliferating_rp = 2/(24*60);
input.proliferating_rm_microns = 2;
input.proliferating_rd = input.proliferating_rp * 1e-2;

input.quiescent_rm_microns = 0.2;
input.quiescent_rd = input.proliferating_rd;
input.quiescent_chemotaxis = true;

input.necrotic_rm_microns = 0.01;
input.necrotic_rd = 1e-5;

input.allow_necrosis = true;

input.proliferating_initialization_locations = "uniform";
input.proliferating_p0 = 0.1;

input.oxygen_degradation = 0.1;
input.oxygen_diffusion_microns = 1e5;

input.oxygen_pk_solver = "matrix_exponential";
input.oxygen_exchange_solver = "matrix_exponential";
input.oxygen_internal_agent_solver = "euler_direct";
% input.oxygen_global_solver = "matlab_solver";
input.oxygen_pde_solver_global_method = "euler_direct";

input.continue_without_agents = true;

input.quiescence_as_event = true;
input.quiescence_threshold = 7.5;
input.quiescence_rate_max = 0.01;

input.oxygen_use_internal_agent_ode = true;
input.proliferating_k_aer = 1e-6;
input.quiescent_k_aer = 1e-6;
input.proliferating_oxygen_secretion_rate = 5;
input.quiescent_oxygen_secretion_rate = 5;

% input.use_regions_for_diffusion = false; % if true, then solve PDE for global method and do cell reactions separately using global method
input.oxygen_global_method_in_one_ode = false;
input.oxygen_global_solver = "matlab_solver";

input.bv_location = "max_shell";
input.oxygen_regions_type = "concentric_circles";

input.use_quorum = false;
input.quorum_max_local_pde_dt = 1e-2; % physicell uses 0.01 minutes
input.quorum_pde_solver_global_method = "matrix_exponential";
input.quorum_global_method_in_one_ode = false;
input.quorum_regions_type = "rectangular_grid";

input.proliferating_quorum_export_rate = 0.01;

%% choose regions to run

input.oxygen_rectangle_compartment_size = [input.TME_size(1:end-1),1];
% input.oxygen_num_concentric_circle_regions = 50;

% input.quorum_rectangle_compartment_size = [5,5];
% input.quorum_num_concentric_circle_regions = 20;

%% prepare runs and do runs
[input,out,use_parpool] = runCohort(input,struct("oxygen_num_concentric_circle_regions",[Inf;70]),nsamps,min_parfor_num);

%% make plots
figs = plotOutputPhenotypeSwitching(out,input,nsamps,summary_flags);

%% Save data
if summary_flags.save_output_by_default
    fileFormatBase = [base_name,'%03d'];
    num = 1; % vary to avoid overwriting
    filename = sprintf(fileFormatBase,num);
    save(['data/',filename,'.mat'],'-regexp','^(?!(figs)$).','-v7.3') % do not save figures to .mat file % do not save figures to .mat file
    output_formats = {'png','svg'};
    for fi = 1:numel(figs)
        fig_name = sprintf('figs/%s_%s',filename,figs(fi).name);
        savefig(figs(fi).handle,fig_name)
        for ofi = 1:length(output_formats)
            figs(fi).handle.Position = figs(fi).print_position;
            fig_name = sprintf('figs/%s/%s_%s',output_formats{ofi},filename,figs(fi).name);
            print(figs(fi).handle,fig_name,sprintf('-d%s',output_formats{ofi}))
        end
    end
end