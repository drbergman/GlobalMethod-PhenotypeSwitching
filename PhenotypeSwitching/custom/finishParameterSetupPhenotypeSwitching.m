function input = finishParameterSetupPhenotypeSwitching(input)

input.ndims = length(input.TME_size);

if input.allow_necrosis
    input.ntypes = 3;
    input.type_names(3) = "necrotic";
else
    input.ntypes = 2;
end

for type_ind = 1:input.ntypes
    input.event_pars.rp(type_ind,1) = input.(sprintf("%s_rp",input.type_names(type_ind)));
    input.event_pars.rm(type_ind,1) = input.(sprintf("%s_rm_microns",input.type_names(type_ind))) / input.lattice_h;
    input.event_pars.rd(type_ind,1) = input.(sprintf("%s_rd",input.type_names(type_ind)));
end

for type_ind_substrate = 1:numel(input.oxygen_types)
    input.oxygen_uptake_rate(type_ind_substrate) = input.(sprintf("%s_oxygen_uptake_rate",input.oxygen_types(type_ind_substrate)));
    input.oxygen_secretion_rate(type_ind_substrate) = input.(sprintf("%s_oxygen_secretion_rate",input.oxygen_types(type_ind_substrate)));
end

input.type_dict = containers.Map(input.type_names,1:input.ntypes);

input.save_dt = input.Tend/input.n_saves;
 
%% use aerobic respiration inside agents to metabolize oxygen
if input.oxygen_use_internal_agent_ode
    for type_ind_substrate = numel(input.oxygen_types):-1:1
        temp(type_ind_substrate) = -6*input.(sprintf("%s_k_aer",input.oxygen_types(type_ind_substrate)))*input.glucose;
    end
%     temp = -6*input.k_aer*input.glucose;
    input.oxygen_internal_agent_ode = @(x,i) reshape(temp(i).*reshape(x,[],length(i)).^6,[],1);
end

%% biased movement
bias_temp = input.lattice_h * input.grad_to_bias_ec50;
input.grad_to_bias = @(grad) 1 - 1./(1+sqrt(sum(grad.^2,1))./bias_temp);


%% quorum factor
if input.use_quorum
    input.substrate_names(2) = "quorum";
    input.nsubstrates = 2;
    quorum_bias_temp = input.lattice_h * input.grad_to_bias_ec50_quorum;
    input.grad_to_bias_quorum = @(grad) 1 - 1./(1+sqrt(sum(grad.^2,1))./quorum_bias_temp);

    for type_ind_substrate = 1:numel(input.quorum_types)
        input.quorum_export_rate(type_ind_substrate) = input.(sprintf("%s_quorum_export_rate",input.quorum_types(type_ind_substrate)));
    end
else
    input.nsubstrates = 1;
end

%% quiescence as event parameters
input.event_pars.quiescence_as_event = input.quiescence_as_event;
input.quiescence_as_event = input.quiescence_as_event;

if input.event_pars.quiescence_as_event
    switch input.quiescence_rate_function
        case "ramp_up"
            input.event_pars.quiescence_rate = @(is) quiescenceRateRampUp(is,input.quiescence_rate_max,input.quiescence_saturation,input.quiescence_threshold);
        otherwise
            error("have not created a quiescence rate function based on %s",input.quiescence_rate_function)
    end
end


