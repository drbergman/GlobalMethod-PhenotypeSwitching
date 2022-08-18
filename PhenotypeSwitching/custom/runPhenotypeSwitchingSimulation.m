function out = runPhenotypeSwitchingSimulation(input)

% code up abm with tumor cells trying to get enough oxygen to not become
% quiescent
%% load parameters from input
TME_size = input.TME_size;
save_dt = input.save_dt;
Tend = input.Tend;
BV = input.BV;
quiescence_as_event = input.quiescence_as_event;
event_pars = input.event_pars;
quiescence_threshold = input.quiescence_threshold;
substrates = input.substrates;
substrate_time_advance = input.substrate_time_advance;

ntypes = input.ntypes;
nsubstrates = input.nsubstrates;

allow_necrosis = input.allow_necrosis;
%%
t = 0;
next_substrate_t = 0;
save_t = save_dt;

%% setup probabilities
[L,num_cells,cells_Lind] = initializeAgents(input,TME_size,ntypes);

last_cell_id = 0;
for type_ind = ntypes:-1:1
    cell_ids{type_ind} = reshape(last_cell_id + (1:num_cells(type_ind)),[],1);
    last_cell_id = last_cell_id + num_cells(type_ind);
end

internalized_substrates(1:ntypes,1:nsubstrates) = {zeros(0,1)};
internalized_substrates{1,1} = zeros(num_cells(1),1);
substrates = initializeMicroenvironment(input,substrates);

if event_pars.quiescence_as_event
    quiescence_rates = event_pars.quiescence_rate(internalized_substrates{1});
else
    quiescence_rates = [];
end


%% setup outputs
out.t = 0;
out.num_cells = num_cells;
out.cell_types_per_compartment = L;
out.prolifs = zeros(ntypes,1);
out.contact_inhibitions = zeros(ntypes,1);
out.spontaneous_apoptosis = zeros(ntypes,1);

out.C = {substrates(:).concentration};
out.internalized_substrates = internalized_substrates;

if input.track_movement
    out.paths = repmat(struct("path",[]),sum(num_cells),1);
    for type_ind = 1:2
        for ci = 1:num_cells(type_ind)
%             out.paths(cell_ids{type_ind}(ci)).type = type_ind;
            out.paths(cell_ids{type_ind}(ci)).path = [type_ind;0;cells_Lind{type_ind}(ci)];
        end
    end
else
    out.paths = [];
end

counter.prolifs = zeros(ntypes,1);
counter.contact_inhibitions = zeros(ntypes,1);
counter.spontaneous_apoptosis = zeros(ntypes,1);

while t<Tend && sum(num_cells)>0

    if t>=next_substrate_t
        assert(numel(cells_Lind{1})==numel(cell_ids{1}))

        [substrates,internalized_substrates,L,cells_Lind,num_cells(1:2),cell_ids] = ...
            substrateUpdatePhenotypeSwitching(substrates,internalized_substrates,...
            min(t+substrate_time_advance,Tend),TME_size,quiescence_as_event,L,...
            cells_Lind,BV,...
            quiescence_threshold,cell_ids);

        next_substrate_t = min([substrates.t]);

        if event_pars.quiescence_as_event
            quiescence_rates = event_pars.quiescence_rate(internalized_substrates{1});
        end
        assert(numel(cells_Lind{1})==numel(cell_ids{1}))

    end

    [event,current_type,cell_ind,dt] = eventSelectionPhenotypeSwitching(event_pars,num_cells,quiescence_rates);

    cell_Lind = cells_Lind{current_type}(cell_ind);

    t = t+dt;

    if t>Tend
        t = Tend;
    else
        switch event
            case 1 % proliferation

                assert(numel(cells_Lind{1})==numel(cell_ids{1}))

                [L,cells_Lind{current_type},...
                    internalized_substrates(current_type,:),num_cells,...
                    counter,prolif_success,cell_ids{current_type},last_cell_id] = performProliferation(L,...
                    cells_Lind{current_type},...
                    internalized_substrates(current_type,:),num_cells,...
                    cell_Lind,cell_ind,TME_size,...
                    current_type,counter,cell_ids{current_type},last_cell_id,true);

                if event_pars.quiescence_as_event && current_type==1 && prolif_success
                    quiescence_rates([cell_ind,num_cells(1)]) = event_pars.quiescence_rate(internalized_substrates{1}([cell_ind,num_cells(1)]));
                end


                assert(numel(cells_Lind{1})==numel(cell_ids{1}))

                if input.track_movement && prolif_success % last_cell_id ~= length(out.paths)
                    out.paths(end+1).path = [current_type;t;cells_Lind{current_type}(end)];
                end
    
            case 2 % movement
                if input.quiescent_chemotaxis && current_type==2
                    if substrates(1).method == "local" && isempty(substrates(1).solver.regions)
                        substrates(1).solver.regions = reshape(1:numel(substrates(1).concentration),TME_size);
                    end
                    grad = computeGradient(substrates(1).concentration(substrates(1).solver.regions),cell_Lind,TME_size,input.ndims);
                    if rand() < input.grad_to_bias(grad)
                    else
                        grad = NaN;
                    end
                elseif input.use_quorum && current_type==1
                    if substrates(2).method == "local" && isempty(substrates(2).solver.regions)
                        substrates(2).solver.regions = reshape(1:prod(TME_size),TME_size);
                    end
                    grad = computeGradient(substrates(2).concentration(substrates(2).solver.regions),cell_Lind,TME_size,input.ndims);
                    if rand() < input.grad_to_bias_quorum(grad)
                        disp('')
                    else
                        grad = NaN;
                    end
                else
                    grad = NaN;
                end
                [L,cells_Lind] = performMove(L,cells_Lind,cell_Lind,...
                    cell_ind,TME_size,current_type,...
                    grad);

                if input.track_movement
                    if any(out.paths(cell_ids{current_type}(cell_ind)).path([1,3],end)~=[current_type;cells_Lind{current_type}(cell_ind)])
                        out.paths(cell_ids{current_type}(cell_ind)).path(:,end+1) = [current_type;t;cells_Lind{current_type}(cell_ind)];
                    end
                end

            case 3 % death
                [L,cells_Lind{current_type},...
                    internalized_substrates(current_type,:),num_cells,...
                    substrates,counter] = ...
                    performDeath(L,cells_Lind{current_type},...
                    internalized_substrates(current_type,:),num_cells,...
                    cell_Lind,cell_ind,current_type,...
                    substrates,counter);

                if allow_necrosis && current_type==2
                    L(cell_Lind) = 3;
                    cells_Lind{3}(end+1,1) = cell_Lind;
                    cell_ids{3}(end+1,1) = cell_ids{current_type}(cell_ind);
                    num_cells(3) = num_cells(3)+1;
                end

                if event_pars.quiescence_as_event && current_type == 1
                    quiescence_rates(cell_ind) = [];
                end

                if input.track_movement
                    out.paths(cell_ids{current_type}(cell_ind)).path(:,end+1) = [-current_type;t;cell_Lind];
                end

                cell_ids{current_type}(cell_ind) = [];

            case 4 % quiescence
                [L,cells_Lind,internalized_substrates,num_cells,cell_ids] = ...
                    performQuiescence(L,cells_Lind,internalized_substrates,...
                    num_cells,cell_Lind,cell_ind,cell_ids);
                    
                if event_pars.quiescence_as_event
                    quiescence_rates(cell_ind) = [];
                end

            otherwise
                error('no event done')
        end
    end

    if t>=save_t
        [out,save_t] = savePhenotypeSwitching(out,t,num_cells,L,counter,substrates,internalized_substrates,save_t,save_dt,Tend);
    end

end

if t<Tend && isfield(input,"continue_without_agents") && input.continue_without_agents
    
    while t<Tend
        dt = save_t-t;
        t = t+dt;
        if t>=next_substrate_t

            [substrates,internalized_substrates,L,cells_Lind,num_cells(1:2)] = ...
                substrateUpdatePhenotypeSwitching(substrates,internalized_substrates,...
                min(t+substrate_time_advance,Tend),TME_size,quiescence_as_event,L,...
                cells_Lind,BV,...
                quiescence_threshold);

            next_substrate_t = min([substrates.t]);
        end
        
        if t>=save_t
            [out,save_t] = savePhenotypeSwitching(out,t,num_cells,L,counter,substrates,internalized_substrates,save_t,save_dt,Tend);
        end
    end

end

out.num_cells = out.num_cells';
out.cell_types_per_compartment = permute(out.cell_types_per_compartment,[ndims(L)+1,1:ndims(L)]); % so time varies along first dimension
out.prolifs = out.prolifs';
out.contact_inhibitions = out.contact_inhibitions';
out.spontaneous_apoptosis = out.spontaneous_apoptosis';

for si = 1:nsubstrates
    out.C{si} = permute(out.C{si},[ndims(out.C{si}),1:(ndims(out.C{si})-1)]); % so time varies along first dimension
end
out.internalized_substrates = permute(out.internalized_substrates,[ndims(out.internalized_substrates),1:(ndims(out.internalized_substrates)-1)]); % so time varies along first dimension

if input.track_movement
    for type_ind = 1:ntypes
        for ci = 1:num_cells(type_ind)
            out.paths(cell_ids{type_ind}(ci)).path(:,end+1) = [type_ind;Tend;cells_Lind{type_ind}(ci)];
        end
    end
end
