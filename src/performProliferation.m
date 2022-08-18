function [L,cells_Lind,internalized_substrates,num_cells,counter,success_flag,cell_ids,last_cell_id] = ...
    performProliferation(L,cells_Lind,internalized_substrates,num_cells,...
    cell_Lind,cell_ind,TME_size,current_type,...
    counter,cell_ids,last_cell_id,try_one_direction)

if try_one_direction
    [movement_ind,success_flag] = chooseNeighboringCompartment(cell_Lind,TME_size);
    if ~success_flag
        return;
    end
    new_ind = cell_Lind + movement_ind;

    success_flag = L(new_ind)==0;
else
    on_lattice_rel_inds = chooseAllNeighboringCompartment(cell_Lind,TME_size);
    on_lattice_rel_inds = on_lattice_rel_inds(L(cell_Lind + on_lattice_rel_inds)==0);
    if isempty(on_lattice_rel_inds)
        success_flag = false;
    elseif length(on_lattice_rel_inds)==1
        new_ind = cell_Lind + on_lattice_rel_inds;
        success_flag = true;
    else
        new_ind = cell_Lind + on_lattice_rel_inds(randi(length(on_lattice_rel_inds)));
        success_flag = true;
    end
end

if success_flag
    L(new_ind) = current_type;
    cells_Lind = [cells_Lind;new_ind];

    if ~isempty(cell_ids)
        last_cell_id = last_cell_id + 1;
        cell_ids(end+1) = last_cell_id; % doing it in this order saves me adding one twice
    end

    num_cells(current_type) = num_cells(current_type)+1;
    for si = 1:numel(internalized_substrates)
        if ~isempty(internalized_substrates{si})
            internalized_substrates{si}([cell_ind,end+1],1) = internalized_substrates{si}(cell_ind)/2;
        end
    end
    counter.prolifs(current_type) = counter.prolifs(current_type)+1;
else
    counter.contact_inhibitions(current_type) = counter.contact_inhibitions(current_type)+1;
end
