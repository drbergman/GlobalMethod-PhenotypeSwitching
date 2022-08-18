function [L,num_cells,substrates,internalized_substrates,cells_Lind,counter] = ...
    performAttack(cell_Lind,TME_size,L,substrates,...
    cells_Lind,num_cells,internalized_substrates,...
    counter,target_type)

[target_relative_ind,success_flag] = chooseNeighboringCompartment(cell_Lind,TME_size);
if ~success_flag
    return;
end
target_Lind = cell_Lind + target_relative_ind;
attack_success = L(target_Lind)==target_type;

if attack_success
    L(target_Lind) = 0;
    target_ind = find(cells_Lind{target_type}==target_Lind,1);
    num_cells(target_type) = num_cells(target_type)-1;
    for si = 1:numel(internalized_substrates)
        if ~isempty(internalized_substrates{si})
            C_increase = internalized_substrates{si}(target_ind);
            if substrates(si).method=="global"
                C_ind = substrates(si).solver.regions(cells_Lind{target_type}(target_ind));
                C_increase = C_increase/substrates(si).solver.region_volumes(C_ind);
            else
                C_ind = cells_Lind{target_type}(target_ind);
            end
            substrates(si).concentration(C_ind) = substrates(si).concentration(C_ind) + C_increase;
            internalized_substrates{si}(target_ind,:) = [];
        end
    end
    
    cells_Lind{target_type}(target_ind,:) = [];
    counter.killed(target_type) = counter.killed(target_type)+1;
end