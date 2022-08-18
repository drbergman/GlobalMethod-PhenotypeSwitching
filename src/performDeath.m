function [L,cells_Lind,internalized_substrates,num_cells,substrates,counter] = ...
    performDeath(L,cells_Lind,internalized_substrates,num_cells,cell_Lind,...
    cell_ind,current_type,substrates,counter)

L(cell_Lind) = 0;
num_cells(current_type) = num_cells(current_type)-1;

for si = 1:numel(internalized_substrates)
    if ~isempty(internalized_substrates{si})
        concentration_increase = internalized_substrates{si}(cell_ind);
        if substrates(si).method=="global"
            concentration_ind = substrates(si).solver.regions(cell_Lind);
            concentration_increase = concentration_increase/substrates(si).solver.region_volumes(concentration_ind); % spread the increase across the entire region 
        else
            concentration_ind = cell_Lind;
        end
        substrates(si).concentration(concentration_ind) = substrates(si).concentration(concentration_ind) + concentration_increase;
%         
%         if substrates(si).method=="global"
%             regions = substrates(si).solver.regions;
%             region_volumes = substrates(si).solver.region_volumes;
%             substrates(si).concentration(regions(cell_Lind)) = substrates(si).concentration(regions(cell_Lind)) + internalized_substrates{si}(cell_ind)/region_volumes(regions(cell_Lind));
%         else
%             substrates(si).concentration(cell_Lind) = substrates(si).concentration(cell_Lind) + internalized_substrates{si}(cell_ind);
%         end
        internalized_substrates{si}(cell_ind,:) = [];
    end
end

cells_Lind(cell_ind,:) = [];

counter.spontaneous_apoptosis(current_type) = counter.spontaneous_apoptosis(current_type)+1;