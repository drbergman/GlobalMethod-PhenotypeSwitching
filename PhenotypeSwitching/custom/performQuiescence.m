function [L,cells_Lind,internalized_substrates,num_cells,cell_ids] = ...
    performQuiescence(L,cells_Lind,internalized_substrates,num_cells,...
    cell_Lind,cell_ind,cell_ids)

L(cell_Lind) = 2;

cells_Lind{2} = [cells_Lind{2};cell_Lind];
cells_Lind{1}(cell_ind,:) = [];

cell_ids{2} = [cell_ids{2};cell_ids{1}(cell_ind)];
cell_ids{1}(cell_ind,:) = [];

internalized_substrates{2} = [internalized_substrates{2};internalized_substrates{1}(cell_ind)];
internalized_substrates{1}(cell_ind,:) = [];

num_cells(1:2) = num_cells(1:2)+[-1;1];