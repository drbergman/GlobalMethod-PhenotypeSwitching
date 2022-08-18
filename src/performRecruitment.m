function [L,cells_Lind,num_cells,counter] = performRecruitment(BV_sum,BV_cumsum,L,current_type,cells_Lind,num_cells,counter)

new_ind = find(BV_sum*rand()<=BV_cumsum,1); % try putting new cell here
if L(new_ind) == 0
    L(new_ind) = current_type;
    cells_Lind{current_type} = [cells_Lind{current_type};new_ind];
    num_cells(current_type) = num_cells(current_type)+1;

    counter.recruits(current_type) = counter.recruits(current_type)+1;
end