function [L,cells_Lind] = performMove(L,cells_Lind,cell_Lind,cell_ind,...
    TME_size,current_type,grad)

if isnan(grad) % move without bias
    [movement_ind,success_flag] = chooseNeighboringCompartment(cell_Lind,TME_size);
    if ~success_flag
        return;
    end

else % move with bias

    if length(TME_size)==2

        directions = [-1,0;...
            0,-1;...
            1,0;...
            0,1];

        weights = directions * grad;
        positive_weights_log = weights > 0;
        if sum(positive_weights_log)==1
            temp_ind = find(positive_weights_log);
        else
            weights(~positive_weights_log) = 0;
            temp_ind = randsample(1:4,1,true,weights);
        end
        switch temp_ind % handle cases where the gradient movement takes the cell off the grid
            case 1
                if mod(cell_Lind,TME_size(1))==1
                    return;
                end
            case 2
                if cell_Lind <= TME_size(1)
                    return;
                end
            case 3
                if mod(cell_Lind,TME_size(1))==0
                    return;
                end
            case 4
                if cell_Lind + TME_size(1) > prod(TME_size)
                    return;
                end
        end
        movement_ind = directions(temp_ind,:)*[1;TME_size(1)];

    elseif length(TME_size)==3
        [xind,yind,zind] = ind2sub(TME_size,cell_Lind);

        directions = [-1,0,0;...
            0,-1,0;...
            0,0,-1;...
            1,0,0;...
            0,1,0;...
            0,0,1];

        weights = directions * grad;
        positive_weights_log = weights > 0;
        if sum(positive_weights_log)==1
            temp_ind = find(positive_weights_log);
        else
            weights(~positive_weights_log) = 0;
            temp_ind = randsample(1:6,1,true,weights);
        end
        switch temp_ind % handle cases where the gradient movement takes the cell off the grid
            case 1
                if xind==1
                    return;
                end
            case 2
                if yind == 1
                    return;
                end
            case 3
                if zind == 1
                    return;
                end
            case 4
                if xind == TME_size(1)
                    return;
                end
            case 5
                if yind == TME_size(2)
                    return;
                end
            case 6
                if zind == TME_size(3)
                    return;
                end
        end
        movement_ind = directions(temp_ind,:)*[1;TME_size(1)*[1;TME_size(2)]];

    end



end

new_ind = cell_Lind + movement_ind;

if L(new_ind)==0
    L(cell_Lind) = 0;
    L(new_ind) = current_type;
    cells_Lind{current_type}(cell_ind) = new_ind;
end