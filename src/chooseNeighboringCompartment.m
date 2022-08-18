function [relative_ind,success_flag] = chooseNeighboringCompartment(cell_Lind,TME_size)

g = rand();
relative_ind = 0;
if length(TME_size)==2
    if g < 0.5 % then choose left or right
        relative_ind = 2*randi([0,1])-1;
        switch mod(cell_Lind,TME_size(1))
            case 0 % then cell is at right boundary of compartment grid
                if relative_ind == 1 || TME_size(1)==1
                    success_flag = false;
                    return; % attempted to attack off compartment grid in x direction
                    % Also, if only one compartment in x direction, leaving this way means nothing happened
                end
            case 1 % then cell is at left boundary of compartment grid
                if relative_ind == -1
                    success_flag = false;
                    return; % attempted to attack off compartment grid in x direction
                end
        end
    else % then choose up or down with equal probability
        relative_ind = TME_size(1)*(2*randi([0,1])-1);
        if cell_Lind + relative_ind < 1 || cell_Lind + relative_ind >= prod(TME_size)
            success_flag = false;
            return; % attempted to attack off compartment grid in y direction
        end
    end
elseif length(TME_size)==3
    [xind,yind,zind] = ind2sub(TME_size,cell_Lind);
    if g < 0.3333333333333333 % then attack left or right
        dx = 2*randi([0,1])-1;
        if any(xind+dx==[0,TME_size(1)+1]) % attacking off the grid
            success_flag = false;
            return; % attempted to attack off compartment grid in x direction
            % Also, if only one compartment in x direction, leaving this way means nothing happened
        else
            relative_ind = dx;
        end
    elseif g < 0.6666666666666666 % then attack up or down with equal probability
        dy = 2*randi([0,1])-1;
        if any(yind+dy==[0,TME_size(2)+1]) % attacking off the grid
            success_flag = false;
            return; % attempted to attack off compartment grid in y direction
        else
            relative_ind = dy*TME_size(1);
        end
    else
        dz = 2*randi([0,1])-1;
        if any(zind+dz==[0,TME_size(3)+1]) % attacking off the grid
            success_flag = false;
            return; % attempted to attack off compartment grid in z direction
        else
            relative_ind = dz*TME_size(1)*TME_size(2);
        end
    end
end
success_flag = true;