function on_lattice_rel_inds = chooseAllNeighboringCompartment(cell_Lind,TME_size)

switch length(TME_size)
    case 2
        [x(1),x(2)] = ind2sub(TME_size,cell_Lind);
        on_lattice_rel_inds = [-1,0;...
            0,-1;...
            1,0;...
            0,1];

    case 3
        [x(1),x(2),x(3)] = ind2sub(TME_size,cell_Lind);
        on_lattice_rel_inds = [-1,0,0;...
            0,-1,0;...
            0,0,-1;...
            1,0,0;...
            0,1,0;...
            0,0,1];
end

for di = 1:length(TME_size)
    switch x(di)
        case 1
            on_lattice_rel_inds(on_lattice_rel_inds(:,di)==-1,:) = [];
        case TME_size(di)
            on_lattice_rel_inds(on_lattice_rel_inds(:,di)==1,:) = [];
    end
end

on_lattice_rel_inds = on_lattice_rel_inds * cumprod([1;TME_size(1:end-1)']);