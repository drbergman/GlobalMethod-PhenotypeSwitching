function [L,num_cells,cells_Lind] = initializeAgents(input,TME_size,ntypes)


L = zeros(TME_size,"uint8");
N = prod(TME_size);
cells_Lind = cell(ntypes,1);
num_cells = zeros(ntypes,1);


initialization_dict = containers.Map(["uniform","center","set_location"],[3,2,1]);  


for type_ind = ntypes:-1:1

    type = input.type_names(type_ind);

    if ~isfield(input,sprintf("%s_initialization_locations",type))
        agent_initialization_location(type_ind) = "uniform";
    else
        agent_initialization_location(type_ind) = input.(sprintf("%s_initialization_locations",type));
    end

    initialization_vals(type_ind) = initialization_dict(agent_initialization_location(type_ind));

end

[~,type_order] = sort(initialization_vals);

for i = 1:ntypes
   
    type = input.type_names(type_order(i));


    switch agent_initialization_location(type_order(i))
        case "uniform"
            p0 = input.(sprintf("%s_p0",type));
            
            if i>1 % compensate for some cells already having been placed
                R = find(L==0); % remaining available sites
                nR = length(R); % number of remaining sites
                p0 = p0*N/nR;
            else
                R = 1:N;
                nR = N;
            end

            num_cells(type_order(i)) = binornd(nR,p0);
            cells_Lind{type_order(i)} = reshape(R(randperm(nR,num_cells(type_order(i)))),[],1);


        case "center"
            p0 = input.(sprintf("%s_p0",type));
            for di = length(TME_size):-1:1
                coords{di} = 1:TME_size(di);
            end
            locs = allCombos(coords{:},'matlab');
            if i>1 % compensate for some cells already having been placed
                locs(L(:)~=0,:) = [];
%                 if input.ndims==2
%                     locs(L(sub2ind(TME_size,locs(:,1),locs(:,2)))~=0,:) = [];
%                 elseif input.ndims==3
%                     locs(L(sub2ind(TME_size,locs(:,1),locs(:,2),locs(:,3)))~=0,:) = [];
%                 end
            end
            dists = sum((locs-0.5*(1+TME_size)).^2,2);
            [~,order] = sort(dists,1,"ascend");
            starting_inds = order(1:round(p0*N));
            num_cells(type_order(i)) = length(starting_inds);
            if length(TME_size)==2
                cells_Lind{type_order(i)} = sub2ind(TME_size,locs(starting_inds,1),locs(starting_inds,2));
            elseif length(TME_size)==3
                cells_Lind{type_order(i)} = sub2ind(TME_size,locs(starting_inds,1),locs(starting_inds,2),locs(starting_inds,3));
            end

        case "set_location"
            all_inds = reshape(1:N,TME_size);
            xx = cell(input.ndims,1);
            [xx{:}] = ind2sub(TME_size,all_inds);
            distance_to_center = zeros(TME_size);
            disc_center = input.(sprintf("%s_disc_center",type)) / input.lattice_h + ones(1,input.ndims);
            for di = 1:input.ndims
                distance_to_center = distance_to_center + (xx{di} - disc_center(di)).^2;
            end
            cells_Lind{type_order(i)} = find(distance_to_center <= input.(sprintf("%s_disc_radius",type))^2/input.lattice_h^2);
            num_cells(type_order(i)) = length(cells_Lind{type_order(i)});

    end

    L(cells_Lind{type_order(i)}) = type_order(i);

end

