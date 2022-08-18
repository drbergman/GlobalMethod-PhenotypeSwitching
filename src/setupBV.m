function BV = setupBV(input)

bv_location = input.bv_location;
TME_size = input.TME_size;

switch bv_location

    case "floor"
        BV = false(TME_size);
        if length(TME_size)==2
            BV(:,1) = true;
        elseif length(TME_size)==3
            BV(:,:,1) = true;
        else
            error('not ready for a TME with these dimensions')
        end

    case "outside"
        BV = true(TME_size);
        if length(TME_size)==2
            BV(2:end-1,2:end-1) = false;
        elseif length(TME_size)==3
            BV(2:end-1,2:end-1,2:end-1) = false;
        end

    case "lines"
        BV = false(TME_size);
        ymid = round((TME_size(2)+1)/2);
        dely = 10; % 10 cell widths = 200 microns
        yinds = ymid + [0:-dely:-TME_size(2),dely:dely:TME_size(2)];
        yinds(yinds<1 | yinds>TME_size(2)) = [];
        if length(TME_size)==2
            BV(:,yinds) = true;
        else
            zmid = round((TME_size(3)+1)/2);
            delz = 10; % 10 cell widths = 200 microns
            zinds = zmid + [0:-delz:-TME_size(3),delz:delz:TME_size(3)];
            zinds(zinds<1 | zinds>TME_size(3)) = [];
            BV(:,yinds,zinds) = true;
        end

    case "disc"
        BV = false(TME_size);
        all_inds = reshape(1:prod(TME_size),TME_size);
        xx = cell(length(TME_size),1);
        [xx{:}] = ind2sub(TME_size,all_inds);
        distance_to_center = zeros(TME_size);
        for i = 1:length(TME_size)
            distance_to_center = distance_to_center + ((xx{i}-1)*input.lattice_h - input.bv_disc_center(i)).^2; % distance in microns
        end
        BV(distance_to_center <= input.bv_disc_radius^2) = true;

    case "max_shell"
        BV = false(TME_size);
        ndims = length(TME_size);
        radius = 0.5*min(TME_size)-0.5;
        if ndims==2
            [X,Y] = ndgrid(1:TME_size(1),1:TME_size(2));
            P = [X(:),Y(:)];
        else
            [X,Y,Z] = ndgrid(1:TME_size(1),1:TME_size(2),1:TME_size(3));
            P = [X(:),Y(:),Z(:)];
        end
        D = sqrt(sum((P - 0.5*(1+TME_size)).^2,2)); % distance of each lattice point from center of TME
        interior = D<=radius;
        P = P(interior,:);
        if ndims==2
            Pind = boundary(P(:,1),P(:,2));
            ind = sub2ind(TME_size,P(Pind,1),P(Pind,2));
        else
            Pind = boundary(P(:,1),P(:,2),P(:,3));
            ind = sub2ind(TME_size,P(Pind,1),P(Pind,2),P(Pind,3));
        end
        BV(ind) = true;

    case "weighted_shell"
        ntheta = 1000;
        theta = linspace(0,2*pi,ntheta)';
        dtheta = theta(2)-theta(1);
        r = input.shell_radius / input.lattice_h;
        if length(TME_size)==2

            X{1} = r * cos(theta) + 0.5*(1+TME_size(1));
            X{2} = r * sin(theta) + 0.5*(1+TME_size(2));

            weights = r*dtheta;
            
        else
            nphi = 1000;

            phi = linspace(0,pi,nphi);
            dphi = phi(2)-phi(1);

            T = zeros(ntheta*nphi,1);
            P = zeros(ntheta*nphi,1);
            [T(:),P(:)] = ndgrid(theta,phi);

            X{1} = r * cos(T) .* sin(P) + 0.5*(1+TME_size(1));
            X{2} = r * sin(T) .* sin(P) + 0.5*(1+TME_size(2));
            X{3} = r * cos(P) + 0.5*(1+TME_size(3));

            weights = r^2*sin(P)*dtheta*dphi; % dA in spherical coordinates
            
        end
        [coord,~,coord_ind] = unique(round(cat(2,X{:})),'rows');
        BV = reshape(...
                accumarray(...
                    1 + (coord(coord_ind,:)-1) * [1,cumprod(TME_size(1:end-1))]',... % index of integral coordinate each point rounds to
                    weights,... % weight associated with each point
                    [prod(TME_size),1]),... % put it in array with correct number of elements
                TME_size); % reshape to match TME size

    case "A"
        BV = false(TME_size);
        if length(TME_size)==3
            error("not sure what the shape of an 'A' in 3D should be")
        end
        crossbar_height = round(.6*TME_size(2)); % height of crossbar of A
        BV(:,crossbar_height) = true;
        
        left_leg_x_vals = linspace(round(TME_size(1)/3),round(TME_size(1)/2),1000);
        left_leg_y_vals = linspace(1,TME_size(2),1000);
        subs_left = unique(round([left_leg_x_vals',left_leg_y_vals']),'rows');
        inds_left = sub2ind(TME_size,subs_left(:,1),subs_left(:,2));
        right_leg_x_vals = linspace(round(TME_size(1)/2),round(2*TME_size(1)/3),1000);
        right_leg_y_vals = flip(linspace(1,TME_size(2),1000));
        subs_right = unique(round([right_leg_x_vals',right_leg_y_vals']),'rows');
        inds_right = sub2ind(TME_size,subs_right(:,1),subs_right(:,2));
        BV([inds_left;inds_right]) = true;

    case "none"
        BV = false(TME_size);

    otherwise
        error("%s is not an allowable way to initialize blood vessels",bv_location)
end