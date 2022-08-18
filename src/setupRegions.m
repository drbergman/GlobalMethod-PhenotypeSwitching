function regions = setupRegions(S,TME_size,bv_location,BV)


switch S.pars.regions_type
    case "concentric_circles"
        if length(TME_size)==2
            dist_from_origin = sqrt((((1:TME_size(1)) - (TME_size(1)+1)/2)').^2 + (((1:TME_size(2)) - (TME_size(2)+1)/2)).^2);
        elseif length(TME_size)==3
            dist_from_origin = sqrt((((1:TME_size(1)) - (TME_size(1)+1)/2)').^2 + (((1:TME_size(2)) - (TME_size(2)+1)/2)).^2 + ((reshape(1:TME_size(3),1,1,[]) - (TME_size(3)+1)/2)).^2);
        end
        if isfield(S.pars,"radius_transform_fn")
            dist_from_origin = S.pars.radius_transform_fn(dist_from_origin);
        end
        norm_dist_from_origin = dist_from_origin/max(dist_from_origin,[],'all');
        regions = ceil(S.pars.num_concentric_circle_regions*norm_dist_from_origin);

    case "rectangular_grid"
        assert(all(mod(TME_size,S.pars.rectangle_compartment_size)==0))
        x_additions = repelem((1:(TME_size(1)/S.pars.rectangle_compartment_size(1))),S.pars.rectangle_compartment_size(1));
        y_additions = repelem(((1:(TME_size(2)/S.pars.rectangle_compartment_size(2)))-1)*(TME_size(1)/S.pars.rectangle_compartment_size(1)),S.pars.rectangle_compartment_size(2));
        if length(TME_size)==2
            regions = x_additions' + y_additions;
        elseif length(TME_size)==3
            z_additions = repelem(((1:(TME_size(3)/S.pars.rectangle_compartment_size(3)))-1)*(prod(TME_size(1:2))/prod(S.pars.rectangle_compartment_size(1:2))),S.pars.rectangle_compartment_size(3));
            regions = x_additions' + y_additions + reshape(z_additions,1,1,[]);
        end

    case "distance_to_dc"
        all_inds = reshape(1:prod(TME_size),TME_size);
        xx = cell(length(TME_size),1);
        [xx{:}] = ind2sub(TME_size,all_inds);

        yy = cell(length(TME_size),1);
        [yy{:}] = ind2sub(TME_size,S.pars.dc_inds);

        distance_to_dc = zeros(prod(TME_size),numel(S.pars.dc_inds));
        for i = 1:length(TME_size)
            distance_to_dc = distance_to_dc + (xx{i}(:) - yy{i}(:)').^2;
        end
        distance_to_dc = sqrt(min(distance_to_dc,[],2));
        regions = reshape(ceil(distance_to_dc/S.pars.distance_delta),TME_size);


    case "concentric_rectangles"
        regions = zeros(TME_size);
        nregions = ceil(min(TME_size)/2);
        for i = 2:nregions
            if length(TME_size)==2
                regions(i:(end-i+1),i:(end-i+1)) = regions(i:(end-i+1),i:(end-i+1)) - 1;
            else
                regions(i:(end-i+1),i:(end-i+1),i:(end-i+1)) = regions(i:(end-i+1),i:(end-i+1),i:(end-i+1)) - 1;
            end
        end

    otherwise
        regions = reshape(1:prod(TME_size),TME_size);

end

if S.method~="local" && isfield(S.pars,"blood_vessels_in_separate_region") && S.pars.blood_vessels_in_separate_region && S.pars.is_pk

%     BV = setupBV(bv_location,TME_size);

    if S.pars.regions_type=="rectangular_grid" && (bv_location=="floor" || bv_location=="none")
        regions(BV>0) = regions(BV>0)-0.5; % keep them neighboring their original region; could perhaps choose it as +0.5 to keep it tridiagonal? not sure about this
    elseif (any(S.pars.regions_type==["concentric_circles","concentric_rectangles"]) && bv_location=="outside") || (S.pars.regions_type=="concentric_circles" && any(bv_location==["max_shell","weighted_shell"]))
        regions(BV>0) = regions(BV>0)+0.5; % move BV regions further from center
    elseif S.pars.regions_type=="distance_to_dc" && bv_location=="none"
        % should not need to fix this because there are no blood vessels to
        % cause adjustments to regions
    else
        error("Have not yet decided how to handle creating separate regions when blood vessels are at %s and regions are in %s",bv_location,S.pars.regions_type)
    end
%     [~,~,new_regions] = unique(regions);
%     regions = reshape(new_regions,size(regions));

end


% make sure the regions are numbered 1:nregions without skipping any
% numbers so that the region number is the index in the concentration
% vector
[~,~,new_regions] = unique(regions);
regions = reshape(new_regions,size(regions));


%% now take care of DCs if they exist and reset region numbers if need be
if S.method~="local" && S.pars.contains_dcs
    
    regions(S.pars.dc_inds) = floor(median(regions(S.pars.dc_inds))) + 0.5; % set the region index for the dcs to be near the index they originally were
    [~,~,new_regions] = unique(regions);
    regions = reshape(new_regions,size(regions));
end