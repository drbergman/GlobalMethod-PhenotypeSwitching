function substrates = initializeMicroenvironment(input,substrates)

for si = 1:length(substrates)

    if substrates(si).method=="global"
        substrates(si).concentration = zeros(numel(unique(substrates(si).solver.regions)),1);
    else
        substrates(si).concentration = zeros(input.TME_size);
    end
    if isfield(substrates(si).pars,"concentration_initialization")

        switch substrates(si).pars.concentration_initialization

            case "equal_to_circulation"

                substrates(si).concentration = substrates(si).concentration + substrates(si).pars.circulation_concentration;

            case "equal_to_constant"

                substrates(si).concentration = substrates(si).concentration + substrates(si).pars.initial_concentration;

            case "linear_interpolation_bottom_up"

                C_temp = zeros(input.TME_size);

                if length(input.TME_size)==2
                    C_temp = C_temp + substrates(si).pars.circulation_concentration*linspace(1,0,input.TME_size(2));
                elseif length(input.TME_size)==3
                    C_temp = C_temp + substrates(si).pars.circulation_concentration*reshape(linspace(1,0,input.TME_size(3)),1,1,[]);
                end

                if substrates(si).method=="global"
                    substrates(si).concentration = accumarray(substrates(si).solver.regions(:),C_temp(:),[numel(substrates(si).concentration),1])./substrates(si).solver.region_volumes;
                else
                    substrates(si).concentration = C_temp;
                end

                % set C as circulation concentration at BV and decay to 0 at
                % furthest point in TME, interpolate between and average across
                % regions

            case "exponential_decay_from_bv"

                C_temp = zeros(input.TME_size);
                if license("test","image_processing_toolbox")
                    distance = bwdist(input.BV);
                else
                    TME_size = input.TME_size;
                    all_inds = reshape(1:prod(TME_size),TME_size);
                    xx = cell(length(TME_size),1);
                    [xx{:}] = ind2sub(TME_size,all_inds);

                    yy = cell(length(TME_size),1);
                    [yy{:}] = ind2sub(TME_size,find(input.BV));

                    distance = zeros(prod(TME_size),numel(yy{1}));
                    for i = 1:length(TME_size)
                        distance = distance + (xx{i}(:) - yy{i}(:)').^2;
                    end
                    distance = sqrt(min(distance,[],2));
                end

                C_temp(:) = substrates(si).pars.intitial_exp_decay_max * exp(-substrates(si).pars.initial_exp_decay_rate*distance*input.lattice_h);

                if substrates(si).method=="global"
                    substrates(si).concentration = accumarray(substrates(si).solver.regions(:),C_temp(:),[numel(substrates(si).concentration),1])./substrates(si).solver.region_volumes;
                else
                    substrates(si).concentration = C_temp;
                end

            otherwise

                error("No method for concentration initialization = %s",substrates(si).pars.concentration_initialization)

        end

    end
end



