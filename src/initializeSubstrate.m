function S = initializeSubstrate(input,name)

S = emptySubstrate();
S.name = name;
input_pars = fieldnames(input);
input_pars = input_pars(startsWith(input_pars,name));

name_length = strlength(name);
for i = 1:length(input_pars)
    S_name = input_pars{i}(name_length+2:end);
    if ~isfield(S.pars,S_name)
        warning('%s is not an expected field for a substrate',S_name)
    end
    S.pars.(S_name) = input.(input_pars{i});
end

S.pars.ntypes = numel(S.pars.types);
S.pars.type_inds = zeros(S.pars.ntypes,1,"uint8");
for i = 1:S.pars.ntypes
    S.pars.type_inds(i) = input.type_dict(S.pars.types(i));
end
S.pars.diffusion = S.pars.diffusion_microns / (input.lattice_h^2);
S.pars.disc_center_lattice = S.pars.disc_center / input.lattice_h + ones(1,length(input.TME_size));
S.pars.disc_radius_lattice = S.pars.disc_radius / input.lattice_h;

if input.bv_location=="none"
    S.pars.is_pk = false;
end

if numel(S.pars.uptake_rate)==1
    S.pars.uptake_rate = S.pars.uptake_rate * ones(1,S.pars.ntypes);
end
if numel(S.pars.secretion_rate)==1
    S.pars.secretion_rate = S.pars.secretion_rate * ones(1,S.pars.ntypes);
end
if numel(S.pars.export_rate)==1
    S.pars.export_rate = S.pars.export_rate * ones(1,S.pars.ntypes);
end

%% choose method
if S.pars.regions_type=="rectangular_grid"
    if all(S.pars.rectangle_compartment_size==1)
        S.method = "local";
    elseif ~S.pars.use_regions_for_diffusion
        S.method = "hybrid";
    else
        S.method = "global";
    end
elseif S.pars.regions_type=="concentric_circles"
    if S.pars.num_concentric_circle_regions==Inf
        S.method = "local";
    elseif ~S.pars.use_regions_for_diffusion
        S.method = "hybrid";
    else
        S.method = "global";
    end
elseif S.pars.regions_type=="concentric_rectangles"
    warning("no way to make concentric_rectangles use the local method yet")
    if ~S.pars.use_regions_for_diffusion
        S.method = "hybrid";
    else
        S.method = "global";
    end
elseif S.pars.regions_type=="distance_to_dc"
    if S.pars.distance_delta==0
        S.method = "local";
    elseif ~S.pars.use_regions_for_diffusion
        S.method = "hybrid";
    else
        S.method = "global";
    end
end

switch S.method
    case "global"
        S = setupGlobalSolver(S,input.TME_size,input.bv_location,input.BV,input.substrate_time_advance);
    case "hybrid"
        S = setupHybridSolver(S,input.TME_size,input.bv_location,input.BV,input.substrate_time_advance);
    case "local"
        S = setupLocalSolver(S,input.TME_size,input.BV,input.substrate_time_advance);
    otherwise
        error("what even is this?")
end
