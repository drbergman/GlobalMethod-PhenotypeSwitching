function [NumCells,InternSubstrates,AllConcentrations] = normalizeOutputPhenotypeSwitching(out,input,nsubcohorts,nsamps)

%% normalize output
TME_size = input(1).TME_size;
Tend = input(1).Tend;
save_dt = input(1).save_dt;

N = prod(TME_size);
nt = round(Tend/(save_dt)+1);
t = linspace(0,Tend,nt);
NumCells = zeros(nt,nsubcohorts,nsamps,input(1).ntypes);
InternSubstrates = zeros(nt,nsubcohorts,nsamps,2);
AvgConc = zeros(nt,nsubcohorts,nsamps);
ConcPerivasc = zeros(nt,nsubcohorts,nsamps);
for j = 1:(nsubcohorts*nsamps)
    [i,si] = ind2sub([nsubcohorts,nsamps],j);
    for type_ind = 1:input(i).ntypes
        NumCells(:,i,si,type_ind) = interp1(out(i,si).t,out(i,si).num_cells(:,type_ind),t,'linear',out(i,si).num_cells(end,type_ind));
        InternSubstrates(:,i,si,type_ind) = interp1(out(i,si).t,cellfun(@(in) mean(in),out(i,si).internalized_substrates(:,type_ind)),t,'linear',mean(out(i,si).internalized_substrates{end,type_ind}));
    end
    if input(i).substrates(1).method=="global" && ((~isfield(input(i).substrates(1),"solver")) || isempty(input(i).substrates(1).solver.region_volumes))
        input(i).substrates(1).pars.regions_type = input(i).regions_type;
        input(i).substrates(1).pars.rectangle_compartment_size = input(i).oxygen_M;
        
        input(i).substrates(1).solver.regions = setupRegions(input(i).substrates(1),input(i).TME_size,input(i).bv_location);
        input(i).substrates(1).solver.region_volumes = accumarray(input(i).substrates(1).solver.regions(:),1,[numel(sliceof(out(i,si).C{1},1,1)),1]);
        
        %     else
        %         Regions = reshape(1:N,TME_size);
        %     end
        %     Region_volumes = accumarray(Regions(:),1,[numel(sliceof(out(i,si).C,1,1)),1]);
    elseif (input(i).substrates(1).method=="local" && ((~isfield(input(i).substrates(1),"solver")) || isempty(input(i).substrates(1).solver.region_volumes))) || input(i).substrates(1).method=="hybrid"
        input(i).substrates(1).solver.regions = reshape(1:prod(input(i).TME_size),input(i).TME_size);
        input(i).substrates(1).solver.region_volumes = reshape(ones(input(i).TME_size),[],1);
    end
    total_free_substrate = reshape(out(i,si).C{1},[],length(input(i).substrates(1).solver.region_volumes))*input(i).substrates(1).solver.region_volumes;
    AvgConc(:,i,si) = interp1(out(i,si).t,total_free_substrate,t,'linear',total_free_substrate(end))/N;

    if ~isfield(input,"BV")
        BV = setupBV(input(i));
    else
        BV = input(i).BV;
    end
    regions_with_bv = unique(input(i).substrates(1).solver.regions(BV>0));
    conc_at_bv = out(i,si).C{1}(:,regions_with_bv)*input(i).substrates(1).solver.region_volumes(regions_with_bv)/sum(input(i).substrates(1).solver.region_volumes(regions_with_bv));
    ConcPerivasc(:,i,si) = interp1(out(i,si).t,...
        conc_at_bv,...
        t,'linear',conc_at_bv(end));
end
TotalSubstrate = (N*AvgConc + sum(NumCells.*InternSubstrates,4,"omitnan"));
AllConcentrations = cat(4,ConcPerivasc,AvgConc,TotalSubstrate);