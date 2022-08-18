function [out,save_t] = savePhenotypeSwitching(out,t,num_cells,L,counter,substrates,internalized_substrates,save_t,save_dt,Tend)

if any(isnan(internalized_substrates{1})) || any(isnan(internalized_substrates{2}))
    error('nan in internalized substrates')
end

out.t(end+1,1) = t;
out.num_cells(:,end+1) = num_cells;
out.cell_types_per_compartment = cat(ndims(L)+1,out.cell_types_per_compartment,L);
out.prolifs(:,end+1) = counter.prolifs;
out.contact_inhibitions(:,end+1) = counter.contact_inhibitions;
out.spontaneous_apoptosis(:,end+1) = counter.spontaneous_apoptosis;

for si = 1:numel(substrates)
    out.C{si} = cat(ndims(substrates(si).concentration)+1,out.C{si},substrates(si).concentration);
end
out.internalized_substrates(:,:,end+1) = internalized_substrates;

save_t = min(Tend,save_t + save_dt);