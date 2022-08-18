function [event,current_type,cell_ind,dt] = eventSelectionPhenotypeSwitching(event_pars,num_cells,quiescence_rates)

cell_event_rates = num_cells.*[event_pars.rp,event_pars.rm,event_pars.rd];

fixed_rate_cumsum = cumsum(cell_event_rates(:)); % cumsum of fixed rate events
fixed_rate_sum = fixed_rate_cumsum(end); 

if event_pars.quiescence_as_event && num_cells(1)>0
    quiescence_rates_sum = sum(quiescence_rates);
    total_rate_sum = fixed_rate_sum + quiescence_rates_sum;
else
    total_rate_sum = fixed_rate_sum;
end

temp = rand()*total_rate_sum;
ind = find(temp<=fixed_rate_cumsum,1);

if ~isempty(ind)

    [current_type,event] = ind2sub(size(cell_event_rates),ind);

    cell_ind = randi(num_cells(current_type));

else % then a cell is becoming quiescent

    temp = temp-fixed_rate_sum;
    cell_ind = find(temp<=cumsum(quiescence_rates),1);
    current_type = 1;
    if isempty(cell_ind)
        error('should have found a cell to become quiescent here');
    end
    event = 4;

end

dt = -log(rand())/total_rate_sum;