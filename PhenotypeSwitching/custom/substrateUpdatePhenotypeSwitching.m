function [substrates,internalized_substrates,L,cells_Lind,num_cells,cell_ids] = ...
    substrateUpdatePhenotypeSwitching(substrates,internalized_substrates,tf,TME_size,quiescence_as_event,L,...
    cells_Lind,BV,...
    quiescence_threshold,cell_ids)

for si = 1:numel(substrates)

    dt_substrate = tf-substrates(si).t;

    if dt_substrate == 0
        continue;
    end
    
    switch substrates(si).method
        case "global"
            [substrates(si),internalized_substrates(substrates(si).pars.type_inds,si)] = ...
                solve_global(substrates(si),internalized_substrates(substrates(si).pars.type_inds,si),...
                dt_substrate,cells_Lind(substrates(si).pars.type_inds));
        case "local"
            [substrates(si),internalized_substrates(substrates(si).pars.type_inds,si)] = ...
                solve_local(substrates(si),internalized_substrates(substrates(si).pars.type_inds,si),...
                dt_substrate,cells_Lind(substrates(si).pars.type_inds));
        case "hybrid"
            [substrates(si),internalized_substrates(substrates(si).pars.type_inds,si)] = ...
                solve_hybrid(substrates(si),internalized_substrates(substrates(si).pars.type_inds,si),...
                dt_substrate,cells_Lind(substrates(si).pars.type_inds));
        otherwise
            error("not sure what method is being used")

    end

end

% dt_oxygen = tf-substrates(1).t;
% switch substrates(1).method
%     case "global"
%         [substrates(1),internalized_substrates(substrates(1).pars.type_inds,1)] = ...
%             solve_global(substrates(1),internalized_substrates(substrates(1).pars.type_inds,1),...
%             dt_oxygen,cells_Lind,length(TME_size));
%     case "local"
%         [substrates(1),internalized_substrates(substrates(1).pars.type_inds,1)] = ...
%             solve_local(substrates(1),internalized_substrates(substrates(1).pars.type_inds,1),...
%             dt_oxygen,cells_Lind);
%     case "hybrid"
%         [substrates(1),internalized_substrates(substrates(1).pars.type_inds,1)] = ...
%             solve_hybrid(substrates(1),internalized_substrates(substrates(1).pars.type_inds,1),...
%             dt_oxygen,cells_Lind,BV,TME_size);
%     otherwise 
%         error("not sure what method is being used")
% 
% end

% if length(substrates)>=2 % then using quorum factor
%     q_internalized_substrates = cell(2,1); % healthy and quiescent cells export quorum factor and these are always in the model
%     q_internalized_substrates{1} = zeros(size(cells_Lind{1}));
%     q_internalized_substrates{2} = zeros(size(cells_Lind{2}));
%     dt_quorum = tf-substrates(2).t;
%     switch substrates(2).method
%         case "global"
%             [substrates(2),~] = solve_global(substrates(2),...
%                 q_internalized_substrates,dt_quorum,...
%                 cells_Lind(1:2),ndims(TME_size));
%         case "local"
%             [substrates(2),~] = solve_local(substrates(2),...
%                 q_internalized_substrates,dt_quorum,...
%                 cells_Lind(1:2));
%         case "hybrid"
%             [substrates(2),~] = solve_hybrid(substrates(2),...
%                 q_internalized_substrates,dt_quorum,...
%                 cells_Lind(1:2),BV,TME_size);
% 
%             otherwise
%                 error("not sure what method is being used")
% 
%     end
% end

if ~quiescence_as_event
    newly_quiescent_cells_ind = find(internalized_substrates{1}<quiescence_threshold);
else
    newly_quiescent_cells_ind = [];
end
newly_unquiescent_cells_ind = find(internalized_substrates{2}>quiescence_threshold);

if ~isempty(newly_quiescent_cells_ind) || ~isempty(newly_unquiescent_cells_ind)
    L(cells_Lind{1}(newly_quiescent_cells_ind)) = 2;
    L(cells_Lind{2}(newly_unquiescent_cells_ind)) = 1;

    cells_Lind{1} = [cells_Lind{1};cells_Lind{2}(newly_unquiescent_cells_ind)];
    cells_Lind{2} = [cells_Lind{2};cells_Lind{1}(newly_quiescent_cells_ind)];
    cells_Lind{1}(newly_quiescent_cells_ind,:) = [];
    cells_Lind{2}(newly_unquiescent_cells_ind,:) = [];

    cell_ids{1} = [cell_ids{1};cell_ids{2}(newly_unquiescent_cells_ind)];
    cell_ids{2} = [cell_ids{2};cell_ids{1}(newly_quiescent_cells_ind)];
    cell_ids{1}(newly_quiescent_cells_ind,:) = [];
    cell_ids{2}(newly_unquiescent_cells_ind,:) = [];

    internalized_substrates{1} = [internalized_substrates{1};internalized_substrates{2}(newly_unquiescent_cells_ind)];
    internalized_substrates{2} = [internalized_substrates{2};internalized_substrates{1}(newly_quiescent_cells_ind)];
    internalized_substrates{1}(newly_quiescent_cells_ind,:) = [];
    internalized_substrates{2}(newly_unquiescent_cells_ind,:) = [];
end

num_cells = [length(cells_Lind{1});length(cells_Lind{2})];
