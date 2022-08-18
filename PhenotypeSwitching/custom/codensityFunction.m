% function [means,stds,means_by_type,cod_intensities] = ...
function means_by_type = ...
    codensityFunction(L,n,num_types) %,cod_edges,cod_intensity_mat,track_intensities)
% means = [];
% stds = [];
% cod_intensities = [];
num_per_type = zeros(num_types,1);
if islogical(L)
    sz = size(L);
    if num_types>1
        sz = sz(1:end-1);
    end
    nd = length(sz);
    locs = zeros(0,nd);
    types = zeros(0,1);
    for type_ind = 1:num_types

        new_ind = find(sliceof(L,nd+1,type_ind));
        X = cell(nd,1);
        [X{:}] = ind2sub(sz,new_ind);
        new_locs = cat(2,X{:});
        new_types = type_ind * ones(length(new_ind),1);

        locs = [locs;new_locs];
        types = [types;new_types];

        num_per_type(type_ind) = length(new_ind);
    end

else
    sz = size(L);
    nd = length(sz);
    locs = zeros(0,nd);
    types = zeros(0,1);
    for type_ind = 1:num_types

        new_ind = find(L==type_ind);
        X = cell(nd,1);
        [X{:}] = ind2sub(sz,new_ind);
        new_locs = cat(2,X{:});
        new_types = type_ind * ones(length(new_ind),1);

        locs = [locs;new_locs];
        types = [types;new_types];

        num_per_type(type_ind) = length(new_ind);
    end
end

d = locs-reshape(locs',1,nd,[]);
d = squeeze(sqrt(sum(d.^2,2)));
% d_sorted = sort(d);
% if n<=size(d_sorted,1)
%     cod = d_sorted(n,:);
% else
%     cod = d_sorted(end,:);
% end

% cod_intensities = zeros(length(cod_edges),num_types+1,num_types);
% for ti = num_types:-1:1
%     cod_temp = cod(types==ti);
%     means(ti) = mean(cod_temp);
%     stds(ti) = std(cod_temp);
%     if track_intensities
%         counts = histcounts(cod_temp,cod_edges,'Normalization','countdensity');
%         relfreqs = counts/max(counts);
%         cod_intensities(:,end,ti) = cod_intensity_mat*relfreqs';
%     end
% end

means_by_type = NaN(num_types);
for ti = num_types:-1:1
    if num_per_type(ti)==0
        continue;
    end
    d_temp = d(:,types==ti);

    for oi = num_types:-1:1
    
        if num_per_type(oi)==0
            continue;
        end
        d_temp_temp = d_temp(types==oi,:);
        d_temp_temp = sort(d_temp_temp);
        if n<=size(d_temp_temp,1)
            cod = d_temp_temp(n,:);
        else
            cod = d_temp_temp(end,:);
        end
        means_by_type(oi,ti) = mean(cod);

%         if track_intensities
%             counts = histcounts(cod,cod_edges,'Normalization','countdensity');
%             relfreqs = counts/max(counts);
%             cod_intensities(:,oi,ti) = cod_intensity_mat*relfreqs';
%         end
    end
end