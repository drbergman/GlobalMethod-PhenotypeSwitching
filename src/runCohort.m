function [input_to_return,out,use_parpool] = runCohort(input_constant,input_vary,nsamps,min_parfor_num)

subcohort_dims = reshape(cellfun(@(n) size(input_vary.(n),1),fieldnames(input_vary)),1,[]);
nsubcohorts = prod(subcohort_dims);
names_vary = fieldnames(input_vary);

total_runs = nsubcohorts*nsamps;
use_parpool = total_runs>=min_parfor_num;
input = repmat(input_constant,nsubcohorts,1);

if use_parpool
    if exist('F','var')
        delete(F)
        clear F
    end
    ppool = gcp;
    if numel(ppool.FevalQueue.QueuedFutures)>0
        delete(ppool);
        ppool = gcp;
    end
    F(1:total_runs) = parallel.FevalFuture;
else
    F = [];
end

IDX = cell(length(subcohort_dims),1);
for i = 1:nsubcohorts
    [IDX{:}] = ind2sub(subcohort_dims',i);
    for vi = 1:length(subcohort_dims)
        input(i).(names_vary{vi}) = input_vary.(names_vary{vi})(IDX{vi},:);
    end

    for si = 1:nsamps
        if use_parpool
            idx = sub2ind([nsubcohorts,nsamps],i,si);
            F(idx) = parfeval(@startSample,2,input(i));
        else
            [out(i,si),input_temp] = startSample(input(i));
            if si==1
                input_to_return(i) = input_temp;
            end
        end
    end
end

if use_parpool
    for j = 1:total_runs
        [idx,out_temp,input_temp] = fetchNext(F);
        [i,si] = ind2sub([nsubcohorts,nsamps],idx);
        out(i,si) = out_temp;
        if si==1
            input_to_return(i) = input_temp;
        end
    end
end

if length(subcohort_dims)>1
    out = reshape(out,[subcohort_dims,nsamps]);
end