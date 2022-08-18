function dx = globalExchangeFunction(x,secretion_rate,uptake_rate,sz)

x = reshape(x,sz);

dx = ( secretion_rate .* x(:,:,2)  -  uptake_rate .* x(:,:,1) ).*reshape([1,-1],[1,1,2]);

dx = dx(:);

%% what the above code does but with naming the slices of x
% C = x(:,:,1); % ambient concentration in each (region,type) pair
% IS = x(:,:,2); % internalized substrate concentration for each
% (region,type) pair
% 
% sec_rate = secretion_rate * IS;
% up_rate = uptake_rate * C;
