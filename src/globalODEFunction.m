function dx = globalODEFunction(x,p,...
    M_diffusion_and_efflux,agent_density,sz)

dx = zeros(sz);

x = reshape(x,sz);

C = x(:,1);
IS = x(:,2:end);

sec_rate = p.secretion_rate .* IS;
up_rate = p.uptake_rate .* C;

dx(:,1) = M_diffusion_and_efflux * C + ...
    sum((sec_rate-up_rate + p.export_rate).*agent_density,2) + ...
    p.fluid_exchange_rate_all_regions .* p.circulation_concentration;

dx(:,2:end) = up_rate-sec_rate;

if p.use_internal_agent_ode
    dx(:,2:end) = dx(:,2:end) + reshape(p.internal_agent_ode(x(:,2:end),1:p.ntypes),[],p.ntypes);
end

if p.contains_dcs
    dx(p.dc_inds,1) = 0; % don't let DC nodes change value
end

dx = dx(:);