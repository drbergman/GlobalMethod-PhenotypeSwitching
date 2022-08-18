function M = coefficientMatrixRegions(S,ndims)

M = sparse(2*ndims*S.pars.diffusion*(S.solver.M' - diag(ones(S.solver.nregions,1))) - diag(S.pars.degradation + S.pars.fluid_exchange_rate_all_regions));

if S.pars.contains_dcs
    M(S.pars.dc_inds,:) = 0;
end