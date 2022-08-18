function local_concentration = getLocalConcentration(S,cell_Lind)

if S.method=="global"
    local_concentration = S.concentration(S.solver.regions(cell_Lind));
else
    local_concentration = S.concentration(cell_Lind);
end