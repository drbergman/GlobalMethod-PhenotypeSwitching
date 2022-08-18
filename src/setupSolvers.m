function input = setupSolvers(input)

for si = 1:length(input.substrate_names)
    input.substrates(si) = initializeSubstrate(input,input.substrate_names(si));
end