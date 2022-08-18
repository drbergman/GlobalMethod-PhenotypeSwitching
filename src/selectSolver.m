function input = selectSolver(input)

for si = 1:length(input.substrates.names)
    if input.substrates.(input.substrates.names(si)).gm
        if input.substrates.(input.substrates.names(si)).use_regions_for_diffusion
            input.substrates.(input.substrates.names(si)).method = "global";
        else
            input.substrates.(input.substrates.names(si)).method = "hybrid";
        end
    else
        input.substrates.(input.substrates.names(si)).method = "local";
    end

    input.substrates.(input.substrates.names(si)) = chooseTimeSteps(input.substrates.(input.substrates.names(si)),input.substrate_time_advance);
end