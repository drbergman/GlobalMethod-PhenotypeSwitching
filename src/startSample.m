function [out,input] = startSample(input)

start_timer = tic;

input = input.parameter_setup_fn(input);
input.BV = setupBV(input);

input = setupSolvers(input);

out = input.sim_function(input);

out.wall_time = toc(start_timer);