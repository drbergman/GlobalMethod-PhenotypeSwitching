function dx = agentODEFunction(x,secretion_rate,uptake_rate,metabolism_rate)

dx = [secretion_rate*[-1;1] - metabolism_rate*[1;0],uptake_rate*[1;-1]]*x;