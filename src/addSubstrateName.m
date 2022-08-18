function substrates = addSubstrateName(substrates,substrate_name)

if ~isfield(substrates,"names")
    substrates.names = substrate_name;
else
    substrates.names = [substrates.names,substrate_name];
end

substrates.(substrate_name).name = substrate_name;