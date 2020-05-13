%sim_evolve.m
%Jamie Near, 2014.
%
% USAGE:
% d_out = sim_evolve(d_in,H,t)
% 
% DESCRIPTION:
% This function simulates free evolution of the spin system under the 
% effects of chemical shift and scalar coupling.
% 
% INPUTS:
% d_in      = input density matrix structure.
% H         = Hamiltonian operator structure.
% t         = duration of evolution (s)
%
% OUTPUTS:
% d_out     = output density matrix following free evolution.

function d_out = sim_evolve(d_in,H,t)

for n=1:length(H) %JN - loop through the different parts of the spin-system
    p=expm(1i*H(n).HAB*t);
    d_out{n} = p' * d_in{n} * p;  %Faster than doing expm twice.  Suggested by Martin Froeling.
end

end


