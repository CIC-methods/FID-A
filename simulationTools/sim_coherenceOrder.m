%sim_coherenceOrder.m
%Dana Goerzen and Jamie Near, 2021.
%
% USAGE:
% out = sim_coherenceOrder(spinSys);
%
% DESCRIPTION:
% Creates a n x n matrix corresponding to the coherence order of each
% proton in spin system. For use in nulling signal incoherences during
% pulse sequence simulations
%
% INPUTS:
% spinSys     = spin system definition structure.
%
% OUTPUTS:
% out       = n x n coherence order matrix for spin system.

function out=sim_coherenceOrder(spinSys)
out=cell(length(spinSys),1);
for m=1:length(spinSys)
    Nspins=length(spinSys(m).J);

    p0=[0 -1;1 0];
    p=p0;
    for n=1:Nspins-1
        p=kron(ones(2),p)+kron(p0,ones(size(p)));
    end
    out{m,1}=p;
end
out=p;
end
