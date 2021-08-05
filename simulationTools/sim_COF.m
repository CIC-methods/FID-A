%sim_COF.m
%Dana Goerzen and Jamie Near, 2021.
%
% USAGE:
%  d_out=sim_COF(H,d_in,order)
%
% DESCRIPTION:
% This function nulls signal from any undesired coherences in a spin system
% Desired coherences are determined through extended phase graph analysis
% and pulse sequence design.
%
% INPUTS:
% H         = Hamiltonian operator structure.
% d_in      = input density matrix structure
% order     = desired coherence order that you wish to keep signal from
%
% OUTPUTS:
% d_out     = output density matrix structure

function d_out=sim_COF(H,d_in,order)
%initialize mask as permitting all coherences through, then iterate through coherence matrix in Hamiltonian to
%set any values that don't correspond to the desired coherence order to 0.
for n=1:length(H) %JN - Looping through the parts of the spin system:
    mask1=H(n).coherenceOrder==order;    
    %zero any undesired coherences
    d_temp=mask1.*d_in{n};
    d_in{n}=d_temp;
end
d_out=d_in;
end
