%sim_spoil.m
%Jamie Near, 2014.
%
% USAGE:
% d_out = sim_spoil(d_in,H,angle)
% 
% DESCRIPTION:
% This function simulates the effect of a rotation about the z-axis.
% 
% INPUTS:
% d_in      = input density matrix structure.
% H         = Hamiltonian operator structure.
% angle     = Spoil angle (degrees).
%
% OUTPUTS:
% d_out     = output density matrix following z-rotation.

function d_out = sim_spoil(d_in,H,angle)

for m=1:length(H)
    %make vector of spoil angles;
    angle=angle*ones(length(H(m).shifts),1);
    
    %Make spoiler hamiltonian;
    spoil=zeros(2^H(m).nspins,2^H(m).nspins);
    for n=1:H(m).nspins
        spoil=spoil+((angle(n)*pi/180)*H(m).Iz(:,:,n));
    end
    
    %Do matrix multiplication;
    p=expm(1i*spoil);
    d_out{m} = p' * d_in{m} * p;
end


