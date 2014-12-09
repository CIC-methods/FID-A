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

function d_out = sim_spoil(d_in,H,angle)

%make vector of spoil angles;
angle=angle*ones(length(H.shifts),1);

%Make spoiler hamiltonian;
spoil=zeros(2^H.nspins,2^H.nspins); 
for n=1:H.nspins
    spoil=spoil+((angle(n)*pi/180)*H.Iz(:,:,n));
end

%Do matrix multiplication;
d_out = expm(-1i*spoil) * d_in * expm(1i*spoil);


