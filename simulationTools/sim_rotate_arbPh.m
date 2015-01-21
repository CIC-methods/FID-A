%sim_rotate_arbPh.m
%Jamie Near, 2014.
%
% USAGE:
% d_out = sim_rotate_arbPh(d_in,H,angle,ph)
% 
% DESCRIPTION:
% This function simulates the effect of an ideal (instantaneous) rotation
% on the density matrix.  Used in simulation tools.  The phase of the
% rf pulse can be arbitrarily chosen.  To achieve an arbitrary
% phase, this code executes a rotation about z by an angle of -phi (the rf
% pulse phase), then executes a rotation about x by an angle of "angle"
% (the flip angle of the pulse), and then executes a rotation
% back about z by an angle of phi.
% 
% INPUTS:
% d_in      = input density matrix structure.
% H         = Hamiltonian operator structure.
% angle     = RF pulse flip angle (degrees).  If this value is a scalar,
%               then the same flip angle will be applied to all spins in
%               the system.  To apply a different flip angle to the
%               different spins in the system, the angle variable can be a
%               vector of flip angles with length equal to the H.nspins.
% ph        = Phase of rotation (in degrees; ie.  0='x', 90='y');

function d_out = sim_rotate_arbPh(d_in,H,angle,ph)


if length(angle) == 1
    angle=angle*ones(length(H.shifts),1);
elseif length(angle) ~= H.nspins
    error('ERROR:  The length of input variable ''angle'' must be the same as the number of spins in the spin system!  ABORTING!');
else
    %do nothing;
end

rotMat=zeros(2^H.nspins,2^H.nspins); 
Rz=zeros(2^H.nspins,2^H.nspins);
for n=1:H.nspins
        if H.shifts(n)>=30
            theta=0;
        else
            theta=angle(n)*pi/180;
        end
        rotMat=rotMat+theta*H.Ix(:,:,n);
        Rz=Rz+(ph*pi/180)*H.Iz(:,:,n);
end
d_out = ...
    expm(-1i * Rz) * ...
    expm(-1i * rotMat) * ...
    expm(-1i * -Rz) * ...
    d_in * ...
    expm(1i * -Rz) * ...
    expm(1i * rotMat) * ...
    expm(1i * Rz);



    