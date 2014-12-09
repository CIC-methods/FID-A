%sim_rotate.m
%Jamie Near, 2014.
%
% USAGE:
% d_out = sim_rotate(d_in,H,angle,axis)
% 
% DESCRIPTION:
% This function simulates the effect of an ideal (instantaneous) rotation
% on the density matrix.  Used in simulation tools.
% 
% INPUTS:
% d_in      = input density matrix structure.
% H         = Hamiltonian operator structure.
% angle     = RF pulse flip angle (degrees).  If this value is a scalar,
%               then the same flip angle will be applied to all spins in
%               the system.  To apply a different flip angle to the
%               different spins in the system, the angle variable can be a
%               vector of flip angles with length equal to the H.nspins.
% axis      = Axis of rotation ('x', 'y' or 'z'); (A z-rotation technically
%               doesn't correspond to an rf pulse roatation, but it is included here
%               anyway).  

function d_out = sim_rotate(d_in,H,angle,axis)


if length(angle) == 1
    angle=angle*ones(length(H.shifts),1);
elseif length(angle) ~= H.nspins
    error('ERROR:  The length of input variable ''angle'' must be the same as the number of spins in the spin system!  ABORTING!');
else
    %do nothing;
end

rotMat=zeros(2^H.nspins,2^H.nspins); 
for n=1:H.nspins
        if H.shifts(n)>=30
            theta=0;
        else
            theta=angle(n)*pi/180;
        end
        if axis == 'x' || axis == 'X'
            rotMat=rotMat+theta*H.Ix(:,:,n);
        elseif axis == 'y' || axis == 'Y'
            rotMat=rotMat+theta*H.Iy(:,:,n);
        elseif axis == 'z' || axis == 'Z'
            rotMat=rotMat+theta*H.Iz(:,:,n);
        end
end
d_out = expm(-1i*rotMat) * d_in * expm(1i*rotMat);



    