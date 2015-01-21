%sim_excite_arbPh.m
%Jamie Near, 2014.
%
% USAGE:
% d_out = sim_excite_arbPh(H,phase,angle)
% 
% DESCRIPTION:
% This function simulates the effect of an ideal (instantaneous) excitation
% pulse on the density matrix.  Used in simulation tools.  The phase of the
% excitation pulse can be arbitrarily chosen.  To achieve an arbitrary
% phase, this code executes a rotation about z by an angle of -phi (the rf
% pulse phase), then executes a rotation about x by an angle of "angle"
% (the flip angle of the excitation pulse), and then executes a rotation
% back about z by an angle of phi.
% 
% 
% INPUTS:
% d_in      = input density matrix structure.
% H         = Hamiltonian operator structure.
% phase     = Phase of rotation in degrees (ie. 0='x', 90='y', etc);
% angle     = Flip angle of excitation (degrees).  Optional.  Default=90.
% If angle is a scalar, then the same flip angle is applied to all spins 
% in the spin system.  If angle is a vector, then the elements of the 
% vector specify the flip angles to apply to each spin in the system.  In 
% this case, the length of the vector must be the same as the number of 
% spins in the spin system.  

function d_out = sim_excite_arbPh(H,phase,angle)


if nargin<3
    %angle is 90 degrees by default.
    angle=90*ones(length(H.shifts),1);
else
    if length(angle) == 1
        angle=angle*ones(length(H.shifts),1);
    elseif length(angle) ~= H.nspins
        error('ERROR:  The length of input variable ''angle'' must be the same as the number of spins in the spin system!  ABORTING!');
    else
        %do nothing;
    end
end


excite=zeros(2^H.nspins,2^H.nspins); 
Rz=zeros(2^H.nspins,2^H.nspins);
for n=1:H.nspins
        if H.shifts(n)>=30
            alpha=0;
        else
            alpha=angle(n)*pi/180;
        end
        excite=excite+alpha*H.Ix(:,:,n);
        Rz=Rz+(phase*pi/180)*H.Iz(:,:,n);
end
d_out = expm(-1i * Rz) * ...
        expm(-1i * excite) * ...
        expm(-1i * -Rz) * ...
        H.Fz * ...
        expm(1i * -Rz) * ...
        expm(1i * excite) * ...
        expm(1i * Rz);


