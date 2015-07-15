%sim_excite.m
%Jamie Near, 2014.
%
% USAGE:
% d_out = sim_excite(H,axis,angle)
% 
% DESCRIPTION:
% This function simulates the effect of an ideal (instantaneous) excitation
% pulse on the density matrix.  Used in simulation tools.
% 
% INPUTS:
% d_in      = input density matrix structure.
% H         = Hamiltonian operator structure.
% axis      = Axis of rotation ('x' or 'y');
% angle     = Flip angle of excitation (degrees).  Optional.  Default=90.
% If angle is a scalar, then the same flip angle is applied to all spins 
% in the spin system.  If angle is a vector, then the elements of the 
% vector specify the flip angles to apply to each spin in the system.  In 
% this case, the length of the vector must be the same as the number of 
% spins in the spin system.  

function d_out = sim_excite(H,axis,angle)


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
for n=1:H.nspins
        if H.shifts(n)>=30
            alpha=0;
        else
            alpha=angle(n)*pi/180;
        end
        if axis == 'x' || axis == 'X'
            excite=excite+alpha*H.Ix(:,:,n);
        elseif axis == 'y' || axis == 'Y'
            excite=excite+alpha*H.Iy(:,:,n);
        end
end

[U,D]=eig(excite);D=diag(D);
d1=diag(exp(-1i*D));
d2=diag(exp(1i*D));


%d_out = expm(-1i*excite) * H.Fz * expm(1i*excite);
d_out = U*d1*U' * H.Fz * U*d2*U';

