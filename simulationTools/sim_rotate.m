%sim_rotate.m
%Jamie Near, 2014.
%
% USAGE:
% d_out = sim_rotate(d_in,H,anglein,axis)
% 
% DESCRIPTION:
% This function simulates the effect of an ideal (instantaneous) rotation
% on the density matrix.  Used in simulation tools.
% 
% INPUTS:
% d_in      = input density matrix structure.
% H         = Hamiltonian operator structure.
% anglein   = RF pulse flip angle (degrees).  If anglein is a scalar, then 
%             the same flip angle is applied to all spins in the spin 
%             system.  If anglein is a cell array of scalars, each part of 
%             the spin system will have the same flip angle applied to each 
%             of the spins in that part.  Finally, if anglein is a cell 
%             array of vectors, then a unique flip angle can be applied to 
%             every spin in each part of the spin system.  In this case, 
%             the length of the vectors in the cell array must be equal to 
%             the number of spins in the corresponding part of the spin 
%             system. 
% axis      = Axis of rotation ('x', 'y' or 'z'); (A z-rotation technically
%             doesn't correspond to an rf pulse roatation, but it is 
%             included here anyway).  
%
% OUTPUTS:
% d_out     = output density matrix following rf rotation.

function d_out = sim_rotate(d_in,H,anglein,axis)


if isnumeric(anglein)
    if length(anglein)==1 %Anglein is a scalar
        for n=1:length(H)
            angle{n}=anglein*ones(length(H(n).shifts),1);
        end
    else
        error('ERROR:  anglein must be either a numeric scalar, or a cell array.  ABORTING!');
    end
elseif iscell(anglein)
    if all(cellfun('length',anglein)==1) & length(anglein) == length(H)
        %Anglein is a cell array of scalars whose length is the same as the
        %number of parts of the spin system.  I suspect that this option
        %will be rarely used.
        for n=1:length(H)
            angle{n}=anglein{n}*ones(length(H(n).shifts),1);
        end
    elseif length(anglein) == 1 & cellfun('length',anglein) == 1
        %anglein is a cell array with a single element that is a scalar
        for n=1:length(H)
            angle{n}=anglein{1}*ones(length(H(n).shifts),1);
        end
    else  %anglein is a cell array with arrays as elements.  
        % Each array element corresponds to the flip angle for a particular
        % spin in the spin system.
        for n=1:length(H)
            if length(anglein{n}) ~= H(n).nspins
                error('ERROR:  The length of input variable ''anglein'' must be the same as the number of spins in the spin system!  ABORTING!');
            end
        end
        angle=anglein;
    end
end

for m=1:length(H) %JN - loop through the parts of the spin system
    rotMat=zeros(2^H(m).nspins,2^H(m).nspins);
    for n=1:H(m).nspins
        if H(m).shifts(n)>=30
            theta=0;
        else
            theta=angle{m}(n)*pi/180;
        end
        if axis == 'x' || axis == 'X'
            rotMat=rotMat+theta*H(m).Ix(:,:,n);
        elseif axis == 'y' || axis == 'Y'
            rotMat=rotMat+theta*H(m).Iy(:,:,n);
        elseif axis == 'z' || axis == 'Z'
            rotMat=rotMat+theta*H(m).Iz(:,:,n);
        end
    end
    p=expm(1i*rotMat);
    d_out{m} = p'* d_in{m} * p;
end



    