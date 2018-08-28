%sim_Hamiltonian.m
%Robin Simpson and Jamie Near, 2014.
%Kimberly Chan added separate Hamiltonian for J-coupling only, for use
%during shaped rf pulses (sim_shapedRF.m).
%
% USAGE:
% [H,d] = sim_Hamiltonian(sys,Bfield);
%
% DESCRIPTION:
% Creates the nxn Hamiltonian matrix for a spin system which can then be
% used in other functions to simulate NMR experiments.
%
% INPUTS:
% sys     = spin system definition structure.
% Bfield  = magnetic field strength (Tesla).
%
% OUTPUTS:
% H       = n x n Hamiltonian matrix for spin system.
% d       = Equilibrium density matrix.

function [H,d] = sim_Hamiltonian(sys,Bfield)

omega0 = -2*pi*Bfield*42576000;

for n=1:length(sys) %JN - Looping through the parts of the spin system:
    %initialize parameters:
    nspins = length(sys(n).shifts);
    
    %The J matrix and the shifts vector have to be
    %converted into radians and rad/s, respectively:
    H(n).J=sys(n).J*2*pi;
    H(n).shifts=sys(n).shifts;
    H(n).shifts_rads = sys(n).shifts*(omega0/10^6);
    
    
    %add some other "header" information to the Hamiltonian structure;
    H(n).nspins=nspins;
    H(n).Bfield=Bfield;
    
    %Creating the set of basis states
    A=(1:2^nspins);
    B=dec2bin(A-1);
    C=zeros(2^nspins,nspins);
    for k=1:nspins
        C(:,k)=str2num(B(:,k));
    end
    
    set=-1*(C==0)+C;
    
    
    %Creating the basis tensor
    H(n).basis = zeros(2^nspins,2^nspins,2*nspins);
    for k=1:2^nspins
        for m=1:2^nspins
            H(n).basis(k,m,:) = [set(k,:) set(m,:)];
        end
    end
    
    H(n).basisA=H(n).basis(:,:,1:nspins);
    H(n).basisB=H(n).basis(:,:,nspins+1:2*nspins);
    
    H(n).basis = 0.5*H(n).basis;
    H(n).basisA= 0.5*H(n).basisA;
    H(n).basisB= 0.5*H(n).basisB;
    diagonal = eye(2^nspins);
    
    H(n).Fz = zeros(2^nspins);
    H(n).Fy = zeros(2^nspins);
    H(n).Fx = zeros(2^nspins);
    H(n).HAB = zeros(2^nspins,2^nspins);
    H(n).Ix = zeros(2^nspins,2^nspins,nspins);
    H(n).Iy = zeros(2^nspins,2^nspins,nspins);
    H(n).Iz = zeros(2^nspins,2^nspins,nspins);
    
    %Now create the Hamiltonian and H.Fx,y,z
    
    %Create H.Fz
    H(n).Fz = diagonal.*sum(H(n).basisA,3);
    d{n}=H(n).Fz;
    
    %Store the spin-system scaling factor in the H structure, and scale the 
    %density matrix for each part of the spin system:
    H(n).scaleFactor=sys(n).scaleFactor;
    d{n}=d{n}*sys(n).scaleFactor;
    
    %Create individual H.Izs
    H(n).Iz=repmat(diagonal,[1,1,nspins]).*H(n).basisA;
    
    
    %Create H.Fy and H.Fx
    %We also need H.Iy for each spin so create large array containing all
    %the H.Iys piled on top of each other
    %We also need H.Ix for each spin so create large array containing all
    %the H.Ixs piled on top of each other
    for q=1:nspins
        deltaFy = -0.5*1i*ones(2^nspins,2^nspins,nspins);
        deltaFx = 0.5*ones(2^nspins,2^nspins,nspins);
        deltaIy = -0.5i*ones(2^nspins,2^nspins,nspins);
        deltaIx = 0.5*ones(2^nspins,2^nspins,nspins);
        for p=1:nspins
            if p == q
                deltaFy(:,:,q) = deltaFy(:,:,q).*...
                    ((H(n).basis(:,:,p)==H(n).basis(:,:,p+nspins)+1) - ...
                    (H(n).basis(:,:,p)==H(n).basis(:,:,p+nspins)-1));
                
                deltaFx(:,:,q) = deltaFx(:,:,q).*...
                    ((H(n).basis(:,:,p)==H(n).basis(:,:,p+nspins)+1) +...
                    (H(n).basis(:,:,p)==H(n).basis(:,:,p+nspins)-1));
                
            else
                deltaFy(:,:,q) = deltaFy(:,:,q).*...
                    (H(n).basis(:,:,p)==H(n).basis(:,:,p+nspins));
                
                deltaFx(:,:,q) = deltaFx(:,:,q).*...
                    (H(n).basis(:,:,p)==H(n).basis(:,:,p+nspins));
                
            end
            %Also, now start constructing the Hamiltonian
            %First, the z component of JI.S - only appears on the diagonal
            dotzcomp = H(n).J(q,p)*H(n).basis(:,:,q).*H(n).basis(:,:,p);
            H(n).HAB(:,:) = H(n).HAB(:,:) + diagonal(:,:).*dotzcomp;
            H(n).HABJonly=H(n).HAB;  %Added by Kimberly Chan.
            
        end
        
        H(n).Fy(:,:) = H(n).Fy(:,:) + deltaFy(:,:,q);
        H(n).Fx(:,:) = H(n).Fx(:,:) + deltaFx(:,:,q);
        H(n).Iy(:,:,q) = H(n).Iy(:,:,q) + deltaFy(:,:,q);
        H(n).Ix(:,:,q) = H(n).Ix(:,:,q) + deltaFx(:,:,q);
        
    end
    
    
    
    %Now the resonance component - each spin has a resonance frequency
    %omega0-chemshift. This component is resfreq*H.Iz for each spin. In fact
    %we assume we are in a frame rotating at omega0 so we only have to
    %include H.shifts*H.Iz for each spin. This means we are not drastically
    %undersampled as we would otherwise be.
    for e=1:nspins
        rescomp = (H(n).shifts_rads(e))*H(n).basis(:,:,e);
        H(n).HAB(:,:) = H(n).HAB(:,:) + diagonal(:,:).*rescomp;
    end
    
    
    
    %Finally the off-diagonal components which are the x and y components
    %of JI.S
    for t=1:nspins
        for u=1:nspins
            
            deltaterma = ones(2^nspins,2^nspins);
            deltatermb = ones(2^nspins,2^nspins);
            
            for p=1:nspins
                if p==u
                    deltaterma = deltaterma.*...
                        (H(n).basis(:,:,p)==H(n).basis(:,:,p+nspins)+1);
                elseif p==t
                    deltaterma = deltaterma.*...
                        (H(n).basis(:,:,p)==H(n).basis(:,:,p+nspins)-1);
                elseif p~=u && p~=t
                    deltaterma = deltaterma.*...
                        (H(n).basis(:,:,p)==H(n).basis(:,:,p+nspins));
                end
                
            end
            
            for q=1:nspins
                if q==u
                    deltatermb = deltatermb.*...
                        (H(n).basis(:,:,q)==H(n).basis(:,:,q+nspins)-1);
                elseif q==t
                    deltatermb = deltatermb.*...
                        (H(n).basis(:,:,q)==H(n).basis(:,:,q+nspins)+1);
                elseif q~=u && q~=t
                    deltatermb = deltatermb.*...
                        (H(n).basis(:,:,q)==H(n).basis(:,:,q+nspins));
                end
                
            end
            
            H(n).HAB(:,:) = H(n).HAB(:,:) + (deltaterma+deltatermb)*0.5*H(n).J(t,u).*...
                ((H(n).basis(:,:,t)==(H(n).basis(:,:,t+nspins)+1)).*...
                (H(n).basis(:,:,u)==(H(n).basis(:,:,u+nspins)-1))+ ...
                (H(n).basis(:,:,t)==(H(n).basis(:,:,t+nspins)-1)).*...
                (H(n).basis(:,:,u)==(H(n).basis(:,:,u+nspins)+1)));
            
            %Added by Kimberly Chan:
            H(n).HABJonly=H(n).HABJonly + (deltaterma+deltatermb)*0.5*H(n).J(t,u).*...
                ((H(n).basis(:,:,t)==(H(n).basis(:,:,t+nspins)+1)).*...
                (H(n).basis(:,:,u)==(H(n).basis(:,:,u+nspins)-1))+ ...
                (H(n).basis(:,:,t)==(H(n).basis(:,:,t+nspins)-1)).*...
                (H(n).basis(:,:,u)==(H(n).basis(:,:,u+nspins)+1)));
            
            
        end
    end
end




