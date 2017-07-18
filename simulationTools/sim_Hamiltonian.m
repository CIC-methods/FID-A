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

%initialize parameters:
nspins = length(sys.shifts);
omega0 = -2*pi*Bfield*42576000;

%The J matrix and the shifts vector have to be
%converted into radians and rad/s, respectively:
H.J=sys.J*2*pi;
H.shifts=sys.shifts;
H.shifts_rads = sys.shifts*(omega0/10^6);


%add some other "header" information to the Hamiltonian structure;
H.nspins=nspins;
H.Bfield=Bfield;

%Creating the set of basis states
A=(1:2^nspins);
B=dec2bin(A-1);
C=zeros(2^nspins,nspins);
for n=1:nspins
    C(:,n)=str2num(B(:,n));
end

set=-1*(C==0)+C;


%Creating the basis tensor
H.basis = zeros(2^nspins,2^nspins,2*nspins);
for n=1:2^nspins
    for m=1:2^nspins
        H.basis(n,m,:) = [set(n,:) set(m,:)];
    end
end

H.basisA=H.basis(:,:,1:nspins);
H.basisB=H.basis(:,:,nspins+1:2*nspins);

H.basis = 0.5*H.basis;
H.basisA= 0.5*H.basisA;
H.basisB= 0.5*H.basisB;
diagonal = eye(2^nspins);

H.Fz = zeros(2^nspins);
H.Fy = zeros(2^nspins);
H.Fx = zeros(2^nspins);
H.HAB = zeros(2^nspins,2^nspins);
H.Ix = zeros(2^nspins,2^nspins,nspins);
H.Iy = zeros(2^nspins,2^nspins,nspins);
H.Iz = zeros(2^nspins,2^nspins,nspins);

%Now create the Hamiltonian and H.Fx,y,z

%Create H.Fz
H.Fz = diagonal.*sum(H.basisA,3);
d=H.Fz;

%Create individual H.Izs
H.Iz=repmat(diagonal,[1,1,nspins]).*H.basisA;


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
                ((H.basis(:,:,p)==H.basis(:,:,p+nspins)+1) - ...
                (H.basis(:,:,p)==H.basis(:,:,p+nspins)-1));
            
            deltaFx(:,:,q) = deltaFx(:,:,q).*...
                ((H.basis(:,:,p)==H.basis(:,:,p+nspins)+1) +...
                (H.basis(:,:,p)==H.basis(:,:,p+nspins)-1));
            
        else
            deltaFy(:,:,q) = deltaFy(:,:,q).*...
                (H.basis(:,:,p)==H.basis(:,:,p+nspins));
            
            deltaFx(:,:,q) = deltaFx(:,:,q).*...
                (H.basis(:,:,p)==H.basis(:,:,p+nspins));
            
        end
        %Also, now start constructing the Hamiltonian
        %First, the z component of JI.S - only appears on the diagonal
        dotzcomp(:,:) = H.J(q,p)*H.basis(:,:,q).*H.basis(:,:,p);
        H.HAB(:,:) = H.HAB(:,:) + diagonal(:,:).*dotzcomp;
        H.HABJonly=H.HAB;  %Added by Kimberly Chan.
        
    end
    
    H.Fy(:,:) = H.Fy(:,:) + deltaFy(:,:,q);
    H.Fx(:,:) = H.Fx(:,:) + deltaFx(:,:,q);
    H.Iy(:,:,q) = H.Iy(:,:,q) + deltaFy(:,:,q);
    H.Ix(:,:,q) = H.Ix(:,:,q) + deltaFx(:,:,q);
    
end



%Now the resonance component - each spin has a resonance frequency
%omega0-chemshift. This component is resfreq*H.Iz for each spin. In fact
%we assume we are in a frame rotating at omega0 so we only have to
%include H.shifts*H.Iz for each spin. This means we are not drastically
%undersampled as we would otherwise be.
for e=1:nspins
    rescomp(:,:) = (H.shifts_rads(e))*H.basis(:,:,e);
    H.HAB(:,:) = H.HAB(:,:) + diagonal(:,:).*rescomp;
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
                    (H.basis(:,:,p)==H.basis(:,:,p+nspins)+1);
            elseif p==t
                deltaterma = deltaterma.*...
                    (H.basis(:,:,p)==H.basis(:,:,p+nspins)-1);
            elseif p~=u && p~=t
                deltaterma = deltaterma.*...
                    (H.basis(:,:,p)==H.basis(:,:,p+nspins));
            end
            
        end
        
        for q=1:nspins
            if q==u
                deltatermb = deltatermb.*...
                    (H.basis(:,:,q)==H.basis(:,:,q+nspins)-1);
            elseif q==t
                deltatermb = deltatermb.*...
                    (H.basis(:,:,q)==H.basis(:,:,q+nspins)+1);
            elseif q~=u && q~=t
                deltatermb = deltatermb.*...
                    (H.basis(:,:,q)==H.basis(:,:,q+nspins));
            end
            
        end
        
        H.HAB(:,:) = H.HAB(:,:) + (deltaterma+deltatermb)*0.5*H.J(t,u).*...
            ((H.basis(:,:,t)==(H.basis(:,:,t+nspins)+1)).*...
            (H.basis(:,:,u)==(H.basis(:,:,u+nspins)-1))+ ...
            (H.basis(:,:,t)==(H.basis(:,:,t+nspins)-1)).*...
            (H.basis(:,:,u)==(H.basis(:,:,u+nspins)+1)));
        
        %Added by Kimberly Chan:
        H.HABJonly=H.HABJonly + (deltaterma+deltatermb)*0.5*H.J(t,u).*...
            ((H.basis(:,:,t)==(H.basis(:,:,t+nspins)+1)).*...
            (H.basis(:,:,u)==(H.basis(:,:,u+nspins)-1))+ ...
            (H.basis(:,:,t)==(H.basis(:,:,t+nspins)-1)).*...
            (H.basis(:,:,u)==(H.basis(:,:,u+nspins)+1)));
        
        
    end
end

        
        
  

