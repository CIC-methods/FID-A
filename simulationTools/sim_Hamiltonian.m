%sim_Hamiltonian.m
%Robin Simpson and Jamie Near, 2014.
%Kimberly Chan added separate Hamiltonian for J-coupling only, for use
%during shaped rf pulses (sim_shapedRF.m).
%Dana Goerzen added spin system coherence order matrix field for simulation
%coherence selection. 
%Brenden Kadota refactored basis, HAB, HABJonly calculations
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
    % add coherence order field to Hamiltonian structure
    H(n).coherenceOrder=sim_coherenceOrder(sys(n));

    %The J matrix and the shifts vector have to be
    %converted into radians and rad/s, respectively:
    H(n).J=sys(n).J*2*pi;
    H(n).shifts=sys(n).shifts;
    H(n).shifts_rads = sys(n).shifts*(omega0/10^6);

    %add some other "header" information to the Hamiltonian structure;
    H(n).nspins=nspins;
    H(n).Bfield=Bfield;

    %Preallocate space for varibles
    H(n).Fz = zeros(2^nspins);
    H(n).Fy = zeros(2^nspins);
    H(n).Fx = zeros(2^nspins);
    H(n).HAB = zeros(2^nspins,2^nspins);
    H(n).Ix = zeros(2^nspins,2^nspins,nspins);
    H(n).Iy = zeros(2^nspins,2^nspins,nspins);
    H(n).Iz = zeros(2^nspins,2^nspins,nspins);
    H(n).HABJonly = zeros(2^nspins,2^nspins);


    %Basis states
    I0=[1 0;0 1];
    Ix=0.5*[0 1;1 0];
    Iy=(1i/2)*[0 1;-1 0];
    Iz=(1/2)*[-1 0;0 1];

    %loop through all spins
    for spin = 1:nspins

        if spin == 1
            %if spin is one, set temp to be Ix, Iy, Iz
            temp_Ix = Ix;
            temp_Iy = Iy;
            temp_Iz = Iz;
        else
            %else set to be I0
            temp_Ix = I0;
            temp_Iy = I0;
            temp_Iz = I0;
        end

        %loop through spins again starting at 2
        for pos = 2:nspins

            if spin == pos
                temp_Ix = kron(temp_Ix, Ix);
                temp_Iy = kron(temp_Iy, Iy);
                temp_Iz = kron(temp_Iz, Iz);
            else
                temp_Ix = kron(temp_Ix, I0);
                temp_Iy = kron(temp_Iy, I0);
                temp_Iz = kron(temp_Iz, I0);
            end
        end
        %Safe final temp value to Ix, Iy, Iz based on spin number
        H(n).Ix(:,:,spin) = temp_Ix;
        H(n).Iy(:,:,spin) = temp_Iy;
        H(n).Iz(:,:,spin) = temp_Iz;
    end
    %Sum along the third dimension to get F 
    H(n).Fx = sum(H(n).Ix, 3);
    H(n).Fy = sum(H(n).Iy, 3);
    H(n).Fz = sum(H(n).Iz, 3);

    d{n} = H(n).Fz;
    %Store the spin-system scaling factor in the H structure, and scale the
    %density matrix for each part of the spin system:
    H(n).scaleFactor=sys(n).scaleFactor;
    d{n}=d{n}*sys(n).scaleFactor;

    %Multiply shift rads by Iz
    for spin = 1:nspins
        H(n).HAB = H(n).HAB + H(n).shifts_rads(spin)*H(n).Iz(:,:,spin);
    end


    %For every J(i,j) multiply by IS of spins I and J.
    for i = 1:nspins
        for j = i+1:nspins
            IS = H(n).Iy(:,:,i)*H(n).Iy(:,:,j) + H(n).Iz(:,:,i)*H(n).Iz(:,:,j) + H(n).Ix(:,:,i)*H(n).Ix(:,:,j);
            JIS = sys(n).J(i,j)*2*pi*IS;
            H(n).HABJonly = H(n).HABJonly + JIS;
            H(n).HAB = H(n).HAB + JIS;
        end
    end
end
end





