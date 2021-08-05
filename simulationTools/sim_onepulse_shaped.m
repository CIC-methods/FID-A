% sim_onepulse_shaped.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out = sim_onepulse_shaped(n,sw,Bfield,linewidth,sys,RF,tp,phCyc,dfdx,G)
% 
% DESCRIPTION:
% %This function simulates the effect of a frequency selective or slice 
% selective excitation, followed immediately by the acquisition window.  
% This is mainly an exercise to see if I can get slice selective 
% excitation working.
% 
% Note that when simulating a frequency selective pulse, it is okay to
% specify only 8 arguments (no gradient needs to be specified).  If the
% 9th argument, G, is specified and is non-zero, then a slice selective
% pulse is assumed.  
% 
% INPUTS:
% n         = number of points in fid/spectrum
% sw        = desired spectral width in [Hz]
% Bfield    = main magnetic field strength in [T]
% linewidth = linewidth in [Hz]
% sys       = spin system definition structure
% RF        = RF pulse definition structure (obtain using 'io_loadRFwaveform.m')
% tp        = RF pulse duration in [ms]
% phCyc     = Phase of excitation rf pulse in [degrees].
% dfdx      = if simulating a frequency selective pulse, this argument 
%             should be the frequency offset [Hz].  If simulating a slice
%             selective pulse, this argument should be the position offset [cm].
% G         = gradient strength for slice-selective pulse [G/cm];
%
% OUTPUTS:
% out       = simulated spectrum, in FID-A structure format, using pulse-acquire 
%             sequence.


function out = sim_onepulse_shaped(n,sw,Bfield,linewidth,sys,RF,tp,phCyc,dfdx,G)

if nargin>9 && G~=0
    simType='g';
else
    simType='f';
end

%Set water to centre
centreFreq=4.65;
for k=1:length(sys)
    sys(k).shifts=sys(k).shifts-centreFreq;
end

%Calculate Hamiltonian matrices and starting density matrix.
[H,d]=sim_Hamiltonian(sys,Bfield);


%BEGIN PULSE SEQUENCE************
if simType=='g'
    d=sim_shapedRF(d,H,RF,tp,90,90+phCyc,dfdx,G);   %slice selective excitation
elseif simType=='f'
    d=sim_shapedRF(d,H,RF,tp,90,90+phCyc,dfdx);     %frequency selective excitation
end
[out,dout]=sim_readout(d,H,n,sw,linewidth,90);      %Readout along y (90 degree phase);
%END PULSE SEQUENCE**************

%Correct the ppm scale:
out.ppm=out.ppm-(4.65-centreFreq);

%Fill in structure header fields:
out.seq='onepulse';
out.te=tp/2;
out.sim='shaped';

%Additional fields for compatibility with FID-A processing tools.
out.sz=size(out.specs);
out.date=date;
out.dims.t=1;
out.dims.coils=0;
out.dims.averages=0;
out.dims.subSpecs=0;
out.dims.extras=0;
out.averages=1;
out.rawAverages=1;
out.subspecs=1;
out.rawSubspecs=1;
out.flags.writtentostruct=1;
out.flags.gotparams=1;
out.flags.leftshifted=0;
out.flags.filtered=0;
out.flags.zeropadded=0;
out.flags.freqcorrected=0;
out.flags.phasecorrected=0;
out.flags.averaged=1;
out.flags.addedrcvrs=1;
out.flags.subtracted=1;
out.flags.writtentotext=0;
out.flags.downsampled=0;
out.flags.isFourSteps=0;


