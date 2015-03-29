%sim_steam.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% out = sim_steam(n,sw,Bfield,linewidth,sys,te,tm)
% 
% DESCRIPTION:
% Simulate the STEAM sequence using ideal (instantaneous) RF pulses.  To 
% remove unwanted coherences, a 4 step phase cycle is automatically
% performed, with the first and third rf pulses being cycled by 0, 90 180,
% and 270 degrees.  THIS CODE IS NOT TESTED.  RESULTS MAY NOT BE ACCURATE!!
% 
% INPUTS:
% n         = number of points in fid/spectrum
% sw        = desired spectral width in [Hz]
% Bfield    = main magnetic field strength in [T]
% linewidth = linewidth in [Hz]
% sys       = spin system definition structure
% te        = echo time in [s]
% tm        = mixing time in [s]

function out = sim_steam(n,sw,Bfield,linewidth,sys,te,tm)

%Set water to centre
centreFreq=4.65;
sys.shifts=sys.shifts-4.65;

%Calculate Hamiltonian matrices and starting density matrix.
[H,d]=sim_Hamiltonian(sys,Bfield);

ph=[0:90:270];
out=struct();
figure; hold;
for m=1:4
    %BEGIN PULSE SEQUENCE************
    d=sim_excite_arbPh(H,ph(m));                    %EXCITE
    d=sim_evolve(d,H,te/2);                         %Evolve by te/2
    d=sim_rotate(d,H,-90,'x');                      %Second 90 degree pulse about x' axis.
    d=sim_evolve(d,H,tm);                           %Evolve by TM delay
    d=sim_rotate_arbPh(d,H,90,ph(m));               %Final 90 degree pulse about x' axis.
    d=sim_evolve(d,H,te/2);                         %Evolve by te/2
    [out_temp,dout]=sim_readout(d,H,n,sw,linewidth,90); %Readout along +y' (90 degree phase);
%END PULSE SEQUENCE**************
out=op_addScans(out,out_temp);
plot(out_temp.ppm,out_temp.specs,'color',[0.5,m/4,m/4]);
end

out=op_ampScale(out,0.25); %scale down to account for four phase cycles;

%Correct the ppm scale:
out.ppm=out.ppm-(4.65-centreFreq);

%Fill in structure header fields:
out.seq='steam';
out.te=te;
out.tm=tm;
out.sim='ideal';

%Additional fields for compatibility with FID-A processing tools.
out.sz=size(out.specs);
out.date=date;
out.dims.t=1;
out.dims.coils=0;
out.dims.averages=0;
out.dims.subSpecs=0;
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
out.flags.isISIS=0;














 