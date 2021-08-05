% run_simMegaPressShapedRefoc_fast.m
% Jamie Near, McGill University 2018.
% 
% USAGE:
% [outDIFF,outOFF,outON]=run_simMegaPressShapedRefoc_fast(spinSys);
% 
% DESCRIPTION:
% This script simulates a MEGA-PRESS experiment with fully shaped 
% refocusing pulses.  Phase cycling of the refocusing pulses is performed 
% and simulations are run at various locations in space to account for the 
% within-voxel spatial variation of the GABA signal.  Summation across 
% spatial positions is performed.  To achieve faster perfomance compared to 
% the original 'run_simPressShapedRefoc.m' function, this code uses the 
% method described by Yan Zhang et al. Med Phys 2017;44(8): 4169-78.  Some 
% additional accelration is currently performed by parallelization of the 
% x- and y- spatial for loops.  To enable the use of the MATLAB parallel 
% computing toolbox, initialize the multiple worked nodes using "matlabpool 
% size X", where "X" is the number of available processing nodes.  If the 
% parallel processing toolbox is not available, then replace the "parfor" 
% loop with a "for" loop.
% 
% INPUTS:
% To run this script, there is technically only one input argument:
% spinSys           = spin system to simulate 
%
% However, the user should also edit the following parameters as desired
% before running the function:
% refocWaveform         = name of refocusing pulse waveform.
% refTp                 = duration of refocusing pulses[ms]
% Npts                  = number of spectral points
% sw                    = spectral width [Hz]
% Bfield                = magnetic field strength [Tesla]
% lw                    = linewidth of the output spectrum [Hz]
% thkX                  = slice thickness of x refocusing pulse [cm]
% thkY                  = slice thickness of y refocusing pulse [cm]
% fovX                  = full simulation FOV in the x direction [cm]
% fovY                  = full simulation FOV in the y direction [cm]
% nX                    = number of spatial grid points to simulate in x-direction
% nY                    = number of spatial grid points to simulate in y-direction
% taus                  = vector of pulse sequence timings  [ms]
% spinSys               = spin system to simulate
% editFlipON            = cell array containing the vectors of edit-ON pulse flip angles for each spin in each part of the spin system.
% editFlipOFF           = cell array containing the vectors of edit-OFF pulse flip angles for each spin in each part of the spin system.
% refPhCyc1             = vector of phase cycling steps for 1st refocusing pulse [degrees]
% refPhCyc2             = vector of phase cycling steps for 2nd refocusing pulse [degrees]
%
% OUTPUTS:
% outON             = Simulated MEGA-PRESS edit-ON spectrum, summed over
%                     all positions.
% outOFF            = Simulated MEGA-PRESS edit-OFF spectrum, summed over
%                     all positions.
% outDIFF           = Simulated MEGA-PRESS difference spectrum, summed over
%                     all positions.

function [outDIFF,outOFF,outON]=run_simMegaPressShapedRefoc_fast(spinSys);

% ************INPUT PARAMETERS**********************************
refocWaveform='sampleRefocPulse.pta'; %name of refocusing pulse waveform.
refTp=5; %duration of refocusing pulses[ms]
Npts=2048; %number of spectral points
sw=2000; %spectral width [Hz]
Bfield=3; %magnetic field strength [Tesla]
lw=2; %linewidth of the output spectrum [Hz]
thkX=3; %slice thickness of x refocusing pulse [cm]
thkY=3; %slice thickness of y refocusing pulse [cm]
fovX=5; %size of the full simulation Field of View in the x-direction [cm]
fovY=5; %size of the full simulation Field of View in the y-direction [cm]
nX=8; %Number of grid points to simulate in the x-direction
nY=8; %Number of grid points to simulate in the y-direction
taus=[5,... %Time from excitation to 1st refoc pulse [ms]
    17,...  %Time from 1st refoc pulse to 1st editing pulse [ms]
    17,...  %Time from 1st editing pulse to 2nd refoc pulse [ms]
    17,...  %Time from 2nd refoc pulse to 2nd editing pulse [ms]
    12];    %Time from 2nd editing pulse to ADC onset [ms]
editFlipON{1}=[0 0 180 180 0 0];  %Each element of this cell refers to a different part of the spin system. 
editFlipOFF{1}=[0 0 0 0 0 0];  %Each element of this cell refers to a different part of the spin system. 
refPhCyc1=[0,90]; %phase cycling steps for 1st refocusing pulse [degrees]
refPhCyc2=[0,90]; %phase cycling steps for 2nd refocusing pulse [degrees]
% ************END OF INPUT PARAMETERS**********************************

%set up spatial grid
x=linspace(-fovX/2,fovX/2,nX); %X positions to simulate [cm]
y=linspace(-fovY/2,fovY/2,nY); %y positions to simulate [cm]

%Load RF waveforms
refRF=io_loadRFwaveform(refocWaveform,'ref',0);

gamma=42577000; %gyromagnetic ratio

%Load spin systems
load spinSystems
sys=eval(['sys' spinSys]);

%Resample refocusing RF pulse from 400 pts to 100 pts to reduce
%computational workload
refRF=rf_resample(refRF,100);

Gx=(refRF.tbw/(refTp/1000))/(gamma*thkX/10000); %[G/cm]
Gy=(refRF.tbw/(refTp/1000))/(gamma*thkY/10000); %[G/cm]

%Initialize structures:
d_ON_temp=cell(length(x),length(refPhCyc1));
d_OFF_temp=cell(length(x),length(refPhCyc1));
d_ON=cell(length(refPhCyc1));
d_OFF=cell(length(refPhCyc1));

%loop through space: If you are using the parfor loops below, and you are
%using an older version of MATLAB (e.g.R2012), don't forget to initialize
%the parallel processing toolbox workers using 'matlabpool open N' (for N
%workers, 12 max).  I don't think this is necessary for newer version of
%MATLAB.

%First loop through all x-positions, simulating only the first refocusing
%pulse.
%for X=1:length(x)  %Use this if you don't have the MATLAB parallel processing toolbox
parfor X=1:length(x)  %Use this if you do have the MATLAB parallel processing toolbox
        for RP1=1:length(refPhCyc1)
                disp(['Executing X-position ' num2str(X) ' of ' num2str(length(x)) ', '...
                    'First Refoc phase cycle ' num2str(RP1) ' of ' num2str(length(refPhCyc1)) '!!!']);
                d_ON_temp{X}{RP1}=sim_megapress_shapedRefoc_fastRef1(Bfield,taus,sys,...
                    refRF,refTp,editFlipON,Gx,x(X),refPhCyc1(RP1));
                d_OFF_temp{X}{RP1}=sim_megapress_shapedRefoc_fastRef1(Bfield,taus,sys,...
                    refRF,refTp,editFlipOFF,Gx,x(X),refPhCyc1(RP1));
    end %end of 1st refocusing phase cycle loop.
end %end of spatial loop (parfor) in x direction.

%calculate the average density matrix (Doing this inside a separate for 
%loop because I couldn't figure out how to do this inside the parfor loop): 
for X=1:length(x)
    for RP1=1:length(refPhCyc1)
        d_ON{RP1}=sim_dAdd(d_ON{RP1},d_ON_temp{X}{RP1});
        d_OFF{RP1}=sim_dAdd(d_OFF{RP1},d_OFF_temp{X}{RP1});
    end
    
end

% %Initialize structures:
outON_temp=cell(length(y),length(refPhCyc1),length(refPhCyc2));
outOFF_temp=cell(length(y),length(refPhCyc1),length(refPhCyc2));
outON=struct([]);
outOFF=struct([]);

%Now loop through y direction (second refoc pulse only);
%for Y=1:length(y)  %Use this if you don't have the MATLAB parallel processing toolbox
parfor Y=1:length(y)  %Use this if you do have the MATLAB parallel processing toolbox
    for RP1=1:length(refPhCyc1)
        for RP2=1:length(refPhCyc2)
            disp(['Executing Y-position ' num2str(Y) ' of ' num2str(length(y)) ', '...
                'First Refoc phase cycle ' num2str(RP1) ' of ' num2str(length(refPhCyc1)) ', '...
                'Second Refoc phase cycle ' num2str(RP2) ' of ' num2str(length(refPhCyc2)) '!!!']);
            outON_temp{Y}{RP1}{RP2}=sim_megapress_shapedRefoc_fastRef2(d_ON{RP1},Npts,sw,Bfield,lw,taus,sys,...
                refRF,refTp,editFlipON,Gy,y(Y),refPhCyc2(RP2));
            outOFF_temp{Y}{RP1}{RP2}=sim_megapress_shapedRefoc_fastRef2(d_OFF{RP1},Npts,sw,Bfield,lw,taus,sys,...
                refRF,refTp,editFlipOFF,Gy,y(Y),refPhCyc2(RP2));
        end %end of 2nd refocusing phase cycle loop.
    end %end of 1st refocusing phase cycle loop.
end %end of spatial loop (parfor) in y direction.

%Now combine the outputs;  Again, doing this inside a separate for loop
%becuase I can't figure out how to do this inside the parfor loop:
for Y=1:length(y)
    for RP1=1:length(refPhCyc1)
        for RP2=1:length(refPhCyc2)
            outON=op_addScans(outON,outON_temp{Y}{RP1}{RP2},xor(RP1-1,RP2-1));
            outOFF=op_addScans(outOFF,outOFF_temp{Y}{RP1}{RP2},xor(RP1-1,RP2-1));
        end
    end
end

%For consistent scaling across different shaped simulations, we need to :
%1.  Scale down by the total number of simulations run (since these were
%    all added together).
numSims=(nX*nY*length(refPhCyc1)*length(refPhCyc2));
outON=op_ampScale(outON,1/numSims);
outOFF=op_ampScale(outOFF,1/numSims);

 %2.  Scale by the total size of the simulated region, relative to the size
%    of the voxel.
if fovX>thkX
    voxRatio=(thkX*thkY)/(fovX*fovY);
else
    voxRatio=1;
end
outON=op_ampScale(outON,1/voxRatio);
outOFF=op_ampScale(outOFF,1/voxRatio); 

outDIFF=op_subtractScans(outON,outOFF);

end %END OF MAIN FUNCTION:  Nested functions below.



%Nested Function #1
function d = sim_megapress_shapedRefoc_fastRef1(Bfield,taus,sys,refPulse,refTp,editFlip,Gx,dx,refPh1,flipAngle)
% 
% USAGE:
% out = sim_megapress_shapedRefoc_fastRef1(Bfield,taus,sys,refPulse,refTp,editFlip,Gx,dx,refPh1,flipAngle)
% 
% DESCRIPTION:
% This function simulates only the first bit of the MEGA-PRESS experiment, 
% up to the beginning of the second refocusing pulse.  The excitation is
% simulated as an instantaneous rotation, and the refocusing pulses are 
% simulated as shaped rotations.
%
% This code is designed to be used in highly-accelerated shaped simulations,
% using the method described by Yan Zhang et al. Med Phys 2017;44(8): 
% 4169-78.
%
% This code enables the choice of the phase of the refocusing pulses.  
% This enables phase cycling of the refocusing pulses by repeating 
% simulations with different refocusing pulse phases.  A four step 
% phase cycling scheme is typically sufficient, where both refocusing 
% pulses are phase cycled by 0 and 90 degrees, and the phase are combined 
% in the following way:
% 
% signal = ([0 90] - [0 0]) + ([90 0] - [90 90]);
% 
% where, in [X Y], X is the phase of the first refocusing pulse and Y is
% the phase of the second refocusing pulse
% 
% Finally, this code simulates the spectrum at a given point in space (x),
% given the values of the slice selection gradient (Gx).  In order
% to fully simulate the MEGA-PRESS experiment, you have to run this
% simulation many times at various points in space (x), followed by 
% sim_press_shaped_fastRef2.m, at all points in space (y).  
% 
% INPUTS:
% Bfield    = main magnetic field strength in [T]
% taus(1)     = time in [ms] from 90 to 1st 180
% taus(2)     = time in [ms] from 1st 180 to 1st edit pulse
% taus(3)     = time in [ms] from 1st edit pulse to 2nd 180
% taus(4)     = time in [ms] from 2nd 180 to 2nd edit pulse
% taus(5)     = time in [ms] from 2nd edit pulse to ADC
%               FOR MEGA-PRESS on SIEMENS SYSTEM:  taus=[4.545,12.7025,21.7975,12.7025,17.2526];
% sys        = Metabolite spin system definition structure;
% refPulse   = RF pulse definition structure for refoc pulses (obtain using 'io_loadRFwaveform.m')
% refTp      = duration of refocusing pulse in [ms]
% editFlip   = cell array containing the vectors of editing flip angles in [degrees] at chemical shifts corresponding to 'shifts' for each part of the spin system.
% Gx         = gradient strength for first selective refocusing pulse [G/cm]
% dx         = position offset in x-direction (corresponding to first refocusing pulse) [cm]
% refPh1     = the phase of the first refocusing pulse in [degrees];
% flipAngle  = the flip angle of the refocusing pulse.  
%
% OUTPUTS:
% d          = output density matrix

if nargin<10
    flipAngle=180;
end
    
%Set 3ppm GABA resonance to centre
centreFreq=3;
for k=1:length(sys)
    sys(k).shifts=sys(k).shifts-centreFreq;
end

%Calculate Hamiltonian matrices and starting density matrix.
[H,d]=sim_Hamiltonian(sys,Bfield);

%Calculate new delays by subtracting the pulse durations from the taus
%vector;
delays=zeros(size(taus));
delays(1)=taus(1)-(refTp/2);
delays(2)=taus(2)-(refTp/2);
delays(3)=taus(3)-(refTp/2);
delays(4)=taus(4)-(refTp/2);
delays(5)=taus(5);
if sum(delays<0)
    error(['ERROR! The following taus are too short: ' num2str(find(delays<0)) '.']);
end


%BEGIN PULSE SEQUENCE************
d=sim_excite(d,H,'x');                                          %EXCITE
d=sim_evolve(d,H,delays(1)/1000);                               %Evolve by delays(1)
d=sim_shapedRF(d,H,refPulse,refTp,flipAngle,90+refPh1,dx,Gx);   %1st shaped 180 degree refocusing pulse
d=sim_evolve(d,H,delays(2)/1000);                               %Evolve by delays(2)
d=sim_rotate(d,H,editFlip,'y');                                 %1st editing pulse rotation
d=sim_evolve(d,H,delays(3)/1000);                               %Evolve by delays(3)
%END PULSE SEQUENCE**************


%After running this many times along x, the density matrices should be
%averaged, and then the average density matrix should be passed through
%'sim_megapress_shaped_fastRef2' at various different y-positions. 


end
        





%Nested Function #2
function out = sim_megapress_shapedRefoc_fastRef2(d,n,sw,Bfield,linewidth,taus,sys,refPulse,refTp,editFlip,Gy,dy,refPh2,flipAngle)
%
% USAGE:
% out = sim_press_shapedRefoc_fastRef2(d,n,sw,Bfield,linewidth,taus,sys,refPulse,refTp,editFlip,Gy,dy,refPh2,flipAngle)
% 
% DESCRIPTION:
% This function takes a starting density matrix, and simulates only the 
% last bit of the MEGA-PRESS experiment, from the second refocusing pulse 
% to the end.  The refocusing pulses are simulated as shaped rotations.
%
% This code is designed to be used in highly-accelerated shaped simulations,
% using the method described by Yan Zhang et al. Med Phys 2017;44(8): 
% 4169-78.
%
% This code enables the choice of the phase of the refocusing pulses.  
% This enables phase cycling of the refocusing pulses by repeating 
% simulations with different refocusing pulse phaseS.  A four step 
% phase cycling scheme is typically sufficient, where both refocusing 
% pulses are phase cycled by 0 and 90 degrees, and the phase are combined 
% in the following way:
% 
% signal = ([0 90] - [0 0]) + ([90 0] - [90 90]);
% 
% where, in [X Y], X is the phase of the first refocusing pulse and Y is
% the phase of the second refocusing pulse
% 
% Before running this code, you need to run sim_megapress_shapedRefoc_fast1.m 
% to generate the density matrix resulting from the first part of the 
% megapress sequence over a range of x positions.  Then you need to take 
% the average density matrix across x.  Finally, this code simulates the 
% spectrum at a given point in space (y), given the values of the slice 
% selection gradient (Gy).  In order to fully simulate the MEGA-PRESS 
% experiment, you have to run this simulation many times at various points 
% in space (y).  
% 
% INPUTS:
% Bfield    = main magnetic field strength in [T]
% taus(1)     = time in [ms] from 90 to 1st 180
% taus(2)     = time in [ms] from 1st 180 to 1st edit pulse
% taus(3)     = time in [ms] from 1st edit pulse to 2nd 180
% taus(4)     = time in [ms] from 2nd 180 to 2nd edit pulse
% taus(5)     = time in [ms] from 2nd edit pulse to ADC
%               FOR MEGA-PRESS on SIEMENS SYSTEM:  taus=[4.545,12.7025,21.7975,12.7025,17.2526];
% sys        = Metabolite spin system definition structure;
% refPulse   = RF pulse definition structure for refoc pulses (obtain using 'io_loadRFwaveform.m')
% refTp      = duration of refocusing pulse in [ms]
% editFlip   = cell array containing the vectors of editing flip angles in [degrees] at chemical shifts corresponding to 'shifts' for each part of the spin system.
% Gx         = gradient strength for first selective refocusing pulse [G/cm]
% dx         = position offset in x-direction (corresponding to first refocusing pulse) [cm]
% refPh1     = the phase of the first refocusing pulse in [degrees];
% flipAngle  = the flip angle of the refocusing pulse.  
%
% OUTPUTS:
% d          = output density matrix

if nargin<14
    flipAngle=180;
end
    
%Set 3ppm GABA resonance to centre
centreFreq=3;
for k=1:length(sys)
    sys(k).shifts=sys(k).shifts-centreFreq;
end

%Calculate Hamiltonian matrices and starting density matrix.
[H]=sim_Hamiltonian(sys,Bfield);

%Calculate new delays by subtracting the pulse durations from the taus
%vector;
delays=zeros(size(taus));
delays(1)=taus(1)-(refTp/2);
delays(2)=taus(2)-(refTp/2);
delays(3)=taus(3)-(refTp/2);
delays(4)=taus(4)-(refTp/2);
delays(5)=taus(5);
if sum(delays<0)
    error(['ERROR! The following taus are too short: ' num2str(find(delays<0)) '.']);
end


%BEGIN PULSE SEQUENCE************
d=sim_shapedRF(d,H,refPulse,refTp,flipAngle,90+refPh2,dy,Gy);  %2nd shaped 180 degree refocusing pulse
d=sim_evolve(d,H,delays(4)/1000);                          %Evolve by delays(4)
d=sim_rotate(d,H,editFlip,'y');                        %2nd editing pulse rotation
d=sim_evolve(d,H,delays(5)/1000);                          %Evolve by delays(5)
[out,~]=sim_readout(d,H,n,sw,linewidth,90);           %Readout along y (90 degree phase);
%END PULSE SEQUENCE**************

%Correct the ppm scale:
out.ppm=out.ppm-(4.65-centreFreq);

%Fill in structure header fields:
out.seq='megapress';
out.te=sum(taus);
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



end


