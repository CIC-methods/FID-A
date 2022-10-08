% run_simSteamShaped_fast.m
% Jamie Near, Sunnybrook Research Institute 2022.
% 
% USAGE:
% out=run_simSteamShaped_fast(spinSys);
% 
% DESCRIPTION:
% This script simulates a STEAM experiment with fully shaped RF 
% pulses. Coherence order filtering is employed to eliminate unwanted
% coherences.  Furthermore, simulations are run at various locations in 
% space to account for the within-voxel spatial variation of the metabolite 
% signal.  Summation across spatial positions is performed. The MATLAB 
% parallel computing toolbox (parfor loop) was used to accelerate the 
% simulations.  Acceleration is currently performed in the direction of the 
% slice selective pulse along the x-direction, but this can be changed.  Up 
% to a factor of 12 acceleration can be achieved using this approach.  If 
% the parallel processing toolbox is not available, then replace
% the "parfor" loop with a "for" loop.
% 
% INPUTS:
% To run this script, there is technically only one input argument:
% spinSys           = spin system to simulate 
%
% However, the user should also edit the following parameters as 
% desired before running the function:
% RFWaveform        = name of rf pulse waveform used for both (2nd & 3rd) selective 90 degree pulses.
% Tp                = duration of rf pulses[ms]
% Bfield            = Magnetic field strength in [T]
% Npts              = number of spectral points
% sw                = spectral width [Hz]
% Bfield            = magnetic field strength [Tesla]
% lw                = linewidth of the output spectrum [Hz]
% thkX              = slice thickness of x refocusing pulse [cm]
% thkY              = slice thickness of y refocusing pulse [cm]
% fovX              = full simulation FOV in the x direction [cm]
% fovY              = full simulation FOV in the y direction [cm]
% nX                = number of spatial grid points to simulate in x-direction
% nY                = number of spatial grid points to simulate in y-direction
% taus              = vector of pulse sequence timings  [ms]
%
% OUTPUTS:
% out               = Simulation results, summed over all space.

function out=run_simSteamShaped_fast(sys)
tic
% ************INPUT PARAMETERS**********************************
rfWaveform='ex40.b4_384_14.pta'; %name of RF pulse waveform.
Tp=1.920; %duration of RF pulses[ms]
flipAngle=90;  %Flip Angle of the RF pulses [degrees]
Npts=8192; %number of spectral points
sw=6000; %spectral width [Hz]
Bfield=2.89; %magnetic field strength [Tesla]
lw=1; %linewidth of the output spectrum [Hz]
thkX=2.5; %slice thickness of x RF pulse [cm]
thkY=2.5; %slice thickness of y RF pulse [cm]
fovX=5; %size of the full simulation Field of View in the x-direction [cm]
fovY=5; %size of the full simulation Field of View in the y-direction [cm]
nX=48; %Number of grid points to simulate in the x-direction
nY=48; %Number of grid points to simulate in the y-direction
tau1=6; %TE for STEAM sequence [ms]
tau2=32; %TM for STEAM sequence [ms]
centreFreq=2.3; %Centre frequency of simulation [ppm]
% ************END OF INPUT PARAMETERS**********************************

%set up spatial grid
x=linspace(-fovX/2,fovX/2,nX); %X positions to simulate [cm]
y=linspace(-fovY/2,fovY/2,nY); %y positions to simulate [cm]

%Load RF waveform
RF=io_loadRFwaveform(rfWaveform,'exc',0);

gamma=42577000; %gyromagnetic ratio

%Load spin systems
%load spinSystems
%sys=eval(['sys' spinSys]);

%If length of RF pulse is >200 pts, resample to 100 pts to reduce
%computational workload
if size(RF.waveform,1)>200
    RF=rf_resample(RF,100);
end

Gx=(RF.tbw/(Tp/1000))/(gamma*thkX/10000); %[G/cm]
Gy=(RF.tbw/(Tp/1000))/(gamma*thkY/10000); %[G/cm]

%Initialize structures:
d_temp=cell(1,1);
d=cell(1,1);

%loop through space: If you are using the parfor loops below, and you are 
%using an older version of MATLAB (e.g.R2012), don't forget to initialize 
%the parallel processing toolbox workers using 'matlabpool open N' (for N 
%workers, 12 max).  I don't think this is necessary for newer versions of 
%MATLAB.  

%First loop through all x-positions, simulating only the first refocusing
%pulse.  
%First loop through x direction (first refoc pulse only);

%for X=1:length(x)  %Use this if you don't have the MATLAB parallel processing toolbox
parfor X=1:length(x)  %Use this if you have the MATLAB parallel processing toolbox
        disp(['Executing X-position ' num2str(X) ' of ' num2str(length(x)) '!!!']);
        d_temp{X}=sim_steam_shaped_fastRF1(Bfield,sys,tau1,tau2,RF,Tp,x(X),Gx,flipAngle,centreFreq);
end

%calculate the average density matrix (Doing this inside a separate for 
%loop because I couldn't figure out how to do this inside the parfor loop): 
for X=1:length(x)
        d{1}=sim_dAdd(d{1},d_temp{X});
end

% %Initialize structures:
out_temp=cell(length(y),1);
out=struct([]);

%Now loop through y direction (second refoc pulse only);
%for Y=1:length(y) %Use this if you don't have the MATLAB parallel processing toolbox
parfor Y=1:length(y) %Use this if you do have the MATLAB parallel processing toolbox
%            disp(['Executing Y-position ' num2str(Y) ' of ' num2str(length(y)) '!!!']);
            out_temp{Y}=sim_steam_shaped_fastRF2(d{1},Npts,sw,Bfield,lw,sys,tau1,tau2,...
                RF,Tp,y(Y),Gy,flipAngle,centreFreq);
end

%Now combine the outputs;  Again, doing this inside a separate for loop
%becuase I can't figure out how to do this inside the parfor loop:
for Y=1:length(y) 
            out=op_addScans(out,out_temp{Y});
end

%For consistent scaling across different shaped simulations, we need to :
%1.  Scale down by the total number of simulations run (since these were
%    all added together.
numSims=(nX*nY);
out=op_ampScale(out,1/numSims);

%2.  Scale by the total size of the simulated region, relative to the size
%    of the voxel.
voxRatio=(thkX*thkY)/(fovX*fovY);
out=op_ampScale(out,1/voxRatio);

toc
end








%Nested Function #1
function d = sim_steam_shaped_fastRF1(Bfield,sys,TE,TM,RF,tp,dx,Gx,flipAngle,centreFreq)
% 
% USAGE:
% d = sim_steam_shaped_fastRF1(Bfield,sys,TE,TM,RF,tp,dx,Gx,flipAngle,centreFreq)
% 
% DESCRIPTION:
% This function simulates only the first bit of the STEAM experiment, up to 
% the beginning of the third RF pulse.  The excitation is
% simulated as an instantaneous rotation, and the refocusing pulse is
% simulated as a shaped rotation.
%
% This code is designed to be used in highly-accelerated shaped simulations,
% using the method described by Yan Zhang et al. Med Phys 2017;44(8): 
% 4169-78.
% 
% This code employs coherence order filtering to select only desired
% signal.  Finally, this code simulates the spectrum at a given point in 
% space (x), given the values of the slice selection gradient (Gx).  In 
% order to fully simulate the STEAM experiment, you have to run this
% simulation many times at various points in space (x), followed by 
% sim_steam_shaped_fastRF2.m, at all points in space (y).  
% 
% INPUTS:
% Bfield    = main magnetic field strength in [T]
% sys       = spin system definition structure
% TE        = echo time in [ms].
% TM        = mixing time in [ms].
% RF        = RF pulse definition structure for RF pulses (obtain using 'io_loadRFwaveform.m')
% tp        = RF pulse duration in [ms]
% dx        = position offset in x-direction (corresponding to 2nd STEAM pulse) [cm]
% Gx        = gradient strength for 2nd STEAM pulse [G/cm]
% flipAngle = flip angle of refocusing pulses [degrees] (Optional.  Default = 90 deg)
%
% OUTPUTS:
% out       = simulated spectrum, in FID-A structure format, using STEAM 
%             sequence.
if nargin<10
    centreFreq=2.3;
    if nargin<9
        flipAngle=90;
    end
end
  
%In the steam sequence, it can be common to use an asymmetric RF pulse for the
%90-degree pulse waveform.  In this case, it is conventional for the 2nd
%and 3rd rf pulses to be time-reversed versions of eachother, with the 2nd
%pulse (RF1) being a max-phase pulse, and the 3rd pulse (RF2) being a 
%min-phase pulse (i.e. the long tails of both pulse occurring during the TM 
%period so that the TE is minimized).  Here, check if the RF pulse is 
%asymmetric and if so, make sure that the 2nd and 3rd pulses are max-phase 
%and min-phase, respectively:
if RF.rfCentre>0.5
    RF1=rf_timeReverse(RF);
    RF2=RF;
else
    RF1=RF;
    RF2=rf_timeReverse(RF);
end

%Check that the TE and TM values are not too short
if TE<(RF1.rfCentre*tp*2)
    error('ERROR:  TE cannot be less than duration of RF pulse! ABORTING!!');
end
if TM<(RF2.rfCentre*tp*2)
    error('ERROR:  TM cannot be less than duration of RF pulse! ABORTING!!');
end

%Set centre frequency
for k=1:length(sys)
    sys(k).shifts=sys(k).shifts-centreFreq;
end

%Calculate Hamiltonian matrices and starting density matrix.
[H,d]=sim_Hamiltonian(sys,Bfield);

%Calculate new delays by subtracting the pulse duration from tau1 and tau2;
delays=zeros(2);
delays(1)=TE-(RF1.rfCentre*tp*2);
delays(2)=TM-(RF2.rfCentre*tp*2);
if sum(delays<0)
    error(['ERROR! The following taus are too short: ' num2str(find(delays<0)) '.']);
end

%BEGIN PULSE SEQUENCE************
d=sim_excite(d,H,'x');                                    %EXCITE
d=sim_COF(H,d,1);                                       %Keep only +1-order coherences
d=sim_evolve(d,H,delays(1)/2000);                         %Evolve by delays(1)/2
d=sim_gradSpoil(d,H,[Gx,0,0],[dx,0,0],tp*RF1.rfCentre);     %Prewind gradient for 2nd 90 degree pulse (Not sure why, but this only works when Gx is positive.  Intiutively, Gx amplitude should be the opposite of the slice select gradient (i.e. -Gx), but this does not seem to work). 
d=sim_shapedRF(d,H,RF1,tp,flipAngle,90,dx,Gx);             %1st shaped 90 degree selective pulse
d=sim_COF(H,d,0);                                       %Keep only 0-order coherences
d=sim_evolve(d,H,(delays(2))/1000);                       %Evolve by delays(2)
%END PULSE SEQUENCE**************

%After running this many times along x, the density matrices should be
%averaged, and then the average density matrix should be passed through
%'sim_press_shaped_fastRef2' at various different y-positions. 


end









%Nested Function #2
function out = sim_steam_shaped_fastRF2(d,n,sw,Bfield,linewidth,sys,TE,TM,RF,tp,dy,Gy,flipAngle,centreFreq)
%
% USAGE:
% out = sim_steam_shaped_fastRF(d,n,sw,Bfield,linewidth,sys,TE,TM,RF,tp,dy,Gy,flipAngle,centreFreq)
% 
% DESCRIPTION:
% This function simulates only the last bit of the STEAM experiment, from the 
% the beginning of the third STEAM pulse, to the end.  The RF 
% pulse is simulated as a shaped rotation.
%
% This code is designed to be used in highly-accelerated shaped simulations,
% using the method described by Yan Zhang et al. Med Phys 2017;44(8): 
% 4169-78.
 
% This code employs coherence order filtering to select only desired
% signal.
% 
% Finally, this code simulates the spectrum at a given point in space (y),
% given the values of the slice selection gradient (Gy).  In order
% to fully simulate the STEAM experiment, you have to first run
% sim_steam_shaped_fastRF1.m at all points in space (x), followed by 
% this code, at all points in space (y).  
% 
% INPUTS:
% d         = starting density matrix (obtained using 'sim_steam_shaped_fastRF1.m')
% n         = number of points in fid/spectrum
% sw        = desired spectral width in [Hz]
% Bfield    = main magnetic field strength in [T]
% linewidth = linewidth in [Hz]
% sys       = spin system definition structure
% TE        = echo time in [ms].
% TM        = mixing time in [ms].
% RF        = RF pulse definition structure for selective RF pulses (obtain using 'io_loadRFwaveform.m')
% tp        = RF pulse duration in [ms]
% dx        = position offset in x-direction (corresponding to first refocusing pulse) [cm]
% dy        = position offset in y-direction (corresponding to second refocusing pulse) [cm]
% Gx        = gradient strength for 2nd STEAM RF pulse [G/cm]
% Gy        = gradient strength for 3rd STEAM RF pulse [G/cm]
% flipAngle = flip angle of RF pulses [degrees] (Optional.  Default = 90 deg)
%
% OUTPUTS:
% out       = simulated spectrum, in FID-A structure format, using STEAM 
%             sequence.

if nargin<14
    centreFreq=2.3;
    if nargin<13
        flipAngle=90;
    end
end

%In the steam sequence, it can be common to use an asymmetric RF pulse for the
%90-degree pulse waveform.  In this case, it is conventional for the 2nd
%and 3rd rf pulses to be time-reversed versions of eachother, with the 2nd
%pulse (RF1) being a max-phase pulse, and the 3rd pulse (RF2) being a 
%min-phase pulse (i.e. the long tails of both pulse occurring during the TM 
%period so that the TE is minimized).  Here, use the asymmetry factor of 
%the pulse to make sure that the 2nd and 3rd pulses are max-phase and 
%min-phase, respectively:
if RF.rfCentre>0.5
    RF1=rf_timeReverse(RF);
    RF2=RF;
else
    RF1=RF;
    RF2=rf_timeReverse(RF);
end

%Check that the TE and TM values are not too short
if TE<(RF1.rfCentre*tp*2)
    error('ERROR:  TE cannot be less than duration of RF pulse! ABORTING!!');
end
if TM<(RF2.rfCentre*tp*2)
    error('ERROR:  TM cannot be less than duration of RF pulse! ABORTING!!');
end

%Set centre frequency
for k=1:length(sys)
    sys(k).shifts=sys(k).shifts-centreFreq;
end

%Calculate Hamiltonian matrices and starting density matrix.
[H]=sim_Hamiltonian(sys,Bfield);

%Calculate new delays by subtracting the pulse duration from TE and TM;
delays=zeros(2);
delays(1)=TE-(RF1.rfCentre*tp*2);
delays(2)=TM-(RF2.rfCentre*tp*2);
if sum(delays<0)
    error(['ERROR! The following delays are too short: ' num2str(find(delays<0)) '.']);
end

%BEGIN PULSE SEQUENCE************
d=sim_shapedRF(d,H,RF2,tp,flipAngle,90,dy,Gy);          %2nd shaped 90 degree selective pulse
d=sim_gradSpoil(d,H,[0,Gy,0],[0,dy,0],tp*RF1.rfCentre);   %Rewind gradient for 3rd 90 degree pulse (Not sure why, but this only works when Gy is positive.  Intiutively, Gx amplitude should be the opposite of the slice select gradient (i.e. -Gy), but this does not seem to work).
d=sim_COF(H,d,-1);                                      %Keep only -1 coherences
d=sim_evolve(d,H,delays(1)/2000);                            %Evolve by delays(1)/2
[out,~]=sim_readout(d,H,n,sw,linewidth,90);      %Readout along y (90 degree phase);
%END PULSE SEQUENCE**************

%Correct the ppm scale:
out.ppm=out.ppm-(4.65-centreFreq);

%Fill in structure header fields:
out.seq='steam';
out.te=TE;
out.tm=TM;
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
