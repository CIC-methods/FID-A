% 2021 EDIT: This function is now deprecated. The current version of the shaped semiLASER 
% simulation is run_simSemiLASERShaped_fast.m. The current version employs coherence selection, 
% resulting in a 4x simulation speed increase.
% run_simSemiLASERShaped_fast_phCyc.m
%
% Fast version by Muhammad G Saleh (Johns Hopkins University School of Medicine, 2019)
%
% Dana Goerzen (McGill University, 2019).
%
% USAGE:
% This script runs a spatially resolved, shaped simulation of the sLASER
% experiment, using the sim_SemiLASER_shaped.m function.
% Adjust all parameters as desired and click "Run"
%
% DESCRIPTION:
% Simulates CMRR semi-LASER experiments designed by Oz et al. (2018) using
% shaped adiabatic refocusing pulses along XX YY gradients
% The excitation is simulated as an instantaneous rotation,
% and the adiabatic refocusing pulses are simulated as a shaped rotation.
% This simulation also employs a 4-step phase cycling scheme as follows:

% signal = ([0 0, 0 0] - [0 0, 0 90]) - ([0 90, 0 0] + [0 90, 0 90]);

%
% Finally, this code simulates the spectrum at a given point in space (x,y),
% given the values of the slice selection gradients (Gx, and Gy).  The pulse
% waveform is assumed to be the same for both refocusing pulses.  In order
% to fully simulate the sLASER experiment, you have to run this
% simulation many times at various points in space (x,y), and then add
% together the resulting spectra and scale down by the number of simulations.
%
% Feb 2020 - Jamie Near:  This code now accepts gradient modulated pulses.  
%

%tic;

% INPUTS:
n=4096; %= number of points in fid/spectrum
sw=5000; %= desired spectral width in [Hz]
Bfield=7; %= main magnetic field strength in [T]
lw=2; %= linewidth in [Hz]
load spinSystems.mat; %= spin system definition structure
sys=sysLac;
rfPulse=io_loadRFwaveform('sampleAFPpulse_HS2_R15.RF','inv'); % adiabatic RF pulse shaped waveform
refTp=3.5; %= RF pulse duration in [ms]
flipAngle=180; %= flip angle of refocusing pulses [degrees] (Optional.  Default = 180 deg)
centreFreq=2.3; %= centre frequency of the spectrum in [ppm] (Optional.  Default = 2.3)
thkX=2; %slice thickness of x refocusing pulse [cm]
thkY=2; %slice thickness of y refocusing pulse [cm]
fovX=3; %size of the full simulation Field of View in the x-direction [cm]
fovY=3; %size of the full simulation Field of View in the y-direction [cm]
nX=32; %Number of grid points to simulate in the x-direction
nY=32; %Number of grid points to simulate in the y-direction
te=135;         %sLASER total echo time [ms]
ph1=[0 0 0 0];  %phase cycling scheme of first refocusing pulse
ph2=[0 0 90 90]; %phase cycling scheme of second refocusing pulse
ph3=[0 0 0 0]; %phase cycling scheme of third refocusing pulse
ph4=[0 90 0 90]; %phase cycling scheme of fourth refocusing pulse

% OUTPUTS:
% out       = simulated spectrum, in FID-A structure format, using PRESS
%             sequence.

%set up spatial grid
x=linspace(-fovX/2,fovX/2,nX); %X positions to simulate [cm]
y=linspace(-fovY/2,fovY/2,nY); %y positions to simulate [cm]

gamma=42577000; %gyromagnetic ratio

%Resample refocusing RF pulse from 400 pts to 100 pts to reduce
%computational workload
rfPulse=rf_resample(rfPulse,100);

%sys=sysRef0ppm
if ~rfPulse.isGM
    %Non-gradient modulated pulse - Calculating the x and y gradient 
    %strengths for the desired slice thickness
    Gx=(rfPulse.tbw/(refTp/1000))/(gamma*thkX/10000); %[G/cm]
    Gy=(rfPulse.tbw/(refTp/1000))/(gamma*thkY/10000); %[G/cm]
else
    %Gradient modulated pulse
    %1.  Calculating the unitless scaling factor for the GM waveform.
    Gx=(rfPulse.tthk/(refTp/1000))/thkX;
    Gy=(rfPulse.tthk/(refTp/1000))/thkY;
end

%Initialize structures:
% out_posxy_rpc=cell(length(x),length(y),length(ph1));
out_posx_rpc =cell(length(x),length(ph1));
d=cell(length((ph1))); %Initialize a cell for the dentity matrix, with elements for each phase cycle;

%loop through space: Don't forget to initialize the parallel processing
%toolbox workers using 'matlabpool open N' (for N workers, 12 max).

%if you do not have the parallel computing toolbox, uncomment the first
%for loop and delete "parfor X=1:length(x)"
parfor X=1:length(x)
    for m=1:length(ph1)
        disp(['Executing X-position ' num2str(X) ...
            '; Phase cycle position ' num2str(m) ' of ' num2str(length(ph1)) '!!' ]);
        out_posx_rpc{X}{m}=sim_sLASER_shaped_Ref1(Bfield,sys,te,rfPulse,refTp,x(X),Gx,ph1(m),ph2(m),flipAngle,centreFreq);
%                            sim_sLASER_shaped_Ref1(Bfield,sys,te,RF,       tp,  dx, Gx,ph1,   ph2,  centreFreq)
    end
end

%calculate the average density matrix (Doing this inside a separate for
%loop because I couldn't figure out how to do this inside the parfor loop):
for X=1:length(x)
    for m=1:length(ph1)
        d{m}=sim_dAdd(d{m},out_posx_rpc{X}{m});
    end
end

% %Initialize structures:
out_posy_rpc =cell(length(x),length(ph1));
out=struct([]);

%Now loop through y direction (second refoc pulse only);
%for Y=1:length(y) %Use this if you don't have the MATLAB parallel processing toolbox
parfor Y=1:length(y) %Use this if you do have the MATLAB parallel processing toolbox
    for m=1:length(ph1)
        disp(['Executing Y-position ' num2str(Y) ...
            '; Phase cycle position ' num2str(m) ' of ' num2str(length(ph1)) '!!' ]);
        out_posy_rpc{Y}{m}=sim_sLASER_shaped_Ref2(d{m},n,sw,Bfield,lw,sys,te,rfPulse,refTp,y(Y),Gy,ph3(m),ph4(m),flipAngle,centreFreq);
%                            sim_sLASER_shaped_Ref2(d,   n,sw,Bfield,linewidth,sys,te,RF,       tp, dy,  Gy,ph3,    ph4,  centreFreq)
    end
end

%Now combine the outputs;  Again, doing this inside a separate for loop
%becuase I can't figure out how to do this inside the parfor loop:
for Y=1:length(y)
    for m=1:length(ph1)
        out=op_addScans(out,out_posy_rpc{Y}{m},xor(ph2(m),ph4(m)));
    end
end


%For consistent scaling across different shaped simulations, we need to :
%1.  Scale down by the total number of simulations run (since these were
%    all added together.
numSims=(nX*nY*length(ph1));
out=op_ampScale(out,1/numSims);

%2.  Scale by the total size of the simulated region, relative to the size
%    of the voxel.
voxRatio=(thkX*thkY)/(fovX*fovY);
out=op_ampScale(out,1/voxRatio);

%Plot output
figure(1),plot(out.ppm,out.specs*exp(1i*-0*pi/180)),xlim([0 4.5])

%time_sim=toc




%%%%%NESTED FUNCTIONS BELOW%%%%%%%%%

%% Simulate in X-direction only
function d = sim_sLASER_shaped_Ref1(Bfield,sys,te,RF,tp,dx,Gx,ph1,ph2,flipAngle,centreFreq)

if nargin<11
    centreFreq=2.3;
    if nargin<10
        flipAngle=180;
    end
end

%Check if this is a gradient modulated pulse.  If so, set Gx equal to zero:
if RF.isGM
    %Scale the GM waveform by the factor Gx and then set Gx equal to zero:
    RF=rf_scaleGrad(RF,Gx);
    Gx=0;
end

if (te/4)<(tp/1000)
    error('ERROR: the duration of the refocusing pulse cannot be longer than a quarter of the echo time! ABORTING!!');
end

%initialize evolution times
tau1=(te/4-tp)/2;
tau2=te/4-tp;

%Set water to centre
for k=1:length(sys)
    sys(k).shifts=sys(k).shifts-centreFreq;
end

%Calculate Hamiltonian matrices and starting density matrix.
[H,d]=sim_Hamiltonian(sys,Bfield);

%BEGIN sLASER PULSE SEQUENCE************
d=sim_excite(d,H,'x');                                  %EXCITE instantaneously
d=sim_evolve(d,H,tau1/1000);                            %Evolve by tau1
d=sim_shapedRF(d,H,RF,tp,flipAngle,ph1,dx,Gx);          %1st shaped 180 degree adiabatic refocusing pulse along X gradient
d=sim_evolve(d,H,tau2/1000);                            %Evolve by tau2
d=sim_shapedRF(d,H,RF,tp,flipAngle,ph2,dx,Gx);          %2nd shaped 180 degree adiabatic refocusing pulse along X gradient
d=sim_evolve(d,H,tau2/1000);                            %Evolve by tau2

end

%% Simulate in Y-direction only
function out = sim_sLASER_shaped_Ref2(d,n,sw,Bfield,linewidth,sys,te,RF,tp,dy,Gy,ph3,ph4,flipAngle,centreFreq)

if nargin<15
    centreFreq=2.3;
    if nargin<14
        flipAngle=180;
    end
end

%Check if this is a gradient modulated pulse.  If so, set Gy equal to zero:
if RF.isGM
    %Scale the GM waveform by the factor Gy and then set Gy equal to zero:
    RF=rf_scaleGrad(RF,Gy);
    Gy=0;
end

if (te/4)<(tp/1000)
    error('ERROR: the duration of the refocusing pulse cannot be longer than a quarter of the echo time! ABORTING!!');
end

%initialize evolution times
tau1=(te/4-tp)/2;
tau2=te/4-tp;

%Set water to centre
for k=1:length(sys)
    sys(k).shifts=sys(k).shifts-centreFreq;
end

%Calculate Hamiltonian matrices and starting density matrix.
[H]=sim_Hamiltonian(sys,Bfield);

%BEGIN sLASER PULSE SEQUENCE************
d=sim_shapedRF(d,H,RF,tp,flipAngle,ph3,dy,Gy);          %3rd shaped 180 degree adiabatic refocusing pulse along Y gradient
d=sim_evolve(d,H,tau2/1000);                            %Evolve by tau2
d=sim_shapedRF(d,H,RF,tp,flipAngle,ph4,dy,Gy);          %4th shaped 180 degree adiabatic refocusing pulse along Y gradient
d=sim_evolve(d,H,tau1/1000);                            %Evolve by tau1

[out,~]=sim_readout(d,H,n,sw,linewidth,90);      %Readout along +y axis (90 degree phase);
%END PULSE SEQUENCE**************

%Correct the ppm scale:
out.ppm=out.ppm-(4.65-centreFreq);

%Fill in structure header fields:
out.seq='semi-LASER';
out.te=te;
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
