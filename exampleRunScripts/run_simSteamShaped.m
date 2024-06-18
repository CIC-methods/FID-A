% run_simSteamShaped.m
% Jamie Near, Sunnybrook Research Institute 2022.
% 
% USAGE:
% [out,out_posxy]=run_simSteamShaped(spinSys);
% 
% DESCRIPTION:
% This script simulates a STEAM experiment with fully shaped RF 
% pulses. Coherence order filtering is employed to remove unwanted
% coherences.  Simulations are run at various locations in space to account 
% for the within-voxel spatial variation of the metabolite signal.  Summation 
% across spatial positions is performed. The MATLAB parallel computing toolbox 
% (parfor loop) was used to accelerate the simulations.  Acceleration 
% is currently performed in the direction of the slice selective pulse along
% the x-direction, but this can be changed.  Up to a factor of 12 acceleration
% can be achieved using this approach.  To enable the use of the MATLAB
% parallel computing toolbox, initialize the multiple worked nodes using
% "matlabpool size X" where "X" is the number of available processing
% nodes.  If the parallel processing toolbox is not available, then replace
% the "parfor" loop with a "for" loop.
% 
% INPUTS:
% To run this script, there is technically only one input argument:
% spinSys           = spin system to simulate 
%
% However, the user should also edit the following parameters as 
% desired before running the function:
% rfWaveform        = name of refocusing pulse waveform.
% Tp                = duration of refocusing pulses[ms]
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
% out_posxy         = Simulation results, spatially resolved.
% out               = Simulation results, summed over all space.

function [out,out_posxy]=run_simSteamShaped(spinSys)
% ************INPUT PARAMETERS**********************************
rfWaveform='sampleExcPulse.pta'; %name of RF pulse waveform.
Tp=3; %duration of refocusing pulses[ms]
flipAngle=90; %RF pulse flip angle [degrees]
Npts=8192; %number of spectral points
sw=6000; %spectral width [Hz]
Bfield=6.98; %magnetic field strength [Tesla]
lw=1; %linewidth of the output spectrum [Hz]
thkX=2.5; %slice thickness of x RF pulse [cm]
thkY=2.5; %slice thickness of y RF pulse [cm]
fovX=5; %size of the full simulation Field of View in the x-direction [cm]
fovY=5; %size of the full simulation Field of View in the y-direction [cm]
nX=9; %Number of grid points to simulate in the x-direction
nY=9; %Number of grid points to simulate in the y-direction
tau1=135; %TE for STEAM sequence [ms]
tau2=32; %TM for STEAM sequence [ms]
centreFreq=2.7; %Centre frequency of refocusing pulses [ppm]
% ************END OF INPUT PARAMETERS**********************************

%set up spatial grid
x=linspace(-fovX/2,fovX/2,nX); %X positions to simulate [cm]
y=linspace(-fovY/2,fovY/2,nY); %y positions to simulate [cm]

%Load RF waveform
RF=io_loadRFwaveform(rfWaveform,'exc',0);


gamma=42577000; %gyromagnetic ratio

%Load spin systems
load spinSystems
sys=eval(['sys' spinSys]);

%If length of RF pulse is >200 pts, resample to 100 pts to reduce
%computational workload
if size(RF.waveform,1)>200
    RF=rf_resample(RF,100);
end

Gx=(RF.tbw/(Tp/1000))/(gamma*thkX/10000); %[G/cm]
Gy=(RF.tbw/(Tp/1000))/(gamma*thkY/10000); %[G/cm]

%Initialize structures:
out_posxy=cell(length(x),length(y));
out=struct([]);

%loop through space: Don't forget to initialize the parallel processing
%toolbox workers using 'matlabpool open N' (for N workers, 12 max).

%for X=1:length(x);
parfor X=1:length(x)
    for Y=1:length(y)
                disp(['Executing X-position ' num2str(X) ' of ' num2str(length(x)) ', '...
                    'Y-position ' num2str(Y) ' of ' num2str(length(y))])
                out_posxy{X}{Y}=sim_steam_shaped(Npts,sw,Bfield,lw,sys,tau1,tau2,...
                    RF,Tp,x(X),y(Y),Gx,Gy,flipAngle,centreFreq);
         
        out=op_addScans(out,out_posxy{X}{Y});
        
    end %end of spatial loop (parfor) in y direction.
end %end of spatial loop (parfor) in x direction.
        
%For consistent scaling across different shaped simulations, we need to :
%1.  Scale down by the total number of simulations run (since these were
%    all added together.
numSims=(nX*nY);
out=op_ampScale(out,1/numSims);

%2.  Scale by the total size of the simulated region, relative to the size
%    of the voxel.
voxRatio=(thkX*thkY)/(fovX*fovY);
out=op_ampScale(out,1/voxRatio);


%Uncomment this part to view the results of the spatially resolved simulation:
%  if length(x)>1 && length(y)>1
%      ppmmin=1.0;
%      ppmmax=1.6;
%      sim_make2DSimPlot(out_posxy,ppmmin,ppmmax);
%  else
%      plot(out.ppm,out.specs);
%  end
% op_plotspec(out)



       
