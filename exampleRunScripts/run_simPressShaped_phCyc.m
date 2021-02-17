% 2021 EDIT: This function is now deprecated. The current version of the shaped
% PRESS simulation is run_simPressShaped.m. The current version employs coherence
% selection, resulting in a 4x simulation speed increase.

% run_simPressShaped_phCyc.m
% Jamie Near, McGill University 2015.
% 
% USAGE:
% This script is run simply by editing the input parameters and then
% clicking "Run".
% 
% DESCRIPTION:
%
% This script simulates a PRESS experiment with fully shaped refocusing 
% pulses.  Phase cycling of refocusing pulses is performed.  Furthermore, 
% simulations are run at various locations in space to account for the 
% within-voxel spatial variation of the metabolite signal.  Summation 
% across phase cycles and spatial positions is performed.  As a result of 
% the phase cycling and spatially resolved simulations, this code takes a 
% long time to run.  Therefore, the MATLAB parallel computing toolbox 
% (parfor loop) was used to accelerate the siumulations.  Acceleration 
% is currently performed in the direction of the slice selective pulse along
% the x-direction, but this can be changed.  Up to a factor of 12 acceleration
% can be achieved using this approach.  To enable the use of the MATLAB
% parallel computing toolbox, initialize the multiple worked nodes using
% "matlabpool size X" where "X" is the number of available processing
% nodes.  If the parallel processing toolbox is not available, then replace
% the "parfor" loop with a "for" loop.
% 
% INPUTS:
% To run this script, edit the parameters below as desired and then click
% "run":
% refocWaveform     = name of refocusing pulse waveform.
% refTp             = duration of refocusing pulses[ms]
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
% spinSys           = spin system to simulate 
% refPhCyc1         = vector of phase cycling steps for 1st refocusing pulse [degrees]
% refPhCyc2         = vector of phase cycling steps for 2nd refocusing pulse [degrees]
%
% OUTPUTS:
% out_posxy         = Simulation results, spatially resolved.
% out               = Simulation results, summed over all space.

% ************INPUT PARAMETERS**********************************
refocWaveform='sampleRefocPulse.pta'; %name of refocusing pulse waveform.
refTp=3.5; %duration of refocusing pulses[ms]
Npts=2048; %number of spectral points
sw=2000; %spectral width [Hz]
Bfield=3; %magnetic field strength [Tesla]
lw=2; %linewidth of the output spectrum [Hz]
thkX=1.66; %slice thickness of x refocusing pulse [cm]
thkY=1.66; %slice thickness of y refocusing pulse [cm]
fovX=2.4; %size of the full simulation Field of View in the x-direction [cm]
fovY=2.4; %size of the full simulation Field of View in the y-direction [cm]
nX=32; %Number of grid points to simulate in the x-direction
nY=32; %Number of grid points to simulate in the y-direction
tau1=30; %TE1 for first spin echo [ms]
tau2=105; %TE2 for second spin echo [ms]
spinSys='Lac'; %spin system to simulate
refPhCyc1=[0,90]; %phase cycling steps for 1st refocusing pulse [degrees]
refPhCyc2=[0,90]; %phase cycling steps for 2nd refocusing pulse [degrees]
% ************END OF INPUT PARAMETERS**********************************

%set up spatial grid
x=linspace(-fovX/2,fovX/2,nX); %X positions to simulate [cm]
y=linspace(-fovY/2,fovY/2,nY); %y positions to simulate [cm]

%Load RF waveform
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
out_posxy_rpc=cell(length(x),length(y),length(refPhCyc1),length(refPhCyc2));
out_posxy=cell(length(x),length(y));
out=struct([]);

%loop through space: Don't forget to initialize the parallel processing
%toolbox workers using 'matlabpool open N' (for N workers, 12 max).

%for X=1:length(x);
parfor X=1:length(x)
    for Y=1:length(y)
        for RP1=1:length(refPhCyc1)
            for RP2=1:length(refPhCyc2)
                disp(['Executing X-position ' num2str(X) ' of ' num2str(length(x)) ', '...
                    'Y-position ' num2str(Y) ' of ' num2str(length(y)) ', '...
                    'First Refoc phase cycle ' num2str(RP1) ' of ' num2str(length(refPhCyc1)) ', '...
                    'Second Refoc phase cycle ' num2str(RP2) ' of ' num2str(length(refPhCyc2)) '!!!']);
                out_posxy_rpc{X}{Y}{RP1}{RP2}=sim_press_shaped(Npts,sw,Bfield,lw,sys,tau1,tau2,...
                    refRF,refTp,x(X),y(Y),Gx,Gy,refPhCyc1(RP1),refPhCyc2(RP2));
                
                if RP1==1 && RP2==1
                    out_posxy{X}{Y}=out_posxy_rpc{X}{Y}{RP1}{RP2};
                else
                    out_posxy{X}{Y}=op_addScans(out_posxy{X}{Y},out_posxy_rpc{X}{Y}{RP1}{RP2},xor(RP1==length(refPhCyc1),RP2==length(refPhCyc2)));
                end
            end %end of 1st refocusing phase cycle loop
        end %end of 2nd refocusing phase cycle loop.
        
        out=op_addScans(out,out_posxy{X}{Y});
        
    end %end of spatial loop (parfor) in y direction.
end %end of spatial loop (parfor) in x direction.
        
%For consistent scaling across different shaped simulations, we need to :
%1.  Scale down by the total number of simulations run (since these were
%    all added together.
numSims=(nX*nY*length(refPhCyc1)*length(refPhCyc2));
out=op_ampScale(out,1/numSims);

%2.  Scale by the total size of the simulated region, relative to the size
%    of the voxel.
voxRatio=(thkX*thkY)/(fovX*fovY);
out=op_ampScale(out,1/voxRatio);


%Uncomment this part to view the results of the spatially resolved simulation:
% figure 
% if length(x)>1 && length(y)>1
%     ppmmin=1;
%     ppmmax=1.6;
%     sim_make2DSimPlot(out_posxy,ppmmin,ppmmax);
% else
%     plot(out.ppm,out.specs);
% end



       
