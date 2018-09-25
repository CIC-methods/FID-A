% run_simMegaPressShapedRefoc.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% This script is run simply by editing the input parameters and then
% clicking "Run".
% 
% DESCRIPTION:
% This script simulates a MEGA-PRESS experiment with fully shaped 
% refocusing pulses.  Phase cycling of the refocusing pulses is performed 
% and simulations are run at various locations in space to account for the 
% within-voxel spatial variation of the GABA signal.  Summation across 
% spatial positions is performed.  As a result of the phase cycling and 
% spatially resolved simulations, this code takes a long time to run.  
% Therefore, the MATLAB parallel computing toolbox (parfor loop) is used 
% to accelerate the siumulations.  Acceleration is currently performed in 
% the direction of the slice selective pulse along the x-direction, but 
% this can be changed.  Up to a factor of 12 acceleration can be achieved 
% using this approach.  To enable the use of the MATLAB parallel computing 
% toolbox, initialize the multiple worked nodes using "matlabpool size X" 
% where "X" is the number of available processing nodes.  If the parallel 
% processing toolbox is not available, then replace the "parfor" loop with 
% a "for" loop.
% 
% INPUTS:
% To run this script, edit the parameters below as desired and then click
% "run":
% refocWaveform         = name of refocusing pulse waveform.
% refTp                 = duration of refocusing pulses[ms]
% Npts                  = number of spectral points
% sw                    = spectral width [Hz]
% Bfield                = magnetic field strength [Tesla]
% lw                    = linewidth of the output spectrum [Hz]
% thkX                  = slice thickness of x refocusing pulse [cm]
% thkY                  = slice thickness of y refocusing pulse [cm]
% x                     = vector of X positions to simulate [cm]
% y                     = vector of y positions to simulate [cm]
% taus                  = vector of pulse sequence timings  [ms]
% spinSys               = spin system to simulate
% editFlipON            = Cell array of edit-ON pulse flip angles for each 
%                         spin in each part of the spin system.
% editFlipOFF           = Cell array of edit-OFF pulse flip angles for each 
%                         spin in each part of the spin system.
% refPhCyc1             = vector of phase cycling steps for 1st refocusing pulse [degrees]
% refPhCyc2             = vector of phase cycling steps for 2nd refocusing pulse [degrees]
%
% OUTPUTS:
% outON_posxy       = Simulated MEGA-PRESS edit-ON spectrum, spatially resolved. 
% outOFF_posxy      = Simulated MEGA-PRESS edit-OFF spectrum, spatially resolved.
% outDIFF_posxy     = Simulated MEGA-PRESS difference spectrum, spatially resolved.
% outON             = Simulated MEGA-PRESS edit-ON spectrum, summed over
%                     all positions.
% outOFF            = Simulated MEGA-PRESS edit-OFF spectrum, summed over
%                     all positions.
% outDIFF           = Simulated MEGA-PRESS difference spectrum, summed over
%                     all positions.

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
spinSys='GABA'; %spin system to simulate
editFlipON{1}=[0 0 180 180 0 0]; %Cell array of edit-on pulse flip angles.
editFlipOFF{1}=[0 0 0 0 0 0]; %Cell array of edit-off pulse flip angles.
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

[DX,DY]=meshgrid(x,y);

%Initialize structures:
outON_posxy_rpc=cell(length(x),length(y),length(refPhCyc1),length(refPhCyc2));
outOFF_posxy_rpc=cell(length(x),length(y),length(refPhCyc1),length(refPhCyc2));
outON_posxy=cell(length(x),length(y));
outOFF_posxy=cell(length(x),length(y));
outDIFF_posxy=cell(length(x),length(y));
outON=struct([]);
outOFF=struct([]);

%loop through space: Don't forget to initialize the parallel processing
%toolbox workers using 'matlabpool open N' (for N workers, 12 max).

%for X=1:length(x);
parfor X=1:length(x);
    for Y=1:length(y);
        for RP1=1:length(refPhCyc1)
            for RP2=1:length(refPhCyc2)
                disp(['Executing X-position ' num2str(X) ' of ' num2str(length(x)) ', '...
                    'Y-position ' num2str(Y) ' of ' num2str(length(y)) ', '...
                    'First Refoc phase cycle ' num2str(RP1) ' of ' num2str(length(refPhCyc1)) ', '...
                    'Second Refoc phase cycle ' num2str(RP2) ' of ' num2str(length(refPhCyc2)) '!!!']);
                outON_posxy_rpc{X}{Y}{RP1}{RP2}=sim_megapress_shapedRefoc(Npts,sw,Bfield,lw,taus,sys,editFlipON,...
                    refRF,refTp,Gx,Gy,x(X),y(Y),refPhCyc1(RP1),refPhCyc2(RP2));
                outOFF_posxy_rpc{X}{Y}{RP1}{RP2}=sim_megapress_shapedRefoc(Npts,sw,Bfield,lw,taus,sys,editFlipOFF,...
                    refRF,refTp,Gx,Gy,x(X),y(Y),refPhCyc1(RP1),refPhCyc2(RP2));
                
                if RP1==1 && RP2==1
                    outON_posxy{X}{Y}=outON_posxy_rpc{X}{Y}{RP1}{RP2};
                    outOFF_posxy{X}{Y}=outOFF_posxy_rpc{X}{Y}{RP1}{RP2};
                else
                    outON_posxy{X}{Y}=op_addScans(outON_posxy{X}{Y},outON_posxy_rpc{X}{Y}{RP1}{RP2},xor(RP1~=length(refPhCyc1),RP2~=length(refPhCyc2)));
                    outOFF_posxy{X}{Y}=op_addScans(outOFF_posxy{X}{Y},outOFF_posxy_rpc{X}{Y}{RP1}{RP2},xor(RP1~=length(refPhCyc1),RP2~=length(refPhCyc2)));
                end
                outDIFF_posxy{X}{Y}=op_subtractScans(outON_posxy{X}{Y},outOFF_posxy{X}{Y});
            end %end of 1st refocusing phase cycle loop
        end %end of 2nd refocusing phase cycle loop
        
        outON=op_addScans(outON,outON_posxy{X}{Y});
        outOFF=op_addScans(outOFF,outOFF_posxy{X}{Y});
        
    end %end of spatial loop (parfor) in y direction.
end %end of spatial loop (parfor) in x direction.

%For consistent scaling across different shaped simulations, we need to :
%1.  Scale down by the total number of simulations run (since these were
%    all added together.
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
        
sim_make2DSimPlot(outON_posxy,2.75,3.25);
sim_make2DSimPlot(outOFF_posxy,2.75,3.25);
sim_make2DSimPlot(outDIFF_posxy,2.75,3.25);



       