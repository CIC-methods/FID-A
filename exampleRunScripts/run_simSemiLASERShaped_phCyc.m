% run_simSemiLASERShaped_phCyc.m
% Dana Goerzen (McGill University, 2019).
%  
% USAGE:
% This script runs a spatially resolved, shaped simulation of the sLASER 
% experiment, using the sim_sLASER_shaped.m function.
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

% 
% INPUTS:

n=4096; %= number of points in fid/spectrum
sw=5000; %= desired spectral width in [Hz]
Bfield=7; %= main magnetic field strength in [T]
lw=2; %= linewidth in [Hz]
load spinSystems.mat; %= spin system definition structure
sys=sysLac;
rfPulse=io_loadRFwaveform('sampleAFPpulse_HS2_R15.RF','inv'); % gradient_modulated adiabatic RF pulse shaped waveform 
refTp=3.5; %= RF pulse duration in [ms]  
flipAngle=180; %= flip angle of refocusing pulses [degrees] (Optional.  Default = 180 deg)
centreFreq=2.3; %= centre frequency of the spectrum in [ppm] (Optional.  Default = 2.3)
thkX=2; %slice thickness of x refocusing pulse [cm].  If pulse is GM, this must equal thkY.
thkY=2; %slice thickness of y refocusing pulse [cm].  If pulse is GM, this must equal thkY.
fovX=3; %size of the full simulation Field of View in the x-direction [cm]
fovY=3; %size of the full simulation Field of View in the y-direction [cm]
nX=8; %Number of grid points to simulate in the x-direction
nY=8; %Number of grid points to simulate in the y-direction
te=135;         %sLASER total echo time [ms]
ph1=[0 0 0 0];  %phase cycling scheme of first refocusing pulse (in degrees)
ph2=[0 0 90 90]; %phase cycling scheme of second refocusing pulse (in degrees)
ph3=[0 0 0 0]; %phase cycling scheme of third refocusing pulse (in degrees)
ph4=[0 90 0 90]; %phase cycling scheme of fourth refocusing pulse (in degrees)

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
    Gx=(rfPulse.tbw/(refTp/1000))/(gamma*thkX/10000); %[G/cm]
    Gy=(rfPulse.tbw/(refTp/1000))/(gamma*thkY/10000); %[G/cm]
else
    %Gradient modulated pulse
    %1.  Calculating the unitless scaling factor for the GM waveform.
    Gx=(rfPulse.tthk/(refTp/1000))/thkX;
    Gy=(rfPulse.tthk/(refTp/1000))/thkY;
end

%Initialize structures:
out_posxy_rpc=cell(length(x),length(y),length(ph1));
out_posxy=cell(length(x),length(y));
out=struct([]);

%loop through space: Don't forget to initialize the parallel processing
%toolbox workers using 'matlabpool open N' (for N workers, 12 max).

%if you do not have the parallel computing toolbox, uncomment the first
%for loop and delete "parfor X=1:length(x)"
%for X=1:length(x)
parfor X=1:length(x)
    for Y=1:length(y)
        for m=1:length(ph1)
            disp(['Executing X-position ' num2str(X) ' of ' num2str(length(x)) '; Y-position ' num2str(Y) ...
                ' of ' num2str(length(y)) '; Phase cycle position ' num2str(m) ' of ' num2str(length(ph1)) '!!' ]);
            out_posxy_rpc{X}{Y}{m}=sim_sLASER_shaped(n,sw,Bfield,lw,sys,te,...
              rfPulse,refTp,x(X),y(Y),Gx,Gy,ph1(m),ph2(m),ph3(m),ph4(m),flipAngle,centreFreq);
            if m==1 
                 out_posxy{X}{Y}=out_posxy_rpc{X}{Y}{m};
             elseif m==2 || m==3
                out_posxy{X}{Y}=op_addScans(out_posxy{X}{Y},out_posxy_rpc{X}{Y}{m},1);
            else
                out_posxy{X}{Y}=op_addScans(out_posxy{X}{Y},out_posxy_rpc{X}{Y}{m});

            end
        end
             out=op_addScans(out,out_posxy{X}{Y});
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


%Uncomment this part to view the results of the spatially resolved simulation:
% 
% if length(x)>1 && length(y)>1
%     ppmmin=0;
%     ppmmax=6;
%     sim_make2DSimPlot(out_posxy,ppmmin,ppmmax);
% else
%     plot(out.ppm,out.specs);
% end
% op_plotspec(out)
