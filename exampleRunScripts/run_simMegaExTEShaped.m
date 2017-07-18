% run_simMegaExTEShaped.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% This script is run simply by editing the input parameters and then
% clicking "Run".
% 
% DESCRIPTION:
% This script simulates an ExTE-MEGA-PRESS experiment with fully shaped editing 
% and refocusing pulses.  Phase cycling of both the editing and refocusing
% pulses is performed.  Simulations are run at various
% locations in space to account for the within-voxel spatial variation of
% the GABA signal.  Summation across phase cycles and spatial positions is
% performed.  As a result of the phase cycling and spatially resolved simulations, 
% this code takes a long time to run.  Therefore, the MATLAB parallel computing
% toolbox (parfor loop) was used to accelerate the siumulations.  Accelration 
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
% editWaveform      = name of editing pulse waveform.
% editOnFreq        = freqeucny of edit on pulse[ppm]
% editOffFreq       = frequency of edit off pulse[ppm]
% refTp             = duration of refocusing pulses[ms]
% editTp            = duration of editing pulses[ms]
% Bfield            = Magnetic field strength in [T]
% Npts              = number of spectral points
% sw                = spectral width [Hz]
% Bfield            = magnetic field strength [Tesla]
% lw                = linewidth of the output spectrum [Hz]
% thkX              = slice thickness of x refocusing pulse [cm]
% thkY              = slice thickness of y refocusing pulse [cm]
% x                 = vector of X positions to simulate [cm]
% y                 = vector of y positions to simulate [cm]
% taus              = vector of pulse sequence timings  [ms]
% spinSys           = spin system to simulate 
% editPhCyc1        = vector of phase cycling steps for 1st editing pulse [degrees]
% editPhCyc2        = vector of phase cycling steps for 2nd editing pulse [degrees]
% refPhCyc1         = vector of phase cycling steps for 1st refocusing pulse [degrees]
% refPhCyc2         = vector of phase cycling steps for 2nd refocusing pulse [degrees]
%
% OUTPUTS:
% outON_posxy       = Simulated ExTE-MEGA-PRESS edit-ON spectrum, spatially resolved. 
% outOFF_posxy      = Simulated ExTE-MEGA-PRESS edit-OFF spectrum, spatially resolved.
% outDIFF_posxy     = Simulated ExTE-MEGA-PRESS difference spectrum, spatially resolved.
% outON             = Simulated ExTE-MEGA-PRESS edit-ON spectrum, summed over
%                     all positions.
% outOFF            = Simulated ExTE-MEGA-PRESS edit-OFF spectrum, summed over
%                     all positions.
% outDIFF           = Simulated ExTE-MEGA-PRESS difference spectrum, summed over
%                     all positions.

% ************INPUT PARAMETERS**********************************
refocWaveform='sampleRefocPulse.pta'; %name of refocusing pulse waveform.
editWaveform='sampleEditPulse.pta'; %name of editing pulse waveform.
editOnFreq=1.88; %freqeucny of edit on pulse[ppm]
refTp=5.2; %duration of refocusing pulses[ms]
editTp=14; %duration of editing pulses[ms]
Npts=2048; %number of spectral points
sw=2000; %spectral width [Hz]
lw=2; %linewidth of the output spectrum [Hz]
Bfield=2.89; %Magnetic field strength in [T]
thkX=3.5; %slice thickness of x refocusing pulse [cm]
thkY=3.5; %slice thickness of y refocusing pulse [cm]
x=linspace(-2.0125,2.0125,12); %X positions to simulate [cm]
y=linspace(-2.0125,2.0125,12); %y positions to simulate [cm]
taus_inv=...  %Timing for J-Inverted scan
        [4.9,...    %time from excitation to 1st refoc pulse [ms]
        67.4885,... %time from 1st refoc pulse to 1st editing pulse [ms]
        33.5115,... %time from 1st editing pulse to 2nd refoc pulse [ms]
        33.4885,... %time from 2nd refoc pulse to 2nd editing pulse [ms]
        62.6115];   %time from 2nd editing pulse to ADC onset [ms]
taus_ref=... %Timing for J-refocused scan
        [4.9,...    %time from excitation to 1st refoc pulse [ms]
        50.4885,... %time from 1st refoc pulse to 1st editing pulse [ms]
        50.5115,... %time from 1st editing pulse to 2nd refoc pulse [ms]
        50.4885,... %time from 2nd refoc pulse to 2nd editing pulse [ms]
        45.6115];   %time from 2nd editing pulse to ADC onset [ms]
spinSys='GABA'; %spin system to simulate
centreFreq=3.0; %Centre frequency of MR spectrum [ppm]
editPhCyc1=[0 90]; %phase cycling steps for 1st editing pulse [degrees]
editPhCyc2=[0 90]; %phase cycling steps for 2nd editing pulse [degrees]
refPhCyc1=[0,90]; %phase cycling steps for 1st refocusing pulse [degrees]
refPhCyc2=[0,90]; %phase cycling steps for 2nd refocusing pulse [degrees]
% ************END OF INPUT PARAMETERS**********************************

%Load RF waveforms
refRF=io_loadRFwaveform(refocWaveform,'ref',0);
editRF=io_loadRFwaveform(editWaveform,'inv',0);

gamma=42577000; %gyromagnetic ratio

%Load spin systems
load spinSystems
if strcmp(spinSys,'MM');
    sys=sysGABA;
    sys.shifts(3)=1.7;
    sys.shifts(4)=1.7;
else
    sys=eval(['sys' spinSys]);
end
    
%Resample refocusing RF pulse from 400 pts to 100 pts to reduce
%computational workload
refRF=rf_resample(refRF,100);

%This is the step where the editing pulse waveform (initially a pulse with 
%zero-frequency) is frequency shifted to produce and edit-on and an
%edit-off pulse;
editRFon=rf_freqshift(editRF,editTp,(centreFreq-editOnFreq)*Bfield*gamma/1e6);

Gx=(refRF.tbw/(refTp/1000))/(gamma*thkX/10000); %[G/cm]
Gy=(refRF.tbw/(refTp/1000))/(gamma*thkY/10000); %[G/cm]

[DX,DY]=meshgrid(x,y);

%n=1;
%totalIters=length(x)*length(y)*length(editPhCyc1)*length(editPhCyc2)*length(refPhCyc1)*length(refPhCyc2);

%Initialize structures:
outON_posxy_epc_rpc=cell(length(x),length(y),length(editPhCyc1),length(editPhCyc2),length(refPhCyc1),length(refPhCyc2));
outOFF_posxy_epc_rpc=cell(length(x),length(y),length(editPhCyc1),length(editPhCyc2),length(refPhCyc1),length(refPhCyc2));
outON_posxy_epc=cell(length(x),length(y),length(editPhCyc1),length(editPhCyc2));
outOFF_posxy_epc=cell(length(x),length(y),length(editPhCyc1),length(editPhCyc2));
outON_posxy=cell(length(x),length(y));
outOFF_posxy=cell(length(x),length(y));
outDIFF_posxy=cell(length(x),length(y));
outON=struct([]);
outOFF=struct([]);


%loop through space: Don't forget to initialize the parallel processing
%toolbox workers using 'matlabpool open N' (for N workers, 12 max).

for X=1:length(x);
%parfor X=1:length(x);
    for Y=1:length(y);
        for EP1=1:length(editPhCyc1)
            for EP2=1:length(editPhCyc2)
                for RP1=1:length(refPhCyc1)
                    for RP2=1:length(refPhCyc2)
                        disp(['Executing X-position ' num2str(X) ' of ' num2str(length(x)) ', '...
                            'Y-position ' num2str(Y) ' of ' num2str(length(y)) ', '...
                            'First Edit phase cycle ' num2str(EP1) ' of ' num2str(length(editPhCyc1)) ', '...
                            'Second Edit phase cycle ' num2str(EP2) ' of ' num2str(length(editPhCyc2)) ', '...
                            'First Refoc phase cycle ' num2str(RP1) ' of ' num2str(length(refPhCyc1)) ', '...
                            'Second Refoc phase cycle ' num2str(RP2) ' of ' num2str(length(refPhCyc2)) '!!!']); 
                        outON_posxy_epc_rpc{X}{Y}{EP1}{EP2}{RP1}{RP2}=sim_megapress_shaped(Npts,sw,Bfield,lw,taus_ref,sys,...
                            editRFon,editTp,editPhCyc1(EP1),editPhCyc2(EP2),...
                            refRF,refTp,Gx,Gy,x(X),y(Y),refPhCyc1(RP1),refPhCyc2(RP2));
                        outOFF_posxy_epc_rpc{X}{Y}{EP1}{EP2}{RP1}{RP2}=sim_megapress_shaped(Npts,sw,Bfield,lw,taus_inv,sys,...
                            editRFon,editTp,editPhCyc1(EP1),editPhCyc2(EP2),...
                            refRF,refTp,Gx,Gy,x(X),y(Y),refPhCyc1(RP1),refPhCyc2(RP2));
                    
                        if RP1==1 && RP2==1
                            outON_posxy_epc{X}{Y}{EP1}{EP2}=outON_posxy_epc_rpc{X}{Y}{EP1}{EP2}{RP1}{RP2};
                            outOFF_posxy_epc{X}{Y}{EP1}{EP2}=outOFF_posxy_epc_rpc{X}{Y}{EP1}{EP2}{RP1}{RP2};
                        else
                            outON_posxy_epc{X}{Y}{EP1}{EP2}=op_addScans(outON_posxy_epc{X}{Y}{EP1}{EP2},outON_posxy_epc_rpc{X}{Y}{EP1}{EP2}{RP1}{RP2},xor(RP1==length(refPhCyc1),RP2==length(refPhCyc2)));
                            outOFF_posxy_epc{X}{Y}{EP1}{EP2}=op_addScans(outOFF_posxy_epc{X}{Y}{EP1}{EP2},outOFF_posxy_epc_rpc{X}{Y}{EP1}{EP2}{RP1}{RP2},xor(RP1==length(refPhCyc1),RP2==length(refPhCyc2)));
                        end
                    end %end of 1st refocusing phase cycle loop
                end %end of 2nd refocusing phase cycle loop.
                
                if EP1==1 && EP2==1
                    outON_posxy{X}{Y}=outON_posxy_epc{X}{Y}{EP1}{EP2};
                    outOFF_posxy{X}{Y}=outOFF_posxy_epc{X}{Y}{EP1}{EP2};
                else
                    outON_posxy{X}{Y}=op_addScans(outON_posxy{X}{Y},outON_posxy_epc{X}{Y}{EP1}{EP2});
                    outOFF_posxy{X}{Y}=op_addScans(outOFF_posxy{X}{Y},outOFF_posxy_epc{X}{Y}{EP1}{EP2});
                end
                outDIFF_posxy{X}{Y}=op_subtractScans(outON_posxy{X}{Y},outOFF_posxy{X}{Y});
            end %end of 1st editing phase cycle loop.
        end %end of 2nd editing phase cycle loop.
        
        
        outON=op_addScans(outON,outON_posxy{X}{Y});
        outOFF=op_addScans(outOFF,outOFF_posxy{X}{Y});
        
        
    end %end of spatial loop (parfor) in y direction.
end %end of spatial loop (parfor) in x direction.

outDIFF=op_subtractScans(outON,outOFF);
        
figure 
sim_make2DSimPlot(outON_posxy,2.75,3.25);
figure
sim_make2DSimPlot(outOFF_posxy,2.75,3.25);
figure
sim_make2DSimPlot(outDIFF_posxy,2.75,3.25);



       