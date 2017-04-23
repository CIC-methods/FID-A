function [outON,outOFF,outON_posxy,outOFF_posxy,outDIFF_posxy,outDIFF]=sim_runMegaPressShaped(spinSys)
% sim_runMegaPressShaped.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% This script is run simply by editing the input parameters and then
% clicking "Run".
% 
% DESCRIPTION:
% This script simulates a MEGA-PRESS experiment with fully shaped editing 
% and refocusing pulses.  Phase cycling of both the editing and refocusing
% pulses is performed.  Furthermore, simulations are run at various
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

% ************INPUT PARAMETERS**********************************
refocWaveform='mao400.pta'; %name of refocusing pulse waveform.
editWaveform='gaussEdit256.pta'; %name of editing pulse waveform.
editOnFreq=1.88; %freqeucny of edit on pulse[ppm]
editOffFreq=7.5; %frequency of edit off pulse[ppm]
refTp=5.2; %duration of refocusing pulses[ms]
editTp=17; %duration of editing pulses[ms]
Npts=2048; %number of spectral points
sw=2000; %spectral width [Hz]
Bfield=2.89; %magnetic field strength [Tesla]
lw=2; %linewidth of the output spectrum [Hz]
thkX=3; %slice thickness of x refocusing pulse [cm]
thkY=3; %slice thickness of y refocusing pulse [cm]
%x=linspace(-2.5,2.5,48); %X positions to simulate [cm]
%y=linspace(-2.5,2.5,48); %y positions to simulate [cm]
x=0;
y=0;
taus=[5.2,16.1,17.9,16.1,12.7]; %timing of the pulse sequence [ms]
%spinSys='NAA_acetyl'; %spin system to simulate
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
systems=load('spinSystems');
sys=eval(['systems.sys' spinSys]);

%Resample refocusing and editing RF pulses down to 100 pts to reduce
%computational workload
refRF=rf_resample(refRF,100);
editRF=rf_resample(editRF,100);

%This is the step where the editing pulse waveform (initially a pulse with 
%zero-frequency) is frequency shifted to produce and edit-on and an
%edit-off pulse;
editRFon=rf_freqshift(editRF,editTp,(centreFreq-editOnFreq)*Bfield*gamma/1e6);
editRFoff=rf_freqshift(editRF,editTp,(centreFreq-editOffFreq)*Bfield*gamma/1e6);


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
outDIFF=struct([]);

length(y)
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
                        outON_posxy_epc_rpc{X}{Y}{EP1}{EP2}{RP1}{RP2}=sim_megapress_shaped(Npts,sw,Bfield,lw,taus,sys,...
                            editRFon,editTp,editPhCyc1(EP1),editPhCyc2(EP2),...
                            refRF,refTp,Gx,Gy,x(X),y(Y),refPhCyc1(RP1),refPhCyc2(RP2));
                        outOFF_posxy_epc_rpc{X}{Y}{EP1}{EP2}{RP1}{RP2}=sim_megapress_shaped(Npts,sw,Bfield,lw,taus,sys,...
                            editRFoff,editTp,editPhCyc1(EP1),editPhCyc2(EP2),...
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
                    outDIFF_posxy{X}{Y}=op_subtractScans(outON_posxy{X}{Y},outOFF_posxy{X}{Y});
                else
                    outON_posxy{X}{Y}=op_addScans(outON_posxy{X}{Y},outON_posxy_epc{X}{Y}{EP1}{EP2});
                    outOFF_posxy{X}{Y}=op_addScans(outOFF_posxy{X}{Y},outOFF_posxy_epc{X}{Y}{EP1}{EP2});
                    outDIFF_posxy{X}{Y}=op_subtractScans(outON_posxy{X}{Y},outOFF_posxy{X}{Y});
                end
                
            end %end of 1st editing phase cycle loop.
        end %end of 2nd editing phase cycle loop.
        
        
        outON=op_addScans(outON,outON_posxy{X}{Y});
        outOFF=op_addScans(outOFF,outOFF_posxy{X}{Y});
        
        
        
    end %end of spatial loop (parfor) in y direction.
end %end of spatial loop (parfor) in x direction.

outDIFF=op_subtractScans(outON,outOFF);

save([spinSys '_MP_sim_68ms_v2'],'outON','outOFF','outON_posxy','outOFF_posxy','outDIFF_posxy');    
% figure 
% hold
% for n=1:length(x)
% plot(outON_posxy{n}{1}.ppm,outON_posxy{n}{1}.specs+5*n);
% plot(outOFF_posxy{n}{1}.ppm,outOFF_posxy{n}{1}.specs+5*n);
% end



       