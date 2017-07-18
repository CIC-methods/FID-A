% run_simMegaSpecialShaped.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% This script is run simply by editing the input parameters and then
% clicking "Run".
% 
% DESCRIPTION:
% This script simulates a MEGA-SPECIAL experiment with fully shaped editing 
% and refocusing pulses.  Phase cycling of both the editing and refocusing
% pulses is performed.  Furthermore, simulations are run at various
% locations in space to account for the within-voxel spatial variation of
% the GABA signal.  Summation across phase cycles and spatial positions is
% performed.  As a result of the phase cycling and spatially resolved simulations, 
% this code takes a long time to run.  Therefore, the MATLAB parallel computing
% toolbox (parfor loop) was used to accelerate the siumulations.  Accelration 
% is performed in the direction of the slice selective pulse (along
% the x-direction).  Up to a factor of 12 acceleration can be achieved using 
% this approach.  To enable the use of the MATLAB parallel computing toolbox, 
% initialize the multiple worked nodes using "matlabpool size X" where "X" 
% is the number of available processing nodes.  If the parallel processing 
% toolbox is not available, then replace the "parfor" loop with a "for" loop.
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
% x                 = vector of X positions to simulate [cm]
% taus              = vector of pulse sequence timings  [ms]
% spinSys           = spin system to simulate 
% editPhCyc1        = vector of phase cycling steps for 1st editing pulse [degrees]
% editPhCyc2        = vector of phase cycling steps for 2nd editing pulse [degrees]
% refPhCyc          = vector of phase cycling steps for 1st refocusing pulse [degrees]
%
% OUTPUTS:
% outON_posx        = Simulated MEGA-SPECIAL edit-ON spectrum, spatially resolved. 
% outOFF_posx       = Simulated MEGA-SPECIAL edit-OFF spectrum, spatially resolved.
% outDIFF_posx      = Simulated MEGA-SPECIAL difference spectrum, spatially resolved.
% outON             = Simulated MEGA-SPECIAL edit-ON spectrum, summed over
%                     all positions.
% outOFF            = Simulated MEGA-SPECIAL edit-OFF spectrum, summed over
%                     all positions.
% outDIFF           = Simulated MEGA-SPECIAL difference spectrum, summed over
%                     all positions.

% ************INPUT PARAMETERS**********************************
refocWaveform='sampleRefocPulse.pta'; %name of refocusing pulse waveform.
editWaveform='sampleEditPulse.pta'; %name of editing pulse waveform.
editOnFreq=1.88; %freqeucny of edit on pulse[ppm]
editOffFreq=7.5; %frequency of edit off pulse[ppm]
refTp=5; %duration of refocusing pulses[ms]
editTp=14; %duration of editing pulses[ms]
Npts=2048; %number of spectral points
sw=2000; %spectral width [Hz]
Bfield=3; %magnetic field strength [T]
lw=2; %linewidth of the output spectrum [Hz]
thkX=3; %slice thickness of x refocusing pulse [cm]
x=linspace(-2.2,2.2,12); %X positions to simulate [cm]
taus=[17,...    %Time from excitation pulse to 1st editing pulse [ms]
    17,...      %Time from 1st editing pulse to refoc pulse [ms]
    17,...      %Time from refoc pulse to 2nd editing pulse [ms]
    17];...     %Time from 2nd editing pulse to ADC onset [ms]
spinSys='GABA'; %spin system to simulate
centreFreq=3.0; %Centre frequency of the MR spectrum [ppm]
editPhCyc1=[0 90]; %phase cycling steps for 1st editing pulse [degrees]
editPhCyc2=[0 90 180 270]; %phase cycling steps for 2nd editing pulse [degrees]
refPhCyc=[0,90]; %phase cycling steps for 1st refocusing pulse [degrees]
% ************END OF INPUT PARAMETERS**********************************

%Load RF waveforms
refRF=io_loadRFwaveform(refocWaveform,'ref',0);
editRF=io_loadRFwaveform(editWaveform,'inv',0);

gamma=42577000; %gyromagnetic ratio

%Load spin systems
load spinSystems
sys=eval(['sys' spinSys]);

%Resample refocusing RF pulse from 400 pts to 100 pts to reduce
%computational workload
refRF=rf_resample(refRF,100);

%This is the step where the editing pulse waveform (initially a pulse with 
%zero-frequency) is frequency shifted to produce and edit-on and an
%edit-off pulse;
editRFon=rf_freqshift(editRF,editTp,(centreFreq-editOnFreq)*Bfield*gamma/1e6);
editRFoff=rf_freqshift(editRF,editTp,(centreFreq-editOffFreq)*Bfield*gamma/1e6);


Gx=(refRF.tbw/(refTp/1000))/(gamma*thkX/10000); %[G/cm]

%n=1;
%totalIters=length(x)*length(y)*length(editPhCyc1)*length(editPhCyc2)*length(refPhCyc1)*length(refPhCyc2);

%Initialize structures:
outON_posx_epc_rpc=cell(length(x),length(editPhCyc1),length(editPhCyc2),length(refPhCyc));
outOFF_posx_epc_rpc=cell(length(x),length(editPhCyc1),length(editPhCyc2),length(refPhCyc));
outON_posx_epc=cell(length(x),length(editPhCyc1),length(editPhCyc2));
outOFF_posx_epc=cell(length(x),length(editPhCyc1),length(editPhCyc2));
outON_posx=cell(length(x),1);
outOFF_posx=cell(length(x),1);
outDIFF_posx=cell(length(x),1);
outON=struct([]);
outOFF=struct([]);


%loop through space: Don't forget to initialize the parallel processing
%toolbox workers using 'matlabpool open N' (for N workers, 12 max).

%for X=1:length(x);
parfor X=1:length(x);
    for EP1=1:length(editPhCyc1)
        for EP2=1:length(editPhCyc2)
            for RP=1:length(refPhCyc)
                disp(['Executing X-position ' num2str(X) ' of ' num2str(length(x)) ', '...
                    'First Edit phase cycle ' num2str(EP1) ' of ' num2str(length(editPhCyc1)) ', '...
                    'Second Edit phase cycle ' num2str(EP2) ' of ' num2str(length(editPhCyc2)) ', '...
                    'Refoc phase cycle ' num2str(RP) ' of ' num2str(length(refPhCyc)) '!!!']);
                outON_posx_epc_rpc{X}{EP1}{EP2}{RP}=sim_megaspecial_shaped(Npts,sw,Bfield,lw,taus,sys,...
                    editRFon,editTp,editPhCyc1(EP1),editPhCyc2(EP2),...
                    refRF,refTp,Gx,x(X),refPhCyc(RP));
                outOFF_posx_epc_rpc{X}{EP1}{EP2}{RP}=sim_megaspecial_shaped(Npts,sw,Bfield,lw,taus,sys,...
                    editRFoff,editTp,editPhCyc1(EP1),editPhCyc2(EP2),...
                    refRF,refTp,Gx,x(X),refPhCyc(RP));
                if RP==1
                    outON_posx_epc{X}{EP1}{EP2}=outON_posx_epc_rpc{X}{EP1}{EP2}{RP};
                    outOFF_posx_epc{X}{EP1}{EP2}=outOFF_posx_epc_rpc{X}{EP1}{EP2}{RP};
                else
                    outON_posx_epc{X}{EP1}{EP2}=op_subtractScans(outON_posx_epc{X}{EP1}{EP2},outON_posx_epc_rpc{X}{EP1}{EP2}{RP});
                    outOFF_posx_epc{X}{EP1}{EP2}=op_subtractScans(outOFF_posx_epc{X}{EP1}{EP2},outOFF_posx_epc_rpc{X}{EP1}{EP2}{RP});
                end
            end %end of refocusing phase cycle loop.
            if EP1==1 && EP2==1
                outON_posx{X}=outON_posx_epc{X}{EP1}{EP2};
                outOFF_posx{X}=outOFF_posx_epc{X}{EP1}{EP2};
            else
                outON_posx{X}=op_addScans(outON_posx{X},outON_posx_epc{X}{EP1}{EP2});
                outOFF_posx{X}=op_addScans(outOFF_posx{X},outOFF_posx_epc{X}{EP1}{EP2});
            end
            outDIFF_posx{X}=op_subtractScans(outON_posx{X},outOFF_posx{X});
        end %end of 1st editing phase cycle loop.
    end %end of 2nd editing phase cycle loop.
    outON=op_addScans(outON,outON_posx{X});
    outOFF=op_addScans(outOFF,outOFF_posx{X});
end %end of spatial loop (parfor) in x direction.
   
outDIFF=op_subtractScans(outON,outOFF);

figure 
hold
for n=1:length(x)
plot(outON_posx{n}{1}.ppm,outON_posx{n}{1}.specs+5*n);
plot(outOFF_posx{n}{1}.ppm,outOFF_posx{n}{1}.specs+5*n);
end



       