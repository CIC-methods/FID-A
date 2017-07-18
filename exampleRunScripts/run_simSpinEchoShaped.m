% run_simSpinEchoShaped.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% This script is run simply by editing the input parameters and then
% clicking "Run".
% 
% DESCRIPTION:
% This script simulates a localized spin-echo experiment with a fully shaped  
% refocusing pulse.  Phase cycling of the refocusing pulse is performed.  
% Furthermore, simulations are run at various locations in space (1-D) to 
% account for the within-voxel spatial variation of the metabolite signal 
% due to J-evolution.  Summation across phase cycles and spatial positions is
% performed.  As a result of the phase cycling and spatially resolved simulations, 
% this code takes a long time to run.  Therefore, the MATLAB parallel computing
% toolbox (parfor loop) was used to accelerate the siumulations.  Accelration 
% is performed in the direction of the slice selective refocusing pulse.  
% Up to a factor of 12 acceleration can be achieved using this approach.  
% To enable the use of the MATLAB parallel computing toolbox, initialize 
% the multiple worked nodes using "matlabpool size X" where "X" is the 
% number of available processing nodes.  If the parallel processing toolbox
% is not available, then replace the "parfor" loop with a "for" loop.
% 
% INPUTS:
% To run this script, edit the parameters below as desired and then click
% "run":
% RFWaveform        = name of refocusing pulse waveform.
% Tp                = duration of refocusing pulses[ms]
% Bfield            = Magnetic field strength in [T]
% Npts              = number of spectral points
% sw                = spectral width [Hz]
% lw                = linewidth of the output spectrum [Hz]
% thk               = slice thickness of refocusing pulse [cm]
% pos               = vector of positions to simulate [cm]
% TE                = Echo-Time  [ms]
% spinSys           = spin system to simulate 
% PhCyc             = vector of phase cycling steps for refocusing pulse [degrees]
%
% OUTPUTS:
% out_pos           = Simulation results, spatially resolved.
% out               = Simulation results, summed over all space.


% ************INPUT PARAMETERS**********************************
RFWaveform='sampleRefocPulse.pta'; %name of refocusing pulse waveform.
Tp=5; %duration of refocusing pulses[ms]
Npts=2048; %number of spectral points
sw=2000; %spectral width [Hz]
lw=2; %linewidth of the output spectrum [Hz]
Bfield=3; %Magnetic field strength in [T]
thk=3; %slice thickness of x refocusing pulse [cm]
pos=linspace(-3,3,12); %X positions to simulate [cm]
TE=68; %Echo time of the pulse sequence [ms]
spinSys='Cr'; %spin system to simulate
phCyc=[0,90]; %phase cycling steps for 1st refocusing pulse [degrees]
% ************END OF INPUT PARAMETERS**********************************

%Load RF waveforms
RF=io_loadRFwaveform(RFWaveform,'ref',0);

gamma=42577000; %gyromagnetic ratio

%Load spin systems
load spinSystems
sys=eval(['sys' spinSys]);

%Resample refocusing RF pulse from 400 pts to 100 pts to reduce
%computational workload
RF=rf_resample(RF,100);

Gx=(RF.tbw/(Tp/1000))/(gamma*thk/10000); %[G/cm]

%n=1;
%totalIters=length(x)*length(y)*length(editPhCyc1)*length(editPhCyc2)*length(refPhCyc1)*length(refPhCyc2);

%Initialize structures:
out_pos_rpc=cell(length(pos),length(phCyc));
out_pos=cell(length(pos),1);
out=struct([]);

%loop through space: Don't forget to initialize the parallel processing
%toolbox workers using 'matlabpool open N' (for N workers, 12 max).

%for X=1:length(x);
parfor X=1:length(pos);
                for RP=1:length(phCyc)
                        disp(['Executing X-position ' num2str(X) ' of ' num2str(length(pos)) ', '...
                            'Refoc phase cycle ' num2str(RP) ' of ' num2str(length(phCyc)) '!!!']); 
                        out_pos_rpc{X}{RP}=sim_spinecho_shaped(Npts,sw,Bfield,lw,sys,TE,RF,Tp,Gx,pos(X),phCyc(RP));
                        if RP==1
                            out_pos{X}=out_pos_rpc{X}{RP};
                        else
                            out_pos{X}=op_subtractScans(out_pos{X},out_pos_rpc{X}{RP});
                        end
                    end %end of 1st refocusing phase cycle loop
        out=op_addScans(out,out_pos{X});
end %end of spatial loop (parfor) in x direction.
        






       