
% run_simSemiLASERShaped_fast.m
% Dana Goerzen and Jamie Near, McGill University 2021.
% Fast version by Muhammad G Saleh (Johns Hopkins University School of Medicine, 2019)
% USAGE:
% out=run_simPressShaped_fast(spinSys);
% 
% DESCRIPTION:
% This script simulates a PRESS experiment with fully shaped refocusing 
% pulses. Coherence order filtering is employed to only simulate desired signals
% This results in a 4x speed up compared to phase cycling (see deprecated run_simSemiLASERShaped_fast_phCyc.m)
% Furthermore, simulations are run at various locations in space to account for the 
% within-voxel spatial variation of the metabolite signal.  Summation 
% across spatial positions is performed. The MATLAB parallel computing toolbox 
% (parfor loop) was used to accelerate the simulations.  Acceleration 
% is currently performed in the direction of the slice selective pulse along
% the x-direction, but this can be changed.  Up to a factor of 12 acceleration
% can be achieved using this approach. To achieve 
% faster perfomance compared to the original 'run_simSemiLASER_shaped.m' function,
% this code uses the method described by Yan Zhang et al. Med Phys 2017;44(8): 
% 4169-78.  Some additional acceleration is currently performed using parfor 
% loops in both x and y directions.  To enable the use of the MATLAB
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
%
% OUTPUTS:
% out               = Simulation results, summed over all space.


function out=run_simSemiLASERShaped_fast(sys)

% INPUTS:
n=4096; %= number of points in fid/spectrum
sw=5000; %= desired spectral width in [Hz]
Bfield=7; %= main magnetic field strength in [T]
lw=2; %= linewidth in [Hz]
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
te=105;         %semiLASER total echo time [ms]
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
out_posx_rpc =cell(length(x),1);
d=cell(1,1); %Initialize a cell for the dentity matrix, with elements for each phase cycle;

%loop through space: Don't forget to initialize the parallel processing
%toolbox workers using 'matlabpool open N' (for N workers, 12 max).

%if you do not have the parallel computing toolbox, uncomment the first
%for loop and delete "parfor X=1:length(x)"
parfor X=1:length(x)
%  for X=1:length(x)
        disp(['Executing X-position ' num2str(X) '!!' ]);
        out_posx_rpc{X}=sim_sLASER_shaped_Ref1(Bfield,sys,te,rfPulse,refTp,x(X),Gx,flipAngle,centreFreq);
%                            sim_sLASER_shaped_Ref1(Bfield,sys,te,RF,       tp,  dx, Gx,ph1,   ph2,  centreFreq)
end

%calculate the average density matrix (Doing this inside a separate for
%loop because I couldn't figure out how to do this inside the parfor loop):
for X=1:length(x)
        d{1}=sim_dAdd(d{1},out_posx_rpc{X});
end

% %Initialize structures:
out_posy_rpc =cell(length(x),1);
out=struct([]);

%Now loop through y direction (second refoc pulse only);
parfor Y=1:length(y) %Use this if you do have the MATLAB parallel processing toolbox
%for Y=1:length(y) %Use this if you don't have the MATLAB parallel processing toolbox
        disp(['Executing Y-position ' num2str(Y) '!!' ]);
        out_posy_rpc{Y}=sim_sLASER_shaped_Ref2(d{1},n,sw,Bfield,lw,sys,te,rfPulse,refTp,y(Y),Gy,flipAngle,centreFreq);
%                            sim_sLASER_shaped_Ref2(d,   n,sw,Bfield,linewidth,sys,te,RF,       tp, dy,  Gy,ph3,    ph4,  centreFreq)
end

%Now combine the outputs;  Again, doing this inside a separate for loop
%becuase I can't figure out how to do this inside the parfor loop:
for Y=1:length(y)
        out=op_addScans(out,out_posy_rpc{Y});
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


end




%%%%%NESTED FUNCTIONS BELOW%%%%%%%%%

%% Simulate in X-direction only
function d = sim_sLASER_shaped_Ref1(Bfield,sys,te,RF,tp,dx,Gx,flipAngle,centreFreq)

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
d=sim_COF(H,d,-1);
d=sim_evolve(d,H,tau1/1000);                            %Evolve by tau1
d=sim_shapedRF(d,H,RF,tp,flipAngle,0,dx,Gx);          %1st shaped 180 degree adiabatic refocusing pulse along X gradient
d=sim_COF(H,d,1);
d=sim_evolve(d,H,tau2/1000);                            %Evolve by tau2
d=sim_shapedRF(d,H,RF,tp,flipAngle,0,dx,Gx);          %2nd shaped 180 degree adiabatic refocusing pulse along X gradient
d=sim_COF(H,d,-1);
d=sim_evolve(d,H,tau2/1000);                            %Evolve by tau2

end

%% Simulate in Y-direction only
function out = sim_sLASER_shaped_Ref2(d,n,sw,Bfield,linewidth,sys,te,RF,tp,dy,Gy,flipAngle,centreFreq)

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
d=sim_shapedRF(d,H,RF,tp,flipAngle,0,dy,Gy);          %3rd shaped 180 degree adiabatic refocusing pulse along Y gradient
d=sim_COF(H,d,1);
d=sim_evolve(d,H,tau2/1000);                            %Evolve by tau2
d=sim_shapedRF(d,H,RF,tp,flipAngle,0,dy,Gy);          %4th shaped 180 degree adiabatic refocusing pulse along Y gradient
d=sim_COF(H,d,-1);
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
out.flags.isISIS=0;
end
