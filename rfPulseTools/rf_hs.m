%rf_goia.m
%Jamie Near, McGill University 2020.  Based on old code from 2007.
%
% USAGE:
% [RF,FM,mv,sc]=rf_goia(N,n,tbw,Tp,trunc,thk)
% 
% DESCRIPTION:
% this funciton creates a GOIA gradient modulated adiabatic pulse, using 
% the method of gradient modulated offset-independent adiabaticity as first
% described by Tannus and Garwood in NMR Biomed 1997; 10:423-434.  N is the 
% number of steps, n
% is the order of the HS pulse, tbw is the time bandwidth product of the
% pulse, Tp is the duration of the pulse and thk is the desired thickness of
% the pulse.
% 
% INPUTS:
% N              = Number of points in RF waveform.
% n              = order of the HS pulse.
% tbw            = Time bandwidth product.
% Tp             = Duration of the RF pulse (ms).
% trunc          = Truncation of the amplitude modulation function.
% thk            = thickness of the slice selective pulse in [cm] (optional). 
%
% OUTPUTS:
% rf             = Output rf waveform for a HS pulse, in FID-A rf pulse 
%                  structure format.
% FM             = Frequency modulation waveform (in Hz).
% mv             = Simulated magnetization vector in three columns (x,y,z) 
%                  as a function of frequency.
% sc             = Frequency scale (in kHz) corresponding to the simulated 
%                  mv vectors.


function [RF,FM,mv,sc]=rf_hs(N,n,tbw,Tp,trunc,thk)


%make sure N is even
if mod(N,2)~=0
    N=N+1;
end

%initialize the time vectors
%ta has N steps from 0 to Tp.
t=[0:Tp/(N-1):Tp];

%tau has N steps from -1 to 1. (useful for defining our AM and GM
%functions.)
tau=t*2/Tp-1;

%create time vector from 0 to 1 that is N/2 in length (useful for creation
%of FM).
tau2=tau(end/2+1:end);

%create truncation factor.
B=asech(trunc);

%Find Bandwith Factor A
bw=tbw/Tp;
A=bw/2;

%find time step size
dt=t(2)-t(1);

%define Gyromagnetic Ratio (Hz/G)
gyro=4257.7;

%First define the AM function:
F1=sech(B*(tau.^n));

%now calculate the FM function based on the assumption of a constant
%gradient by integrating the AM function:
F2=zeros(1,length(F1));
F2(length(F1)/2+1:length(F1))=cumsum(F1(length(F1)/2+1:length(F1)));
F2(1:length(F1)/2)=-F2(end:-1:length(F1)/2+1);
F2=F2/max(F2);
% figure
% subplot(1,2,1)
% plot(tau,F1);
% subplot(1,2,2)
% plot(tau,F2);

%calculate a static gradient value based on desired slice thickness in [cm] 
%and bandwidth in [Hz].  Ans in Gauss/cm.
if nargin > 6
   G=bw/(gyro*thk);
else
    G=0;
end

GM=ones(N,1)*G;
AM=F1;
FM=A*F2;
%+(G*gyro*1);



%create phase modulation function
ph=cumsum(FM)*dt*360;

waveform(:,1)=ph;
waveform(:,2)=AM;
waveform(:,3)=1;
isGM=false;
if nargin > 6
   waveform(:,4)=GM;
   isGM=true;
end


%Now find the b1max required to get full inversion:
%Since the pulse is phase modulated, so we will need to run some test to find
%out the w1max;  To do this, we can plot Mz as a function of w1 and
%find the value of w1 that results in the desired flip angle.
Tp=0.005;
[mv,sc]=bes(waveform,Tp*1000,'b',0,0,5,40000);
plot(sc,mv(3,:));
xlabel('w1 (kHz)');
ylabel('mz');
w1max=input('Input desired w1max in kHz:  ');
w1max=w1max*1000; %convert w1max to [Hz]
tw1=Tp*w1max;

RF.waveform=waveform;
RF.type='inv';
RF.f0=0;
RF.tw1=tw1;
RF.tbw=tbw;
RF.isGM=isGM;
RF.tthk=Tp*thk;  

[sc,mv]=bes(RF.waveform,Tp*1000,'f',tw1/Tp/1000,-tbw/Tp,tbw/Tp,40000);


