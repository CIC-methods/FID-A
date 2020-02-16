%rf_goia.m
%Jamie Near, McGill University 2020.  Based on old code from 2007.
%
% USAGE:
% [RF,FM,mv,sc]=rf_goia(N,Tp,dx,tbw,gmax,xshift)
% 
% DESCRIPTION:
% this funciton creates a GOIA gradient modulated adiabatic pulse, using 
% the method of gradient modulated offset-independent adiabaticity as first
% described by Tannus and Garwood in NMR Biomed 1997; 10:423-434. 
% 
% INPUTS:
% N              = Number of points in RF waveform (must be an even number).
% Tp             = Duration of the RF pulse in [ms].
% dx             = Desired slice thickness in [cm]
% tbw            = Time bandwidth product.
% gmax           = Maximum allowed gradient strength [G/cm] (Optional.
%                  Default = 4 G/cm == 40 mT/m)
% xshift         = Desired shift of the selected slice from isocentre in [cm].
%                  (Optional. Default=0);
%

% OUTPUTS:
% RF             = Output rf waveform for GOIA pulse, in FID-A rf pulse 
%                  structure format.
% FM             = Frequency modulation waveform (in Hz).
% mv             = Simulated magnetization vector in three columns (x,y,z) 
%                  as a function of frequency.
% sc             = Frequency scale (in kHz) corresponding to the simulated 
%                  mv vectors.


function [RF,FM,mv,sc]=rf_goia(N,Tp,dx,tbw,gmax,xshift)

if nargin<7
    xshift=0;
    if nargin<6
        gmax=4;
    end
end

%convert Tp from ms to s;
Tp=Tp/1000;

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
B=asech(0.01);

%Find Bandwith Factor A
bw=tbw/Tp;
A=bw/2;

%find time step size
dt=t(2)-t(1);

%define Gyromagnetic Ratio (Hz/G)
gyro=4257.7;



%First define the AM function:
F1=sech(B*(tau.^4));

%now calculate the FM function based on the assumption of a constant
%gradient by integrating the AM function:
F2a=zeros(1,length(F1));
F2a(length(F1)/2+1:length(F1))=cumsum(F1(length(F1)/2+1:length(F1)));
F2a(1:length(F1)/2)=-F2a(end:-1:length(F1)/2+1);
F2a=F2a/max(F2a);
%plot(tau,F2a);

%calculate a static gradient value based on desired slice thickness and
%bandwidth.  Ans in Gauss/cm.
G=bw/(gyro*dx);

%calculate an x of t function based on a constant gradient using the
%desired slice thickness.
xoft=A*F2a/(gyro*G);

%now offset the x of t function by the amount given by xoffset:
xoft=xoft+xshift;
%figure
%plot(tau,xoft);


%now define the GM function:
F3=(1-0.90*sech(B*(tau.^2)));


%now define the FM function that was derived from the solution to the
%differential equation in Tannus et al, NMR Biomed, 10:423-434 (1997).  The
%differential equation was solved by Jamie Near using Maple.
G=cosh(B*(tau2.^2))./((5*cosh((B*(tau2.^4))+(B*(tau2.^2))))+(5*cosh((B*(tau2.^4))-(B*(tau2.^2))))-(9*cosh((B*(tau2.^4)))));
h=(10*cosh(B*(tau2.^2))-9)./cosh(B*(tau2.^2));

%the above funciton G must be integrated from 0 to t and divided by the 
%expression h to obtain the value of the FM function at time t for t>0. 
%For t<0, simply mirror the function about the origin.

%integrate first half of G and multiply by h.
halfg=cumsum(G);
halfg=halfg*(tau2(2)-tau2(1)).*h;

fullg=[-halfg(end:-1:1) halfg];

F2=fullg/max(fullg);


%now calculate max gradient strength based on slice thickness desired: (assume
%slice located at x=0.
Ag=A.*(1/gyro).*(1/(dx/2));

%GM is the maximum gradient value, and it cannot exceed 4G/cm.  Therefore,
%check if AG is greated than 4, and if so, ask permission to reduce the A
%value of the pulse.  
if Ag >= gmax
    chg = input('Gradient too high for such a narrow slice.  Allow reduction of tbw? ([y] or n):  ','s');
    if isempty(chg)
        chg='y';
    end
    if chg == 'y' || chg == 'Y'
        A = gmax*gyro*(dx/2);
        disp(['Reducing max gradient to' num2str(gmax) ' G/cm.']);
        tbw=2*A*Tp/1000;
    end
end

GM=Ag*F3;
AM=F1;
FM=(A*F2);


%create new FM function based on x of t function, and gradient function:
FM2=xoft.*gyro.*GM;

%If xshift is not zero, then we need to adjust the FM function.  
FM=FM +(GM*gyro*xshift);

%now create phase modulation function using FM
ph=cumsum(FM)*dt*360;

rf(:,1)=ph;
rf(:,2)=AM;
rf(:,3)=1;
rf(:,4)=GM;

%Now find the b1max required to get full inversion:
%Since the pulse is phase modulated, so we will need to run some test to find
%out the w1max;  To do this, we can plot Mz as a function of w1 and
%find the value of w1 that results in the desired flip angle.
Tp_sim=0.005;
[mv,sc]=bes(rf,Tp_sim*1000,'b',0,0,5,40000);
plot(sc,mv(3,:));
xlabel('w1 (kHz)');
ylabel('mz');
w1max=input('Input desired w1max in kHz:  ');
w1max=w1max*1000; %convert w1max to [Hz]
tw1=Tp_sim*w1max;

RF.waveform=rf;
RF.type='inv';
RF.f0=xshift;
RF.tw1=tw1;
RF.tbw='N/A - gradient modulated pulse';
RF.isGM=true;
RF.tthk=Tp*dx;

[sc,mv]=bes(RF.waveform,Tp*1000,'f',tw1/Tp/1000,xshift-(2*dx),xshift+(2*dx),40000);

