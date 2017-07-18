%rf_sinc.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% [rf,AMPINT]=rf_sinc(lobes,n,type);
% 
% DESCRIPTION:
% Create an n point sinc rf waveform.  Waveform can be converted in to siemens
% pta file using rf_writepta.m or into varian/agilent rf file using
% io_writeRF.m.
% 
% INPUTS:
% lobes      = Number of lobes in the sinc pulse.  
% n          = number of points in the rf waveform.
% type       = Type of pulse: 
%               Refocusing = 'ref'
%               Inversion  = 'inv'
%               Excitation = 'exc'
% OUTPUTS:
% rf         = Output rf waveform for sinc shaped rf pulse, in FID-A rf 
%              pulse structure format.
% AMPINT     = Calculated amplitude integral (for use in Siemens .pta files).

function [rf,AMPINT]=rf_sinc(lobes,n,type);

%creates an n-point sinc RF pulse.
%The band will be at df Hz.  Bw will depend on number of lobes;

if nargin<3
    df=0;
end

if lobes<1
    error('ERROR:  sinc pulse must have at least 1 lobe!!  ABORTING!!');
else
    t=linspace(-(0.5+lobes/2),0.5+lobes/2,n);
end

AMfunc=sinc(t);
ph=180*(AMfunc<0);

%filter the gaussians using a cosine filter to minimize any rining
%artefacts.
filts=input('would you like to filter them to minimize ringing? (y or n):  ','s');
if filts=='y' || filts=='Y'
    
    %be sure to look at the plot of the filter.  If the scaling factor
    %before max(t) is too low, and/or the asymmetry factor is too low, the
    %filter will be negative at the tails and will invert  your rf pulse
    %(which is not good).  
    attn=input('by what factor would you like to attenuate the edges of the pulse?  ');
    attnfactor=pi/acos(attn);
    
    fcos=cos(t*pi/(attnfactor*max(t)));
    %plot(t,fcos);

    AMfunc=AMfunc.*fcos; 
end


rfwaveform=zeros(n,3);
rfwaveform(:,1)=ph;
rfwaveform(:,2)=abs(AMfunc);
rfwaveform(:,3)=ones(n,1);

%The pulse is not phase modulated, so we can calculate the w1max:
%find the B1 max of the pulse in [kHz]:

if type=='exc';
    flipCyc=0.25;  %90 degrees is 0.25 cycles
elseif strcmp(type,'ref') || strcmp(type,'inv')
    flipCyc=0.5;   %180 degrees is 0.5 cycles
end

Tp=0.005;  %assume 5ms rf pulse;
intRF=sum(rfwaveform(:,2).*((-2*(rfwaveform(:,1)>179))+1))/length(rfwaveform(:,2));
if intRF~=0
    w1max=flipCyc/(intRF*Tp); %w1max is in [Hz]
else
    w1max=0;
end
tw1=Tp*w1max;

%now it's time to find out the time-bandwidth product:
%First make a high resolution plot the pulse profile over a wide bandwidth:
[mv,sc]=bes(rfwaveform,Tp*1000,'f',w1max/1000,-5,5,100000);

if type=='exc'
    index=find(mv(3,:)<0.5);
    bw=sc(index(end))-sc(index(1));
    %plot(sc(index),mv(3,index),'.-',sc,mv(3,:));
elseif type=='ref'
    index=find(mv(3,:)<0);
    bw=sc(index(end))-sc(index(1));
    %plot(sc(index),mv(3,index),'.-',sc,mv(3,:));
elseif type=='inv'
    index=find(mv(3,:)<0);
    bw=sc(index(end))-sc(index(1));
    %plot(sc(index),mv(3,index),'.-',sc,mv(3,:));
end

%Now make a very high resolution plot the pulse profile over a narrower bandwidth:
[mv,sc]=bes(rfwaveform,Tp*1000,'f',w1max/1000,-bw,bw,100000);

if type=='exc'
    index=find(mv(3,:)<0.5);
    bw=sc(index(end))-sc(index(1));
    %plot(sc(index),mv(3,index),'.-',sc,mv(3,:));
elseif type=='ref'
    index=find(mv(3,:)<0);
    bw=sc(index(end))-sc(index(1));
    %plot(sc(index),mv(3,index),'.-',sc,mv(3,:));
elseif type=='inv'
    index=find(mv(3,:)<0);
    bw=sc(index(end))-sc(index(1));
    %plot(sc(index),mv(3,index),'.-',sc,mv(3,:));
end

tbw=bw*Tp*1000;

rf.waveform=rfwaveform;
rf.type=type;
rf.tw1=tw1;
rf.tbw=tbw;
rf.f0=0;


