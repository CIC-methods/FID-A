% io_loadRFwaveform.m
% Jamie Near, McGill University 2014.
%
% USAGE:
% [RF_struct]=io_loadRFwaveform(filename,type,f0);
% 
% DESCRIPTION:
% Initialize an RF pulse structure to contain an 
% RF Pulse waveform as well as its accompanying header
% information.  This function finds the time-bandwidth
% product (tbw) and the time-w1 product (tw1) of the pulse,
% and stores this information in the header fields of the 
% output RF structure.
% 
% INPUTS:
% filename  = filename of RF pulse waveform text file.  Can be in Siemens
%             format (.pta) or Varian/Agilent format (.RF).
% type      = Excitation ('exc'), Refocusing ('ref') or Inversion ('inv')
% f0        = centre frequency of the rf pulse [Hz].  Optional. Default=0.

function RF_struct=io_loadRFwaveform(filename,type,f0);

if nargin<3
    f0=0;
end


%Now read in the waveform:
if exist(filename)
    if filename(end-3:end)=='.pta'
        disp('Siemens format .pta RF pulse file detected!! Loading waveform now.');
        rf=io_readpta(filename);
    elseif filename(end-2:end)=='.RF'
        disp('Varian/Agilent format .RF RF pulse file detected!! Loading waveform now.');
        rf=io_readRF(filename);
    else
        error('ERROR:  RF Pulse file not recognized.  Aborting!');
    end
else
    error('ERROR:  File not found!  Aborting!');
end


%first store the waveform;
RF_struct.waveform=rf;
RF_struct.type=type;
RF_struct.f0=f0;

%Find out if the pulse is phase modulated.  If it is not, then we can
%determine the time-w1 product of the pulse quite simply.  If it is phase
%modulated (adiabatic, etc) then the determination of the time-w1 product 
%will need to me more interactive.
a=(round(rf(:,1))==180)|(round(rf(:,1))==0);

if sum(a)<length(rf(:,1))
    isPhsMod=true;
else
    isPhsMod=false;
end

%If there are any phase discontinuities in the phase function that are
%equal to a 360 degree jump we can remove these.  This will make it easier
%for rf_resample to do it's job later on:
jumps=diff(RF_struct.waveform(:,1));
jumpsAbs=(abs(jumps)>355 & abs(jumps)<365);  %Assume jumps within this range are exactly = 360 degrees.
jumpIndex=find(jumpsAbs);
for n=1:length(jumpIndex)
    RF_struct.waveform(jumpIndex(n)+1:end,1)=RF_struct.waveform(jumpIndex(n)+1:end,1)-(360*(jumps(jumpIndex(n))/abs(jumps(jumpIndex(n)))));
end




Tp=0.005;  %assume a 5 ms rf pulse;
if ~isPhsMod
    %The pulse is not phase modulated, so we can calculate the w1max:
    %find the B1 max of the pulse in [kHz]:
    if type=='exc'
        flipCyc=0.25; %90 degrees is 0.25 cycles;
    elseif type=='ref'
        flipCyc=0.5;  %180 degress is 0.5 cycles;
    elseif type=='inv'
        flipCyc=0.5;  %180 degrees is 0.5 cycles;
    end
    intRF=sum(rf(:,2).*((-2*(rf(:,1)>179))+1))/length(rf(:,2));
    if intRF~=0
        w1max=flipCyc/(intRF*Tp); %w1max is in [Hz]
    else
        w1max=0;
    end
    tw1=Tp*w1max;
else
    %The pulse is phase modulated, so we will need to run some test to find
    %out the w1max;  To do this, we can plot Mz as a function of w1 and
    %find the value of w1 that results in the desired flip angle.
    [mv,sc]=bes(rf,Tp*1000,'b',f0/1000,0,5,40000);
    plot(sc,mv(3,:));
    xlabel('w1 (kHz)');
    ylabel('mz');
    w1max=input('Input desired w1max in kHz:  ');
    w1max=w1max*1000; %convert w1max to [Hz]
    tw1=Tp*w1max;
end

RF_struct.tw1=tw1;

%now it's time to find out the time-bandwidth product:
%First make a high resolution plot the pulse profile over a wide bandwidth:
[mv,sc]=bes(rf,Tp*1000,'f',w1max/1000,-5+f0/1000,5+f0/1000,100000);

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
[mv,sc]=bes(rf,Tp*1000,'f',w1max/1000,-bw+f0/1000,bw+f0/1000,100000);

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

RF_struct.tbw=bw*Tp*1000;

