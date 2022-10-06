% io_loadRFwaveform.m
% Jamie Near, McGill University 2014.
%
% USAGE:
% [RF_struct]=io_loadRFwaveform(filename,type,off_res);
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
%             format (.pta), Varian/Agilent format (.RF), Bruker format 
%             (.inv, .ref or .exc) or a plain text file (.txt, with two, 
%             three or four columns (amplitude, phase, time-step
%             (optional),and gradient strength (optional).  Filename can 
%             also be the name of a two, three or four column matlab vector 
%             specifing the phase, amplitude, time vectors (optionally) and 
%             gradient (optionally) of an RF pulse waveform.  Note about 
%             Units: RF amplitude waveform [arbitrary units], B1 intensity 
%             [kHz],  phase waveform [degrees], time step vector [arbitrary 
%             units],  gradient waveform [G/cm].
% type      = Excitation ('exc'), Refocusing ('ref') or Inversion ('inv'). 
%             Alternatively, 'type' can specify the exact flip angle [in 
%             degrees] of the pulse.  This can be useful if you want to use 
%             a pulse that does not have an exact flip angle of 90 or 180 
%             degrees (Thanks to Eric Plitman for help with this feature).  
% off_res   = set to "1" if the centre of the RF pulse is not at 0 Hz. 
%             Optional. Default=0.
%
% OUTPUTS:
% RF_struct = RF pulse waveform in FID-A rf pulse structure format.

function RF_struct=io_loadRFwaveform(filename,type,off_res);

if nargin<3
    off_res=0;
end

if isnumeric(type)
    %If 'type' was specified to be exactly 90 degrees, set the type to 'exc';
    if type==90;
        type='exc';
    end
    
    %If 'type' was specified to be exactly 180 degrees, set the type to 'inv'
    %(functionally there is no real difference between 'ref' and 'inv', so I
    %just chose 'Inv';
    if type==180;
        type='inv';
    end
end

%Now read in the waveform:
if isstr(filename)
    if exist(filename)
        if filename(end-3:end)=='.pta'
            disp('Siemens format .pta RF pulse file detected!! Loading waveform now.');
            rf=io_readpta(filename);
        elseif filename(end-2:end)=='.RF'
            disp('Varian/Agilent format .RF RF pulse file detected!! Loading waveform now.');
            rf=io_readRF(filename);
        elseif filename(end-3:end)=='.inv'
            disp('Bruker format .inv RF pulse file detected!! Loading waveform now.');
            rf=io_readRFBruk(filename);
        elseif filename(end-3:end)=='.rfc'
            disp('Bruker format .pta RF pulse file detected!! Loading waveform now.');
            rf=io_readRFBruk(filename);
        elseif filename(end-3:end)=='.exc'
            disp('Bruker format .exc RF pulse file detected!! Loading waveform now.');
            rf=io_readRFBruk(filename);
        elseif filename(end-3:end)=='.txt'
            disp('Basic .txt format RF pulse file detected!! Loading waveform now.');
            rf=io_readRFtxt(filename);
        else
            error('ERROR:  RF Pulse file not recognized.  Aborting!');
        end
    else
        error('ERROR:  File not found!  Aborting!');
    end  
else
    if ismatrix(filename) && ndims(filename)==2 && ( 2 <= size(filename,2) <= 4 )
        disp('Input is an RF waveform already in the matlab workspace.  Loading waveform now.');
        if size(filename,2)==2
            rf=filename;
            rf(:,3)=ones(length(filename(:,1)),1);
        elseif size(filename,2)==3 || size(filename,2)==4
            rf=filename;
        end
    end
end


Tp=0.005;  %assume a 5 ms rf pulse;

%Frequency shifting the pulse to zero if off_res is 1. To be shifted back
%afterwards - PT,2021
if off_res   
    %Initializations
    tmp_w1=0.03;
    tmp_min=1;
    
    %Looping through until a min of <-0.98 (inv/ref) or <0.2 (exc) is reached
    switch type
        case 'exc'
            min_floor = 0.2;
        case {'inv','ref'}
            min_floor = -0.98;
    end
    while tmp_min > min_floor
        [mv,sc]=bes(rf,Tp*1000,'f',tmp_w1,-5,5,1000); %running bes with bare settings
        [tmp_min,tmp_min_pos]=min(mv(3,:)); %finding the min and min position
        freq_shift=sc(tmp_min_pos)*1000;%the freq shift (in Hz) of the pulse @ 5ms
        tmp_w1=tmp_w1+0.02;
    end
    
    %frequency shifting the pulse to 0 Hz
    N=size(rf,1);
    dt=Tp/N;
    t=[0:dt:Tp-dt];
    phaseRamp=t*-freq_shift*360;
    rf(:,1)=rf(:,1)+phaseRamp'; 
end

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
jumps=diff(rf(:,1));
jumpsAbs=(abs(jumps)>355 & abs(jumps)<365);  %Assume jumps within this range are exactly = 360 degrees.
jumpIndex=find(jumpsAbs);
for n=1:length(jumpIndex)
    rf(jumpIndex(n)+1:end,1)=rf(jumpIndex(n)+1:end,1)-(360*(jumps(jumpIndex(n))/abs(jumps(jumpIndex(n)))));
end

%scale amplitude function so that maximum value is 1:
rf(:,2)=rf(:,2)./max(rf(:,2));

if ~isPhsMod 
    %The pulse is not phase modulated, so we can calculate the w1max:
    %find the B1 max of the pulse in [kHz]:
    if isstr(type)
        if type=='exc'
            flipCyc=0.25; %90 degrees is 0.25 cycles;
        elseif type=='ref'
            flipCyc=0.5;  %180 degress is 0.5 cycles;
        elseif type=='inv'
            flipCyc=0.5;  %180 degrees is 0.5 cycles;
        end
    elseif isnumeric(type) %assume that a flip angle (in degrees) was given
        flipCyc=type/360;
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
    [mv,sc]=bes(rf,Tp*1000,'b',0,0,5,40000);
    plot(sc,mv(3,:));
    xlabel('w1 (kHz)');
    ylabel('mz');
    w1max=input('Input desired w1max in kHz (for 5.00 ms pulse):  ');
    w1max=w1max*1000; %convert w1max to [Hz]
    tw1=Tp*w1max;
end

%now it's time to find out the time-bandwidth product:
%First make a high resolution plot the pulse profile over a wide bandwidth:
[mv,sc]=bes(rf,Tp*1000,'f',w1max/1000,-5,5,100000);
if isstr(type)
    if type=='exc'
        index=find(abs(mv(1,:)+1i*mv(2,:))>0.5);
        bw=sc(index(end))-sc(index(1))
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
elseif isnumeric(type)
    mz=cos(type); %Find out the Mz value immediately following the pulse.
    thr=(1+mz)/2;  %Find out the Mz value mid-way between 1 and the mz (half-max):
    index=find(mv(3,:)<thr);  %Find the indices of the corresponding "full width"
    bw=sc(index(end))-sc(index(1));  %Now find the bandwidth at that point ("Full width at half max").  
end


%Now make a very high resolution plot the pulse profile over a narrower bandwidth:
[mv,sc]=bes(rf,Tp*1000,'f',w1max/1000,-bw,bw,100000);
if isstr(type)
    if type=='exc'
        index=find(abs(mv(1,:)+1i*mv(2,:))>0.5);
        bw=sc(index(end))-sc(index(1))
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
elseif isnumeric(type)
    mz=cos(type); %Find out the Mz value immediately following the pulse.
    thr=(1+mz)/2;  %Find out the Mz value mid-way between 1 and the mz (half-max):
    index=find(mv(3,:)<thr);  %Find the indices of the corresponding "full width"
    bw=sc(index(end))-sc(index(1));  %Now find the bandwidth at that point ("Full width at half max").
end

%If the pulse was offset from 0 Hz, revert the previous frequency shift - PT,2021
if off_res
    rf(:,1)=rf(:,1)-phaseRamp';
end

%Now store the output structure;
RF_struct.waveform=rf;
RF_struct.type=type;
RF_struct.f0=0;
if size(rf,2)>3 && any(rf(:,4))
    RF_struct.tbw='N/A - gradient modulated pulse';
    RF_struct.isGM=true;
    RF_struct.tthk=bw*Tp; %This is the time x sliceThickness product for 
                         %gradient modulated pulses.  It is in units [cm.s]
else
    RF_struct.tbw=bw*Tp*1000;
    RF_struct.isGM=false;
    RF_struct.tthk='N/A - frequency selective pulse';
end
RF_struct.tw1=tw1;


