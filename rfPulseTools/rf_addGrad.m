% rf_addGrad.m
% Jamie Near, McGill University 2020.
%
% USAGE:
% RF_out=rf_addGrad(RF_in,grad);
% 
% DESCRIPTION:
% Add a gradient waveform to the input RF pulse. 
% 
% INPUTS:
% RF_in     = Input RF pulse definition structure
% grad      = If grad is a scalar, then rf_addGrad will add a constant
%             gradient with amplitude equal to the value of grad.  If grad 
%             is a vector, then rf_addGrad will add the specified gradient 
%             vector to the RF waveform.  In either case, the gradient 
%             should be specified in units of [G/cm].
% 
% OUTPUTS:
% RF_out    = Output rf waveform following the addition of gradient 
%             waveform.

function RF_out=rf_addGrad(RF_in,grad);

if ~isstruct(RF_in)
    error('ERROR:  the input RF pulse must be in structure format.  Try using rf_readwaveform to convert it!  Aborting.  ');
end

newWaveform=RF_in.waveform;

%Check if the input RF pulse already has a gradient. 
if size(newWaveform,2)>3
    keepGoing='n';
    disp('WARNING:  Input waveform already has a gradient waveform!');
    keepGoing=input('Do you wish to overwrite the existing gradient waveform (y or [n])','s');
    
    if strcmp(keepGoing,'y') || strcmp(keepGoing,'Y');
        disp('OK.  Overwriting existing gradient waveform.');
    elseif strcmp(keepGoing,'n') || strcmp(keepGoing,'N');
        error('ABORTING!!!');
    else
        error('Response not recognized.  ABORTING!!');
    end
end

%Now add the gradient. 
if isscalar(grad)
    %Adding a scalar gradient
    newWaveform(:,4)=grad*ones(size(newWaveform,1),1);
else
    %Adding a gradient waveform
    if size(grad) ~= size(newWaveform,1)
        if size(grad') ~= size(newWaveform,1);
            error('ERROR:  Gradient waveform does not match the length of the input RF pulse waveform!!  ABORTING!!');
        end
        newWaveform(:,4)=grad';
    end
    newWaveform(:,4)=grad';
end


%now it's time to find out the time-bandwidth product:
%First make a high resolution plot the pulse profile over a wide bandwidth:
[mv,sc]=bes(rf,Tp*1000,'f',w1max/1000,-5+f0/1000,5+f0/1000,100000);
if isstr(RF_in.type)
    if RF_in.type=='exc'
        index=find(mv(3,:)<0.5);
        bw=sc(index(end))-sc(index(1));
        %plot(sc(index),mv(3,index),'.-',sc,mv(3,:));
    elseif RF_in.type=='ref'
        index=find(mv(3,:)<0);
        bw=sc(index(end))-sc(index(1));
        %plot(sc(index),mv(3,index),'.-',sc,mv(3,:));
    elseif RF_in.type=='inv'
        index=find(mv(3,:)<0);
        bw=sc(index(end))-sc(index(1));
        %plot(sc(index),mv(3,index),'.-',sc,mv(3,:));
    end
elseif isnumeric(RF_in.type)
    mz=cos(type); %Find out the Mz value immediately following the pulse.
    thr=(1+mz)/2;  %Find out the Mz value mid-way between 1 and the mz (half-max):
    index=find(mv(3,:)<thr);  %Find the indices of the corresponding "full width"
    bw=sc(index(end))-sc(index(1));  %Now find the bandwidth at that point ("Full width at half max").  
end


%Now make a very high resolution plot the pulse profile over a narrower bandwidth:
[mv,sc]=bes(rf,Tp*1000,'f',w1max/1000,-bw+f0/1000,bw+f0/1000,100000);
if isstr(RF_in.type)
    if RF_in.type=='exc'
        index=find(mv(3,:)<0.5);
        bw=sc(index(end))-sc(index(1));
        %plot(sc(index),mv(3,index),'.-',sc,mv(3,:));
    elseif RF_in.type=='ref'
        index=find(mv(3,:)<0);
        bw=sc(index(end))-sc(index(1));
        %plot(sc(index),mv(3,index),'.-',sc,mv(3,:));
    elseif RF_in.type=='inv'
        index=find(mv(3,:)<0);
        bw=sc(index(end))-sc(index(1));
        %plot(sc(index),mv(3,index),'.-',sc,mv(3,:));
    end
elseif isnumeric(RF_in.type)
    mz=cos(RF_in.type); %Find out the Mz value immediately following the pulse.
    thr=(1+mz)/2;  %Find out the Mz value mid-way between 1 and the mz (half-max):
    index=find(mv(3,:)<thr);  %Find the indices of the corresponding "full width"
    bw=sc(index(end))-sc(index(1));  %Now find the bandwidth at that point ("Full width at half max").
end

RF_out=RF_in;
RF_out.waveform=newWaveform;
RF_out.isGM=true;
RF_out.tbw='N/A - gradient modulated pulse';
RF_out.isGM=true;
RF_out.tthk=bw*Tp; %This is the time x sliceThickness product for
                   %gradient modulated pulses.  It is in units [cm.s]


