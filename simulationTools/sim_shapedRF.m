%sim_shapedRF.m
%Jamie Near, 2014.
%
% USAGE:
% d_out = sim_shapedRF(d_in,H,RFstruct,flipAngle,phase,grad,pos)
% 
% DESCRIPTION:
% This function simulates the effect of a shaped rf pulse on the density
% matrix.  The temporal shape of the refocussing pulses is modelled as a 
% series of N instantaneous rotations about the effective RF field, where N 
% is the number of time points in the RF waveform.  The instantaneous 
% effective RF field can be an arbitrary vector, and can be represented in 
% polar coordinates as B(Beff,alpha,zeta), where Beff is the magnitude field, 
% alpha is the polar angle (the angle between the transverse plane and the 
% effective B field), and zeta is the azimuthal angle, which is given by the 
% phase of the RF).  Rotation about the effective B-field is achieved by a 
% composite rotation: Rotate about Y by -alpha, rotate about Z by -zeta, then 
% rotate about X by 2*pi*gamma*Beff*dt, then rotate back about Z by zeta and 
% back about Y by alpha.  
% 
% INPUTS:
% d_in      = input density matrix structure.
% H         = Hamiltonian operator structure.
% RF        = Radiofrequency pulse.  This can be the filename of a Siemens
%             .pta file, or an RF pulse definition structure (obtained using rf_init.m)
% Tp        = Pulse duration in [ms];
% flipAngle = RF pulse flip angle [degrees].
% phase     = Phase of RF pulse [degrees].  Optional.  Default = 0 (x'-axis.  90 degress corresponds to +y' axis)
% grad      = Gradient strength [G/cm]. Optional (for slice selective pulses only).
% pos       = Position of spins relative to voxel centre [cm].  Optional (for slice selective pulses only).


function d_out = sim_shapedRF(d_in,H,RF,Tp,flipAngle,phase,grad,pos)

gamma=42577000;

%Determine whether the input argument RF is the filename of a .pta or .RF 
%file to read in or a RF pulse definition structure, and act accordingly.
if isstr(RF)
    if flipAngle<=110
        type='exc'
    else
        type='inv'
    end
    if RF(end-3:end)=='.pta'
        if exist(RF)
            RF_struct=io_loadRFwaveform(RF,type);
        else
            error('ERROR:  RF pulse file not found! ');
        end
    elseif RF(end-2:end)=='.RF'
        if exist(RF)
            RF_struct=io_loadRFwaveform(RF,type);
        else
            error('ERROR: RF pulse file not found! ');
        end
    end
    w1max=(RF_struct.tw1*flipAngle)/(Tp*90);
elseif isstruct(RF)
    RF_struct=RF;
    w1max=(RF_struct.tw1*flipAngle)/(Tp*180);
end
   
%Determine whether this is a frequency selective pulse or a slice selective
%pulse (depends on whether the gradient and position values were given as input arguments):
if nargin<7
    %disp('No gradient or position provided.  Pulse is Frequency selective.  Using grad=0 and pos=0!');
    grad=0;
    pos=0;
elseif nargin==7
    %disp('No position was provided.  Using pos=0. ');
    pos=0;
else
    %disp('Both grad and pos were provided.  Pulse is slice selective.  ');
end

%Convert gradient value into [T/m], lengths to [m], etc...
grad=grad*0.01;  %convert from G/cm to [T/m]
pos=pos/100;  %convert from cm to [m]

%Define the properties of the refocusing radiofrequency pulse:
Tp=Tp/1000; %convert pulse duration from [ms] to seconds [s]
B1max=w1max*1000/gamma; %convert w1max from [kHz] into Tesla [T]
rfph=RF_struct.waveform(:,1)*pi/180; %convert rf phase from [degrees] into [radians]
rfB1=(RF_struct.waveform(:,2)/max(RF_struct.waveform(:,2)))*B1max; %rf amplitude in units of [T];
dt=RF_struct.waveform(:,3)*Tp/sum(RF_struct.waveform(:,3)); %rf step durations
%plot([1:length(RF_struct.waveform(:,1))],rfB1/B1max);
%pause;

%now create the RF pulse:
Rz=zeros(2^H.nspins,2^H.nspins,length(rfB1));
Ry=zeros(2^H.nspins,2^H.nspins,length(rfB1));
Rx=zeros(2^H.nspins,2^H.nspins,length(rfB1));
Beff=zeros(length(rfB1),H.nspins);  %Magnitude of the effective B1 vector [Tesla]
alpha=zeros(length(rfB1),H.nspins); %Angle between Beff and the transverse plane. [Radians]
zeta=repmat(rfph+(phase*pi/180),1,H.nspins);  %Instantaneous phase of the pulse RF [Radians]
theta=zeros(length(rfB1),H.nspins); %Angle to rotate about Beff1 [Radians]
%Now calculate the composite rotation matrices for the pulse:
for y=1:H.nspins
    %For the refocusing pulses
    Beff(:,y)=sqrt((rfB1.^2)+(((grad*pos)+(H.Bfield*H.shifts(y)/1e6))^2));  %Calculate BeffRef for first refoc pulse
    alpha(:,y)=atan2((grad*pos)+(H.Bfield*H.shifts(y)/1e6),rfB1);  %Calculate alphaRef for first refoc pulse
    theta(:,y)=2*pi*gamma*Beff(:,y).*dt;  %Calculate theta for first refoc pulse
    for n=1:length(rfB1)
        %Populate X- Y- and Z- rotation matrices for each point in RF
        %waveform.
        Rz(:,:,n)=Rz(:,:,n)+zeta(n,y)*H.Iz(:,:,y);
        Ry(:,:,n)=Ry(:,:,n)+alpha(n,y)*H.Iy(:,:,y);
        Rx(:,:,n)=Rx(:,:,n)+theta(n,y)*H.Ix(:,:,y);
    end
end

d=d_in;
%Now do the density matrix evolutions for each RF pulse element:
for n=1:length(rfB1)
    d = expm(-1i*Rz(:,:,n))*...
        expm(-1i*Ry(:,:,n))*...
        expm(-1i*Rx(:,:,n))*...
        expm(-1i*-Ry(:,:,n))*...
        expm(-1i*-Rz(:,:,n))*...    
        d*...                       %start with Mz
        expm(1i*-Rz(:,:,n))*...     %rotate about z by -zeta
        expm(1i*-Ry(:,:,n))*...     %rotate about y by -alpha
        expm(1i*Rx(:,:,n))*...      %rotate about x by theta
        expm(1i*Ry(:,:,n))*...      %rotate about y by alpha
        expm(1i*Rz(:,:,n));         %rotate about z by zeta
end

d_out=d;
    