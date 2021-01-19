%rf_gauss.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% [rf,AMPINT]=rf_gauss(tp,bw,n,type,df);
% 
% DESCRIPTION:
% Create an n point gaussian rf waveform.  Waveform can be converted in to siemens
% pta file using rf_writepta.m or into varian/agilent rf file using
% io_writeRF.m.
% 
% INPUTS:
% tp         = duration of rf pulse in [ms].
% bw         = FWHM of the gaussian inversion profile in the frequency
%              domain [Hz].
% n          = number of points in the rf waveform.
% type       = Type of pulse ('exc','ref' or 'inv').  If 'exc' is chosen,
%              the flip angle will be 90 degrees.  If 'ref' or 'inv' is 
%              chosen, the flip angle will be 180 degrees.  
% df         = frequency of gaussian pulse in [Hz].  (Optional. Default = 0 
%              Hz, which corresponds to the reference frequency of the rf 
%              transmitter.)  

%
% OUTPUTS:
% rf         = Output rf waveform for gaussian rf pulse, in FID-A rf 
%              pulse structure format.
% AMPINT     = Calculated amplitude integral (for use in Siemens .pta files).

function [rf,AMPINT]=rf_gauss(tp,bw,n,type,df);

%Make default df value = 0;
if nargin<5
    df=0;
end

%creates an n-point single banded gaussian RF pulse with duration tp(ms).
%The band will be at df Hz.  Bw is the bandwidth of the selection band in Hz.

tps=tp/1000;
asym=input('would you like to make the pulse asymmetric?  ','s');
asym_factor=0;
if asym=='y'
    asym_factor=input('Enter asymmetry factor   ');
end


t=[-tps/2+(tps/(2*n))-(tps*asym_factor):tps/n:tps/2-(tps/(2*n))-(tps*asym_factor)]';

fwhmf=bw;
fwhmt=1/fwhmf;
c=fwhmt/(2*sqrt(2*log(2)));

gauss1= exp((-t.^2)./(2*c^2));
gauss = exp((-t.^2)./(2*c^2)).*exp(-i*df*2*pi*t);  %shift second gauss by +df

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

    gauss=gauss.*fcos; 
end


%Here we compute the Amplitude integral(AMPINT), which is used by magnetom
%to calculate the transmitter power that is required in order to achieve
%the desired flip angle.  Calculation of AMPINT for dual banded pulses is
%described in the document entitled
%RF_amplitude_integral_in_pulse_sequences.pdf, which can be found in
%/Users/jnear/Documents/RF pulses/ (On jnear's local machine).

gscaled=gauss/max(gauss1);
AI=sum(gscaled);

gauss_scaled=gauss/max(abs(gauss));
AMPINT=AI/max(abs(gauss));

AMPINT=sum(abs(gauss)/max(abs(gauss)));

rfwaveform=zeros(n,3);
rfwaveform(:,1)=phase(gauss_scaled)*180/pi;
rfwaveform(:,2)=abs(gauss_scaled);
rfwaveform(:,3)=ones(n,1);

%The pulse is not phase modulated, so we can calculate the w1max:
%find the B1 max of the pulse in [kHz]:

%First find out the flip angle:
if strcmp(type,'exc')
    flipCyc=0.25; %90 degrees is 0.25 cycles;
elseif strcmp(type,'ref') || strcmp(type,'inv')
    flipCyc=0.5;  %180 degrees is 0.5 cyles;
else
    error('ERROR: type not recognized.  Use ''exc'', ''ref'' or ''inv''.  ABORTING!!');
end

intRF=sum(rfwaveform(:,2).*((-2*(rfwaveform(:,1)>179))+1))/length(rfwaveform(:,2));
if intRF~=0
    w1max=flipCyc/(intRF*tp/1000); %w1max is in [Hz]
else
    w1max=0;
end
tw1=tp*w1max/1000;

%Calculate time-bandwidth product:
tbw=tp*bw/1000;

rf.waveform=rfwaveform;
rf.type=type;
rf.tw1=tw1;
rf.tbw=tbw;
rf.f0=df;
rf.isGM=false;


