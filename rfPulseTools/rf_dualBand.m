%rf_dualBand.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% [rf,AMPINT]=rf_dualBand(tp,df,n,bw,ph,shft)
% 
% DESCRIPTION:
% creates an n-point dual banded gaussian inversion RF pulse with duration 
% tp(ms).  The first band will be at f=0Hz and the second band will be at 
% df Hz.  Bw is the bandwidth of the two selection bands in Hz.
% 
% INPUTS:
% tp         = pulse duration in ms.
% df         = frequency of 2nd gaussian band [Hz].
% n          = number of points in rf waveform.
% bw         = bandwidth of both selection bands [Hz].
% ph         = phase of the second gaussian.
% shft       = frequency shift applied to both bands.
%
% OUTPUTS:
% rf         = Output rf waveform for a dual banded rf pulse, in FID-A rf 
%              pulse structure format.
% AMPINT     = Calculated amplitude integral (for use in Siemens .pta files).

function [rf,AMPINT]=rf_dualBand(tp,df,n,bw,ph,shft)

%convert the pulse duration into seconds;
tps=tp/1000;

%Make a time-scale for the rf pulse.
t=[-tps/2+(tps/(2*n)):tps/n:tps/2-(tps/(2*n))]';

%Calculate the fwhm and c parameter for the first band
fwhmf1=bw;
fwhmt1=1/fwhmf1;
c1=fwhmt1/(2*sqrt(2*log(2)));

%calculate the fwhm and c parameter for the second band
fwhmf2=bw;
fwhmt2=1/fwhmf2;
c2=fwhmt2/(2*sqrt(2*log(2)));

%calculate the time-domain gaussian waveforms for the 1st and 2nd bands.  
gauss1=exp((-t.^2)./(2*c1^2));
gauss2=  exp((-t.^2)./(2*c2^2)).*exp(-i*df*2*pi*t);  %shift second gauss by +df

%add phase to the second band
gauss2=addphase(gauss2,ph);


%Here we compute the Amplitude integral(AMPINT), which is used by magnetom
%to calculate the transmitter power that is required in order to achieve
%the desired flip angle.

gscaled=gauss1/max(gauss1);
AI=sum(gscaled);

dgauss=(gauss1+gauss2).*exp(-i*shft*2*pi*t);

dgauss_scaled=dgauss/max(abs(dgauss));
AMPINT=AI/max(abs(dgauss));

rfwaveform=zeros(n,3);
rfwaveform(:,1)=phase(dgauss_scaled)*180/pi;
rfwaveform(:,2)=abs(dgauss_scaled);
rfwaveform(:,3)=ones(n,1);

%Now find out the w1-max of the pulse:
[mv,sc]=bes(rfwaveform,tp,'b',0,0,5,40000);
plot(sc,mv(3,:));
xlabel('w1 (kHz)');
ylabel('mz');
w1max=input('Input desired w1max in kHz:  ');
w1max=w1max*1000; %convert w1max to [Hz]
tw1=tps*w1max;

%for the time bandwidth product, we can simply calculate:
tbw=2*tps*bw;

rf.waveform=rfwaveform;
rf.type='inv';
rf.tw1=tw1;
rf.tbw=tbw;






