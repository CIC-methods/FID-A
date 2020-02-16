% rf_refocusedComponent.m
% Jamie Near, McGill University 2016.
%
% USAGE:
% [I,ph]=rf_refocusedComponent(RF,tp,flipAngle,fspan);
% 
% DESCRIPTION:
% Calculates the reefocused component and refocused phase of an rf pulse.
% This is done according to the description in Section 5.7.3 of "In vivo 
% NMR Spectroscopy - Principles and Techniques" by Robin A de Graaf.
% Specifically, the RF pulse is simulated twice:  once with the starting
% magnetization along Mx, and once with the startin magnetization along My.
% From these, fxx, fxy, fyy and fyx are calculated and fed into the
% equations for refocused component magnitude (I) and phase (ph).  For a
% plane rotation pulse, the refocused component should be 1 across the
% slice profile, and the phase should be uniform (or linearly varying)
% across the slice profile.  
% 
% INPUTS:
% RF        = RF pulse definition structure
% tp        = pulse duration in [ms] (optional.  Default = 5ms).
% flipAngle = flip angle of pulse [degrees] (optional.  Default = 180 deg).
% fspan     = frequency span in [kHz] or in the case of an rf pulse that 
%             includes a gradient (as indicated by the presence of 4 
%             columns in the RF waveform), the position span in [cm].  
%             (optional.  Default = 10 kHz or 10 cm).
%
% OUTPUTS:
% I         = Refocused component magnitude.
% ph        = Refocused component phase.


function [I,ph]=rf_refocusedComponent(RF,tp,flipAngle,fspan);

if nargin<4
    fspan=10;
    if nargin<3
        flipAngle=180;
        if nargin<2
            tp=5.0;
        end
    end
end

simBW=fspan; %[kHz]  Default simulation bandwidth is 10 kHz;
f0=0; %[kHz]  Centre Frequency of simulation in kHz;
b1max=RF.tw1/tp*(flipAngle/180); %[kHz]  B1max is calculated automatically;
rfph=0;  %[degrees]  Default rf pulse phase is 0 degrees;

%Simulate twice, once starting with [Mx] and once starting with [My];
[mv1,sc]=rf_blochSim(RF,tp,simBW,f0,b1max,rfph,[1 0 0]);
[mv2,sc]=rf_blochSim(RF,tp,simBW,f0,b1max,rfph,[0 1 0]);

%Close open Figures;
close;close;

%Extract fxx, fyy, fxy, and fyx
fxx=mv1(1,:);
fyy=mv2(2,:);
fxy=mv1(2,:);
fyx=mv2(1,:);

%Now calculate refocused component magnitude and phase:
I=0.5*sqrt(((fxx-fyy).^2)+((fxy+fyx).^2));
ph=0.5*atand((fxy+fyx)./(fyy-fxx)); %results in degrees

figure;
set(gcf,'Color','w');
subplot(2,1,1);
plot(sc,I,'LineWidth',1.2);
ylabel('Magnitude','FontSize',20);
title('Refocused Component','FontSize',20);
set(gca,'FontSize',16);
box off;

subplot(2,1,2);
plot(sc,ph,'LineWidth',1.2);
if size(RF,2)<4
    xlabel('Frequency [kHz]','FontSize',20);
else
    xlabel('Position [cm]','FontSize',20);
end
ylabel('Phase [deg]','FontSize',20);
set(gca,'FontSize',16);
box off;
