% rf_blochSim.m
% Jamie Near, McGill University 2014.
%
% USAGE:
% [mv,sc]=rf_blochSim(RF,tp,fspan,f0,peakB1,ph,M0);
% 
% DESCRIPTION:
% Perform a bloch simulation of an RF pulse.  This code simply runs Martyn 
% Klassen's excellent bloch equation simulator.  For more information, see
% help file for bes.m.  (~FID-A/rfPulseTools/mklassenTools/bes.m).
% 
% INPUTS:
% RF        = RF pulse definition structure
% tp        = pulse duration in [ms]
% fspan     = Frequency span in [kHz] (optional.  Default=10kHz)
% f0        = Centre of frequnecy span [kHz] (optional.  Default=0)
% peakB1 	= Peak B1 amplitude in [kHz] (optional.  Default=RF.tw1/tp)
% ph        = Starting phase of the rf pulse [degrees] (optional. Default=0)
% M0        = Starting magnetization [units of M0] (optional. Default=[0,0,1])


function [mv,sc]=rf_blochSim(RF,tp,fspan,f0,peakB1,ph,M0);

if nargin<7
    M0=[0,0,1];
    if nargin<6
        ph=0;
        if nargin<5
            peakB1=RF.tw1/tp;
            if nargin<4
                f0=0;
                if nargin<3
                    fspan=10;
                end
            end
        end
    end
end

[mv,sc]=bes(RF.waveform,tp,'f',peakB1,f0-fspan/2,f0+fspan/2,10000,ph,M0);

figure
subplot(3,1,1),plot(sc,mv(1,:));
subplot(3,1,2),plot(sc,mv(2,:));
subplot(3,1,3),plot(sc,mv(3,:));
