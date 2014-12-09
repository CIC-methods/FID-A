% rf_blochSim.m
% Jamie Near, McGill University 2014.
%
% USAGE:
% [mv,sc]=rf_blochSim(RF,tp,peakB1,sc);
% 
% DESCRIPTION:
% Perform a bloch simulation of an RF pulse.  
% 
% INPUTS:
% RF        = RF pulse definition structure
% tp        = pulse duration in [ms]
% sc        = Frequency span in [kHz] (optional)
% peakB1 	= Peak B1 amplitude in [kHz] (optional)


function [mv,sc]=rf_blochSim(RF,tp,fspan,peakB1,f0);

if nargin<5
    f0=0;
    if nargin<4
        peakB1=RF.tw1/tp;
        if nargin<3
            fspan=10;
        end
    end
end        

[mv,sc]=bes(RF.waveform,tp,'f',RF.tw1/tp,f0-fspan/2,f0+fspan/2,10000);

figure
subplot(3,1,1),plot(sc,mv(1,:));
subplot(3,1,2),plot(sc,mv(2,:));
subplot(3,1,3),plot(sc,mv(3,:));
