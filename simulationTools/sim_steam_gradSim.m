%sim_steam_gradSim.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% Script can be run by pressing "run".
% 
% DESCRIPTION:
% This function runs the sim_steam_spoil function multiple times with
% different spoiler gradient intensities.  The result is a spoiled STEAM
% sequence. 
% 
% INPUTS:
% Initialize the following variables and then click "run":
% spinsys       = Spin system.
% TE            = Echo time [ms].
% TM            = Mixing time [ms].
% N             = Number of 'phase cycles'
%
% OUTPUTS:
% steam         = simulated spectrum, in FID-A structure format, using STEAM 
%                 sequence.
% press         = simulated spectrum, in FID-A structure format, using PRESS 
%                 sequence (for comparison).

% *********INPUT VARIABLES***********
spinsys='Lac';      %Spin system.
TE=20;              %Echo time. [ms]
TM=5;               %Mixing time. [ms]
N=32;               %Number of 'phase cycles'
n=2048;             %Number of spectral points
sw=2000;            %Spectral width
B0=3;               %Magnetic Field Strength
lw=2;               %Linewidth
% ***********************************

eval(['load ' spinsys]);
eval(['J=sys' spinsys '.J;']);
eval(['shifts=sys' spinsys '.shifts;']);
sys.J=J;
sys.shifts=shifts;

steam=sim_steam(n,sw,B0,lw,sys,TE,TM,0);
figure;
hold;
plot(steam.ppm,real(steam.specs));

for spoil=360/N:360/N:(360)-(360/N)
    steam_temp=sim_steam(n,sw,B0,lw,sys,TE,TM,spoil);
    %plot(steam_temp.ppm,steam_temp.specs);
    steam=op_addScans(steam,steam_temp);
end

steam=op_ampScale(steam,1/N);

press=sim_press(n,sw,B0,lw,sys,TE/2,TE/2);

figure;
plot(steam.ppm,real(steam.specs),press.ppm,real(press.specs));
set(gca,'XDir','reverse');
xlim([1 5]);

legend('press','steam');




    
    