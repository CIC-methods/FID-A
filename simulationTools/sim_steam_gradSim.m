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
% TE            = Echo time [s].
% TM            = Mixing time [s].
% N             = Number of 'phase cycles'

% *********INPUT VARIABLES***********
spinsys='Lac';  %Spin system.
TE=0.020;  %Echo time.
TM=0.005;  %Mixing time.
N=32;      %Number of 'phase cycles'
% ***********************************

eval(['load ' spinsys]);
eval(['J=j' spinsys ';']);
eval(['shifts=shifts' spinsys ';']);
sys.J=J;
sys.shifts=shifts;

steam=sim_steam(2048,2000,3,2,sys,TE,TM,0);
figure;
hold;
plot(steam.ppm,steam.specs);

for spoil=360/N:360/N:(360)-(360/N)
    steam_temp=sim_steam(2048,2000,3,2,sys,TE,TM,spoil);
    plot(steam_temp.ppm,steam_temp.specs);
    steam=op_addScans(steam,steam_temp);
end

steam=op_ampScale(steam,1/N);

press=sim_press(2048,2000,3,2,sys,TE/2,TE/2);

figure;
plot(steam.ppm,steam.specs,press.ppm,press.specs);
set(gca,'XDir','reverse');
xlim([1 5]);

%legend('press','steam');




    
    