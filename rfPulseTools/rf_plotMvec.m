% rf_plotMvec.m
% Jamie Near, Sunnybrook Research Institute 2026.
%
% USAGE:
% [~]=rf_plotMvec(mv,sc);
% 
% DESCRIPTION:
% Plot Mx, My, Mxy, and Mz magnetization profiles as a function of
% frequency or position.  
% 
% INPUTS:
% mv        = Magnetization vectors (3 x N) returned by rf_blochSim, where
%             the three dimensions correspond to Mx, My and Mz, and N is 
%             the number of points along the frequency or position axis. 
% sc        = Vector of frequencies in [kHz] or positions in [cm] to plot
%             against.
%
% OUTPUTS:
% N/A      

function []=rf_plotMvec(mv,sc);

%Make a subplot figure with Mx, My, Mxy, and Mz:
%Mx
figure
subplot(4,1,1),plot(sc,mv(1,:),'LineWidth',1.5);
box off;
ylabel('M_x');

%My
subplot(4,1,2),plot(sc,mv(2,:),'LineWidth',1.5);
box off;
ylabel('M_y');

%Mxy
subplot(4,1,3),plot(sc,sqrt(mv(2,:).^2 + mv(1,:).^2),'LineWidth',1.5);
box off;
ylabel('M_x_y');

%Mz
subplot(4,1,4),plot(sc,mv(3,:),'LineWidth',1.5);
box off;
ylabel('M_z');

if size(RF.waveform,2)<4
    xlabel('Frequency (kHz)');
else
    xlabel('Position (cm)');
end

set(gcf,'color','w');
