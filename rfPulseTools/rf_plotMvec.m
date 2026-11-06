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
% mode      = 'f' - Frequency mode (sc is a frequency scale in [kHz])
%             'pos' - Position mode (sc is a position scale in [cm] 
%             (Optional.  Default is 'f'. NOTE: the string 'x' will also 
%             activate position mode).
%
% OUTPUTS:
% N/A      

function []=rf_plotMvec(mv,sc,mode);

if nargin<3
    mode = 'f';
else
    if strcmpi(mode,'f') && strcmpi(mode,'pos') && strcmpi(mode,'x')
        error('ERROR:  Mode not recognized.');
    end
end

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

if strcmpi(mode,'f')
    xlabel('Frequency (kHz)');
elseif strcmpi(mode,'pos') || strcmpi(mode,'x')
    xlabel('Position (cm)');
end

set(gcf,'color','w');
