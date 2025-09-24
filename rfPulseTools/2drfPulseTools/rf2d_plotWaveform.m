% rf2d_plotWaveform.m
% Jamie Near, 2025.
% 
% USAGE:
% out=rf2d_plotWaveform(rf,mode);
% 
% DESCRIPTION:
% Plot the 2DRF pulse waveforms in the time domain.  
% 
% INPUTS:
% rf2d     = RF pulse in FID-A structure format.  
% mode   = Which waveform to plot:  'ph', 'amp', 'gm' or 'all' (optional.  Default = all).
%
% OUTPUTS:
% out    = Figure handle.

function out=rf2d_plotWaveform(rf2d,mode)
cond = 7;
[amp, phase, g_x, g_y, g_z] = rf2d_pulse_readin.read_pulse_for_bloch(cond);


%First make the time vector to use for all plots:
t=linspace(0, rf2d.tp * 1000, length(rf2d.waveform(:,2)));
tgrad = linspace(0, rf2d.tp * 1000, length(g_x));

%Check the mode and plot accordingly:
if strcmp(mode,'ph');
    out=figure;
    set(gcf,'color','w');
    plot(t,rf2d.waveform(:,2),'LineWidth',1.2);
    set(gca,'FontSize',16);
    xlabel('Time (ms)','FontSize',20);
    ylabel('Phase (rad)','FontSize',20);
    title('2DRF Pulse Phase');
elseif strcmp(mode,'amp');
    out=figure;
    set(gcf,'color','w');
    plot(t,rf2d.waveform(:,1),'LineWidth',1.2);
    set(gca,'FontSize',16);
    xlabel('Time (ms)','FontSize',20);
    ylabel('Amplitude (arb units)','FontSize',20);
    title('2DRF Pulse Amplitude');
elseif strcmp(mode,'gm')
    out1=figure;
    set(gcf,'color','w');
    plot(tgrad,g_x,'LineWidth',1.2,'DisplayName','Gradient in X');
    hold on;
    plot(tgrad,g_y,'LineWidth',1.2,'DisplayName','Gradient in Y');
    set(gca,'FontSize',16);
    xlabel('Time (ms)','FontSize',20);
    ylabel('Gradient Amplitude (mT/m)','FontSize',20);
    legend;
    title('2DRF Pulse Gradients');
elseif strcmp(mode,'all');
    out1=figure;
    set(gcf,'color','w');
    plot(tgrad,g_x,'LineWidth',1.2,'DisplayName','Gradient in X');
    hold on;
    plot(tgrad,g_y,'LineWidth',1.2,'DisplayName','Gradient in Y');
    set(gca,'FontSize',16);
    xlabel('Time (ms)','FontSize',20);
    ylabel('Gradient Amplitude (mT/m)','FontSize',20);
    legend;
    title('2DRF Pulse Gradients');

    hold off;

    out2=figure;
    set(gcf,'color','w');
    plot(t,rf2d.waveform(:,1),'LineWidth',1.2);
    set(gca,'FontSize',16);
    xlabel('Time (ms)','FontSize',20);
    ylabel('Amplitude (arb units)','FontSize',20);
    title('2DRF Pulse Amplitude');

    hold off;

    out3=figure;
    set(gcf,'color','w');
    plot(t,rf2d.waveform(:,2),'LineWidth',1.2);
    set(gca,'FontSize',16);
    xlabel('Time (ms)','FontSize',20);
    ylabel('Phase (rad)','FontSize',20);
    title('2DRF Pulse Phase');
else
    error('ERROR:  mode not recognized.  Choose ''ph'', ''amp'', ''gm'', or ''all''.  ABORTING!!');
end
end

        
        
    
  
        
    
    

