% rf_plotWaveform.m
% Jamie Near, McGill University 2020.
% 
% USAGE:
% out=rf_plotWaveform(rf,mode);
% 
% DESCRIPTION:
% Plot the RF pulse waveforms in the time domain.  
% 
% INPUTS:
% rf     = RF pulse in FID-A structure format.  
% mode   = Which waveform to plot:  'ph', 'amp', 'gm' or 'all' (optional.  Default = all).
%
% OUTPUTS:
% out    = Figure handle.

function out=rf_plotWaveform(rf,mode)

if nargin<2
    mode='all';
    
    if nargin<1
        error('ERROR: no input RF pulse specified.  Aborting!!');
    end
    
end

if strcmp(mode,'gm') & ~rf.isGM
    error('ERROR: cannot plot GM function for non-GM pulse.  ABORTING!!');
end


%First make the time vector to use for all plots:
t=cumsum(rf.waveform(:,3));

%Check the mode and plot accordingly:
if strcmp(mode,'ph');
    out=figure;
    set(gcf,'color','w');
    plot(t,rf.waveform(:,1),'LineWidth',1.2);
    set(gca,'FontSize',16);
    xlabel('Time (unitless)','FontSize',20);
    ylabel('Phase (deg.)','FontSize',20);
    box off;
    yl=ylim;
    yrange=(yl(2)-yl(1));
    ylim([yl(1)-0.1*yrange yl(2)+0.1*yrange]);
    xlim([t(1) t(end)]);
elseif strcmp(mode,'amp');
    out=figure;
    set(gcf,'color','w');
    plot(t,rf.waveform(:,2),'LineWidth',1.2);
    set(gca,'FontSize',16);
    xlabel('Time (unitless)','FontSize',20);
    ylabel('Amplitude (arb units)','FontSize',20);
    box off;
    yl=ylim;
    yrange=(yl(2)-yl(1));
    ylim([yl(1)-0.1*yrange yl(2)+0.1*yrange]);
    xlim([t(1) t(end)]);
elseif strcmp(mode,'gm')
    out=figure;
    set(gcf,'color','w');
    plot(t,rf.waveform(:,4),'LineWidth',1.2);
    set(gca,'FontSize',16);
    xlabel('Time (unitless)','FontSize',20);
    ylabel('Gradient (G/cm)','FontSize',20);
    box off;
    yl=ylim;
    yrange=(yl(2)-yl(1));
    ylim([yl(1)-0.1*yrange yl(2)+0.1*yrange]);
    xlim([t(1) t(end)]);
elseif strcmp(mode,'all');
    if rf.isGM
        out=figure('position',[50 50 400 700]);
        set(gcf,'color','w');
        subplot(3,1,1);
        plot(t,rf.waveform(:,2),'LineWidth',1.2);
        set(gca,'FontSize',16);
        xlabel('Time (unitless)','FontSize',20);
        ylabel('Amplitude (arb units)','FontSize',20);
        box off;
        yl=ylim;
        yrange=(yl(2)-yl(1));
        ylim([yl(1)-0.1*yrange yl(2)+0.1*yrange]);
        xlim([t(1) t(end)]);
        
        subplot(3,1,2);
        plot(t,rf.waveform(:,1),'LineWidth',1.2);
        set(gca,'FontSize',16);
        xlabel('Time (unitless)','FontSize',20);
        ylabel('Phase (deg.)','FontSize',20);
        box off;
        yl=ylim;
        yrange=(yl(2)-yl(1));
        ylim([yl(1)-0.1*yrange yl(2)+0.1*yrange]);
        xlim([t(1) t(end)]);
        
        subplot(3,1,3);
        plot(t,rf.waveform(:,4),'LineWidth',1.2);
        set(gca,'FontSize',16);
        xlabel('Time (unitless)','FontSize',20);
        ylabel('Gradient (G/cm)','FontSize',20);
        box off;
        yl=ylim;
        yrange=(yl(2)-yl(1));
        ylim([yl(1)-0.1*yrange yl(2)+0.1*yrange]);
        xlim([t(1) t(end)]);
        
    else
        out=figure('position',[50 50 450 650]);
        set(gcf,'color','w');
        subplot(2,1,1);
        plot(t,rf.waveform(:,2),'LineWidth',1.2);
        set(gca,'FontSize',16);
        xlabel('Time (unitless)','FontSize',20);
        ylabel('Amplitude (arb units)','FontSize',20);
        box off;
        yl=ylim;
        yrange=(yl(2)-yl(1));
        ylim([yl(1)-0.1*yrange yl(2)+0.1*yrange]);
        xlim([t(1) t(end)]);
        
        subplot(2,1,2);
        plot(t,rf.waveform(:,1),'LineWidth',1.2);
        set(gca,'FontSize',16);
        xlabel('Time (unitless)','FontSize',20);
        ylabel('Phase (deg.)','FontSize',20);
        box off;
        yl=ylim;
        yrange=(yl(2)-yl(1));
        ylim([yl(1)-0.1*yrange yl(2)+0.1*yrange]);
        xlim([t(1) t(end)]);
    end
else
    error('ERROR:  mode not recognized.  Choose ''ph'', ''amp'', ''gm'', or ''all''.  ABORTING!!');
end

        
        
    
  
        
    
    

