% op_plotfid.m
% Jamie Near, McGill University 2015.
% 
% USAGE:
% out=op_plotfid(in,tmax,xlab,ylab,tit);
% 
% DESCRIPTION:
% Plot the spectrum in the time domain.
% 
% INPUTS:
% in     = input data in matlab structure format.  This argument can be a
%          MATLAB structure (formatted according to the FID-A structure
%          format), or it can be a cell array, where the elemenets of each
%          cell are FID-A structures.  In the latter case, each spectrum in
%          the structure will be plotted, with the option of a vertical
%          offset.
% tmax   = upper limit of time scale to plot in seconds (optional.  Default = max(in.t)).
% xlab   = Label for the x-axis (optional.  Default = 'Time (sec)');
% ylab   = label for the y-axis (optional.  Default = 'FID Amplitude (arb units)');
% tit    = label for the title of the plot (optional.  Default = '');

function out=op_plotfid(in,tmax,xlab,ylab)

if nargin<5
    tit='';
    if nargin<4
        ylab='FID Amplitude (arb units)';
        if nargin<3
            xlab='Time (sec)';
            if nargin<2
                if isstruct(in)
                    tmax=max(in.t);
                elseif iscell(in)
                    tmax=max(in{1}.t);
                else
                    error('ERROR:  input format not recognized.  Aborting!!');
                end
                    if nargin<1
                        error('ERROR: no input spectrum specified.  Aborting!!');
                    end
            end
        end
    end
end


if isstruct(in)
    out=plot(in.t,real(in.fids));
    xlim([0 tmax]);
    set(gca,'FontSize',16);
    if isempty(ylab)
        set(gca,'YColor','w');
    else
        set(gca,'YColor','k');
    end
    set(gca,'XColor','k');
    set(gca,'Color','w');
    set(gcf,'Color','w');
    box off;
    title(tit);
    xlabel(xlab,'FontSize',20);
    ylabel(ylab,'FontSize',20);
    Fig1Ax1 = get(1, 'Children');
    Fig1Ax1Line1 = get(Fig1Ax1, 'Children');
    set(Fig1Ax1Line1, 'LineWidth', 1.2);
    %set(Fig1Ax1Line1,'MarkerSize', 10);
elseif iscell(in)
    plot(in{1}.t,real(in{1}.fids));
    ylabel('ARB UNITS');
    stagger=input('Please enter the desired vertical spacing of the fids in ARB UNITS:  ');
    close;
    figure
    hold
    colours=distinguishable_colors(length(in));
    for n=1:length(in)
        out=plot(in{n}.t,real(in{n}.fids)+(n-1)*stagger,'Color',colours(n,:));
        xlim([0 tmax]);
        set(gca,'FontSize',16);
        if isempty(ylab)
            set(gca,'YColor','w');
        else
            set(gca,'YColor','k');
        end
        set(gca,'XColor','k');
        set(gca,'Color','w');
        set(gcf,'Color','w');
        box off;
        title(tit);
        xlabel(xlab,'FontSize',20);
        ylabel(ylab,'FontSize',20);
        Fig1Ax1 = get(1, 'Children');
        Fig1Ax1Line1 = get(Fig1Ax1, 'Children');
        set(Fig1Ax1Line1, 'LineWidth', 1.2);
    end
else
    error('ERROR:  Input data format not recognized.  ABORTING!!');
end
