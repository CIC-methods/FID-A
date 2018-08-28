% op_plotspec.m
% Jamie Near, McGill University 2015.
% 
% USAGE:
% out=op_plotspec(in,ppmmin,ppmmax,xlab,ylab,tit);
% 
% DESCRIPTION:
% Plot the MR spectrum in the frequency domain.  
% 
% INPUTS:
% in     = input data in matlab structure format.  This argument can be a
%          MATLAB structure (formatted according to the FID-A structure
%          format), or it can be a cell array, where the elemenets of each
%          cell are FID-A structures.  In the latter case, each spectrum in
%          the structure will be plotted, with the option of a vertical
%          offset.
% ppmmin = lower limit of ppm scale to plot (optional.  Default = 0.2 ppm).
% ppmmax = upper limit of ppm scale to plot (optional.  Default = 5.2 ppm).
% xlab   = Label for the x-axis (optional.  Default = 'Frequency (ppm)');
% ylab   = label for the y-axis (optional.  Default = '');
% tit    = label for the title of the plot (optional.  Default = '');
%
% OUTPUTS:
% out    = Figure handle.

function out=op_plotspec(in,ppmmin,ppmmax,xlab,ylab,tit)

if nargin<6
    tit='';
    if nargin<5
        ylab='';
        if nargin<4
            xlab='Frequency (ppm)';
            if nargin<3
                ppmmax=5.2;
                if nargin<2
                    ppmmin=0.2;
                    if nargin<1
                        error('ERROR: no input spectrum specified.  Aborting!!');
                    end
                end
            end
        end
    end
end


if isstruct(in)
    if ndims(squeeze(in.fids))>2
        disp('More than 1 dimension detected in addition to the frequency dimension.');
        disp('Which dimension would you like to be displayed? ');
        if in.dims.coils
            disp(['     Enter: ''' num2str(in.dims.coils) ''' to display signal from all RF channels.']);
        end
        if in.dims.averages
            disp(['     Enter: ''' num2str(in.dims.averages) ''' to display signal from all averages.']);
        end
        if in.dims.subSpecs
            disp(['     Enter: ''' num2str(in.dims.subSpecs) ''' to display signal from all subspectra.']);
        end
        if in.dims.extras
            disp(['     Enter: ''' num2str(in.dims.extras) ''' to display signal from the ''extras'' dimension.']);
        end
        dim=input('............ ');
    else
        dim=0;
    end
    fignum=figure;
    if ~dim
        out=plot(in.ppm,real(squeeze(in.specs)));
    elseif dim==2
        out=plot(in.ppm,real(squeeze(in.specs(:,:,1))));
    elseif dim==3
        out=plot(in.ppm,real(squeeze(in.specs(:,1,:,1))));
    elseif dim==4
        out=plot(in.ppm,real(squeeze(in.specs(:,1,1,:,1))));
    elseif dim==5 
        out=plot(in.ppm,real(squeeze(in.specs(:,1,1,1,:))));
    end
    xlim([ppmmin ppmmax]);
    set(gca,'XDir','reverse');
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
    Fig1Ax1 = get(fignum, 'Children');
    Fig1Ax1Line1 = get(Fig1Ax1, 'Children');
    set(Fig1Ax1Line1, 'LineWidth', 1.2);
    %set(Fig1Ax1Line1,'MarkerSize', 10);
elseif iscell(in)
    if ndims(squeeze(in{1}.specs))>2
        disp('More than 1 dimension detected in addition to the frequency dimension.');
        disp('Which dimension would you like to be displayed? ');
        if in{1}.dims.coils
            disp(['     Enter: ''' num2str(in{1}.dims.coils) ''' to display signal from all RF channels.']);
        end
        if in{1}.dims.averages
            disp(['     Enter: ''' num2str(in{1}.dims.averages) ''' to display signal from all averages.']);
        end
        if in{1}.dims.subSpecs
            disp(['     Enter: ''' num2str(in{1}.dims.subSpecs) ''' to display signal from all subspectra.']);
        end
        if in{1}.dims.extras
            disp(['     Enter: ''' num2str(in{1}.dims.extras) ''' to display signal from the ''extras'' dimension.']);
        end
        dim=input('............ ');
    else
        dim=0;
    end
    figure;
    if ~dim
        out=plot(in{1}.ppm,real(in{1}.specs));
    elseif dim==2
        out=plot(in{1}.ppm,real(squeeze(in{1}.specs(:,:,1))));
    elseif dim==3
        out=plot(in{1}.ppm,real(squeeze(in{1}.specs(:,1,:,1))));
    elseif dim==4
        out=plot(in{1}.ppm,real(squeeze(in{1}.specs(:,1,1,:,1))));
    elseif dim==5 
        out=plot(in{1}.ppm,real(squeeze(in{1}.specs(:,1,1,1,:))));
    end
    ylabel('ARB UNITS');
    disp('Multiple input spectra detected!! ')
    stagger=input('Please enter the desired vertical spacing of the spectra in ARB UNITS:  ');
    close;
    fignum=figure;
    hold
    colours=distinguishable_colors(length(in));
    for n=1:length(in)
        if ~dim
            out=plot(in{n}.ppm,real(in{n}.specs)+(n-1)*stagger,'Color',colours(n,:));
        elseif dim==2
            out=plot(in{n}.ppm,real(squeeze(in{n}.specs(:,:,1)))+(n-1)*stagger,'Color',colours(n,:));
        elseif dim==3
            out=plot(in{n}.ppm,real(squeeze(in{n}.specs(:,1,:,1)))+(n-1)*stagger,'Color',colours(n,:));
        elseif dim==4
            out=plot(in{n}.ppm,real(squeeze(in{n}.specs(:,1,1,:,1)))+(n-1)*stagger,'Color',colours(n,:));
        elseif dim==5 
            out=plot(in{n}.ppm,real(squeeze(in{n}.specs(:,1,1,1,:)))+(n-1)*stagger,'Color',colours(n,:));
        end
        xlim([ppmmin ppmmax]);
        set(gca,'XDir','reverse');
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
        Fig1Ax1 = get(fignum, 'Children');
        Fig1Ax1Line1 = get(Fig1Ax1, 'Children');
        set(Fig1Ax1Line1, 'LineWidth', 1.2);
    end
else
    error('ERROR:  Input data format not recognized.  ABORTING!!');
end


        
    
    

