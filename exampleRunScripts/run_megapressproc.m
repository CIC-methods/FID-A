% run_megapressproc.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% [out1_diff,out1_sum,out1,outw,coilcombos]=run_megapressproc(filestring,coilcombos,avgAlignDomain,alignSS);
% 
% DESCRIPTION:
% Processing script for Siemens MEGA-PRESS MRS data in .dat format (twix 
% raw data).  Includes combination of reciever channels, removal of bad 
% averages, freqeuncy drift correction, manual alignment of edit-on and 
% edit-off spectra, and leftshifting.
% 
% INPUTS:
% filestring         = String variable for the name of the directory containing
%                     the water suppressed .dat file.  Water unsuppressed
%                     .dat file should be contained in [filestring '_w/'];
% coilcombos         = (Optional).  A structure obtained by running the
%                      op_getcoilcombos function.  This allows the user to 
%                      specify the coil phases and amplitudes as an input, 
%                      rather calculating these from the input data by default.  
% avgAlignDomain     = (Optional) Perform the spectral registration (drift correction) using
%                     the full spectrum ('t'), or only a limited frequency
%                     range ('f').  Default is 'f'.
% alignSS            = (Optional)
%                     0 - Do not align the edit-on and edit-off subspectra (default).
%                     2 - Perform manual alignment of edit-on and edit-off subspectra.
% 
% OUTPUTS:
% out1_diff          = Fully processed difference spectrum.
% out1_sum           = Fully processed sum spectrum.
% out1               = Fully processed edit-on and edit-off subspectra.
% outw               = Fully processed water unsuppressed spectrum. 
% coilcombos         = Estimated coil weights and phases.

function [out1_diff,out1_sum,out1,outw,coilcombos]=run_megapressproc(filestring,coilcombos,avgAlignDomain,alignSS);

if nargin<4
    alignSS=2;
    if nargin<3
        avgAlignDomain='f';
        if nargin<2
            ccGiven=false;
        else
            ccGiven=true;
        end
    end
end

%Find the filename of the first MEGA_GABA dataset
close all
unixString1=[filestring '/*.dat'];
[filename1]=dir(unixString1);
filename1=filename1.name(1:end);

unixStringw=['ls ' filestring(1:end) '_w/*.dat'];
[status,filenamew]=unix(unixStringw);
filenamew=filenamew(1:end-1);

if exist(filenamew);
    water=true
else
    water=false
end


% %read in both datasets:
raw1=io_loadspec_twix([filestring '/' filename1]);
if water
    raww=io_loadspec_twix(filenamew);
end

%first step should be to combine coil channels.  To do this find the coil
%phases from the water unsuppressed data.
if water
    if ~ccGiven
        coilcombos=op_getcoilcombos(raww,1);
    end
    [outw_cc,fidw_pre,specw_pre,phw,sigw]=op_addrcvrs(raww,1,'w',coilcombos);
else
    if ~ccGiven
        coilcombos=op_getcoilcombos(op_averaging(op_combinesubspecs(raw1,'summ')),1);
    end
    
end
[out1_cc,fid1_pre,spec1_pre,ph1,sig1]=op_addrcvrs(raw1,1,'w',coilcombos);
[out1_av_cc,fid1_av_pre,spec1_av_pre]=op_addrcvrs(op_averaging(raw1),1,'w',coilcombos);
raw1_av=op_averaging(raw1);

%plot coil channels before and after phase alignment
figure('position',[0 50 560 420]);
subplot(2,1,1);
plot(raw1_av.ppm,real(raw1_av.specs(:,:,1,1)));xlim([1 5]);
set(gca,'XDir','reverse');
xlabel('Frequency (ppm)');
ylabel('Amplitude (a.u.)');
title('Multi-channel data before phase correction');
subplot(2,1,2);
plot(raw1_av.ppm,real(spec1_av_pre(:,:,1,1)));xlim([1 5]);
set(gca,'XDir','reverse');
xlabel('Frequency (ppm)');
ylabel('Amplitude(a.u.)');
title('Multi-channel data after phase correction');

if water
    figure('position',[0 550 560 420]);
    subplot(2,1,1);
    plot(raww.ppm,real(raww.specs(:,:,1,1)));xlim([1 5]);
    set(gca,'XDir','reverse');
    xlabel('Frequency (ppm)');
    ylabel('Amplitude (a.u.)');
    title('Multi-channel data before phase correction');
    subplot(2,1,2);
    plot(raww.ppm,real(specw_pre(:,:,1,1)));xlim([1 5]);
    set(gca,'XDir','reverse');
    xlabel('Frequency (ppm)');
    ylabel('Amplitude(a.u.)');
    title('Multi-channel data after phase correction');
end
disp('The figure shows multi coil elements, before and after phase alignment.');
disp('Hit "Enter" to continue.');
pause;

%%%%%%%%OPTIONAL REMOVAL OF BAD AVERAGES FROM FIRST DATASET%%%%%%%%%%%%%%%%%%%%
close all;
figure('position',[0 50 560 420]);
subplot(1,2,1);
%if out1_cc.dims.subSpecs==3
    plot(out1_cc.ppm,real(out1_cc.specs(:,:,1)));xlim([1 5]);
%elseif out1_cc.dims.subSpecs==2
%   plot(out1_cc.ppm,squeeze(out1_cc.specs(:,1,:)));xlim([1 5]);
%end
set(gca,'XDir','reverse');
xlabel('Frequency (ppm)');
ylabel('Amplitude(a.u.)');
title('Edit OFF Spectra (All averages)');
subplot(1,2,2);
plot(out1_cc.ppm,real(out1_cc.specs(:,:,2)));xlim([1 5]);
set(gca,'XDir','reverse');
xlabel('Frequency (ppm)');
ylabel('Amplitude(a.u.)');
title('Edit ON Spectra (All averages)');

out1_cc2=out1_cc;
nBadAvgTotal1=0;
nbadAverages1=1;
rmbadav1=input('would you like to remove bad averages from this dataset?  ','s');
close all;
if rmbadav1=='n' || rmbadav1=='N'
    out1_rm=out1_cc;
    nsd1='N/A';
else
    sat1='n'
    while sat1=='n' || sat1=='N'
        nsd1=input('input number of standard deviations (around 3 is recommended): ');
        iter=1;
        nbadAverages1=1;
        nBadAvgTotal1=0;
        out1_cc2=out1_cc;
        while nbadAverages1>0
            %sat1='n'
            %while sat1=='n'||sat1=='N'
            %nsd1=input('input number of standard deviations: ');
            [out1_rm,metric1{iter},badAverages1]=op_rmbadaverages(out1_cc2,nsd1,'t');
            badAverages1;
            nbadAverages1=length(badAverages1)*raw1.sz(raw1.dims.subSpecs);
            %figure;
            %subplot(1,2,1);
            %plot(out1_rm.ppm,out1_rm.specs(:,:,1));xlim([0 5]);
            %subplot(1,2,2);
            %plot(out1_rm.ppm,out1_rm.specs(:,:,2));xlim([0 5]);
            
            %sat1=input('are you satisfied with the removal of the bad averages? ','s');
            %close all;
            
            nBadAvgTotal1=nBadAvgTotal1+nbadAverages1;
            out1_cc2=out1_rm;
            iter=iter+1;
            disp([num2str(nbadAverages1) ' bad averages removed on this iteration.']);
            disp([num2str(nBadAvgTotal1) ' bad averages removed in total.']);
            disp('Press Enter to continue');
            pause;
            close all;
        end
        figure('position',[0 50 560 420]);
        subplot(2,2,1);
        plot(out1_cc.ppm,real(out1_cc.specs(:,:,1)));xlim([1 5]);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)');
        ylabel('Amplitude(a.u.)');
        title('Edit OFF specra prior to Removal of Bad Averages');
        subplot(2,2,2);
        plot(out1_rm.ppm,real(out1_rm.specs(:,:,1)));xlim([1 5]);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)');
        ylabel('Amplitude(a.u.)');
        title('Edit OFF specra after Removal of Bad Averages');
        subplot(2,2,3);  
        plot(out1_cc.ppm,real(out1_cc.specs(:,:,2)));xlim([1 5]);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)');
        ylabel('Amplitude(a.u.)');
        title('Edit ON specra prior to Removal of Bad Averages');
        subplot(2,2,4);
        plot(out1_rm.ppm,real(out1_rm.specs(:,:,2)));xlim([1 5]);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)');
        ylabel('Amplitude(a.u.)');
        title('Edit ON specra after Removal of Bad Averages');
        figure('position',[0 550 560 420]);
        plot([1:length(metric1{1})],metric1{1},'.r',[1:length(metric1{iter-1})],metric1{iter-1},'x');
        xlabel('Scan Number');
        ylabel('Unlikeness Metric');
        legend('Before rmBadAv','Before rmBadAv','After rmBadAv','After rmBadAv');
        title('Metric for rejection of motion corrupted scans');
        sat1=input('are you satisfied with the removal of bad averages? ','s');
                     
    end
end

%write a readme file to record the number of dropped avgs
fid=fopen([filestring '/readme_mp1.txt'],'w+')
fprintf(fid,'Original number of averages: \t%5.6f',raw1.sz(raw1.dims.averages)*raw1.sz(raw1.dims.subSpecs));
fprintf(fid,'\nNumber of bad Averages removed:  \t%5.6f',nBadAvgTotal1);
fprintf(fid,'\nNumber of remaining averages in processed dataset:  \t%5.6f',out1_rm.sz(out1_rm.dims.averages)*raw1.sz(raw1.dims.subSpecs));
fprintf(fid,'\nBad Averages Removal Threshold was:  \t%2.2f',nsd1);
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%NOW ALIGN AVERAGES:  A.K.A. Frequency Drift Correction.
driftCorr=input('\n\nWould you like to perform the frequency drift correction?  ','s');
if driftCorr=='n' || driftCorr=='N'
    out1_av=op_averaging(out1_rm);
    if water
        outw_av=op_averaging(outw_cc);
    end
    fs1=0;
    phs1=0;
else
    if water
        outw_aa=op_alignAverages(outw_cc,0.2,'n');
    end
    sat1='n';
    out1_rm2=out1_rm;
    while sat1=='n' || sat1=='N'
        iter=0;
        iterin=20;
        p1=100;
        fs1cum=zeros(out1_rm.sz(2:end));
        phs1cum=zeros(out1_rm.sz(2:end));
        while (abs(p1(1))>0.0003 && iter<iterin)
            iter=iter+1
            %sat='n';
            %while sat=='n' || sat=='N'
            close all
            tmax=0.25+0.03*randn(1);
            ppmmin=1.6+0.1*randn(1);
            ppmmaxarray=[3.5+0.1*randn(1,2),4+0.1*randn(1,3),5.5+0.1*randn(1,1)];
            ppmmax=ppmmaxarray(randi(6,1));
            switch avgAlignDomain
                case 't'
                    [out1_aa,fs1,phs1]=op_alignAverages(out1_rm2,tmax,'y');
                case 'f'
                    [out1_aa,fs1,phs1]=op_alignAverages_fd(out1_rm2,ppmmin,ppmmax,tmax,'y');
                otherwise
                    error('ERROR: avgAlignDomain not recognized!');
            end
            
            x1=repmat([1:size(fs1,1)]',1,out1_aa.sz(out1_aa.dims.subSpecs));
            p1=polyfit(x1,fs1,1)
            
            fs1cum=fs1cum+fs1;
            phs1cum=phs1cum+phs1;
            
            %plot the data before and after aligning Averages:
            %             figure('position',[0 0 560 420]);
            %             subplot(2,2,1);
            %             plot(out1_rm.ppm,out1_rm.specs(:,:,1));xlim([0 6]);
            %             subplot(2,2,2);
            %             plot(out1_aa.ppm,out1_aa.specs(:,:,1));xlim([0 6]);
            %             subplot(2,2,3);
            %             plot(out1_rm.ppm,out1_rm.specs(:,:,2));xlim([0 6]);
            %             subplot(2,2,4);
            %             plot(out1_aa.ppm,out1_aa.specs(:,:,2));xlim([0 6]);
            
            %             figure('position',[0 500 560 420]);
            %             plot([1:out1_aa.sz(out1_aa.dims.averages)],fs1,...
            %                 '.',x1,polyval(p1,x1));
            
            %             plot([1:out1_aa.sz(out1_aa.dims.averages)],fs1,...
            %                 [1:out2_aa.sz(out2_aa.dims.averages)],fs2,...
            %                 '.',x1,polyval(p1,x1),x2,polyval(p2,x2));
            
            %sat=input('are you satisfied with the drift correction? ','s');
            %end
            %driftCorr=input('would you like to perform drift correction again? ','s');
            if driftCorr=='y' || driftCorr=='Y'
                out1_rm2=out1_aa;
            end
        end
        figure('position',[0 50 560 420]);
        subplot(2,2,1);
        plot(out1_rm.ppm,real(out1_rm.specs(:,:,1)));xlim([1 5]);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)');
        ylabel('Amplitude(a.u.)');
        title('Edit ON specra prior to drift correction');
        subplot(2,2,2);
        plot(out1_aa.ppm,real(out1_aa.specs(:,:,1)));xlim([1 5]);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)');
        ylabel('Amplitude(a.u.)');
        title('Edit ON specra after drift correction');
        subplot(2,2,3);
        plot(out1_rm.ppm,real(out1_rm.specs(:,:,2)));xlim([1 5]);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)');
        ylabel('Amplitude(a.u.)');
        title('Edit OFF specra prior to drift correction');
        subplot(2,2,4);
        plot(out1_aa.ppm,real(out1_aa.specs(:,:,2)));xlim([1 5]);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)');
        ylabel('Amplitude(a.u.)');
        title('Edit OFF specra after drift correction');
        
        figure('position',[0 550 1130 420]);
        plot([1:out1_aa.sz(out1_aa.dims.averages)],fs1cum,'.-');
        xlabel('Scan Number');
        ylabel('Frequency Drift [Hz]');
        legend('Edit-on scans','Edit-off scans','Location','SouthEast');
        legend boxoff;
        title('Estimated Freqeuncy Drift');
        
        figure('position',[600 50 560 420]);
        plot([1:out1_aa.sz(out1_aa.dims.averages)],phs1cum,'.-');
        xlabel('Scan Number');
        ylabel('Phase Drift [Deg.]');
        legend('Edit-on scans','Edit-off scans');
        title('Estimated Phase Drift');

        if water
            if outw_cc.dims.averages>0
                figure('position',[1140 50 560 400]);
                subplot(2,1,1);
                if outw_cc.dims.subSpecs
                    plot(outw_cc.ppm,real(outw_cc.specs(:,:,1)),outw_cc.ppm,real(outw_cc.specs(:,:,2)));xlim([3 7]);
                else
                    plot(outw_cc.ppm,real(outw_cc.specs(:,:)));xlim([3 7]);
                end
                set(gca,'XDir','reverse');
                xlabel('Frequency (ppm)');
                ylabel('Amplitude(a.u.)');
                title('Water Unsuppressed spectrum prior to drift correction');
                subplot(2,1,2);
                if outw_aa.dims.subSpecs
                    plot(outw_aa.ppm,real(outw_aa.specs(:,:,1)),outw_aa.ppm,real(outw_aa.specs(:,:,2)));xlim([3 7]);
                else
                    plot(outw_aa.ppm,real(outw_aa.specs(:,:)));xlim([3 7]);
                end
                set(gca,'XDir','reverse');
                xlabel('Frequency (ppm)');
                ylabel('Amplitude(a.u.)');
                title('Water Unsuppressed spectrum after drift correction');
            end
        end


        sat1=input('are you satisfied with drift correction?  ','s');
        if sat1=='n'
            iter=0;
            p1=100;
            fs1cum=zeros(out1_rm.sz(2:end));
            phs1cum=zeros(out1_rm.sz(2:end));
            fs2cum=zeros(out1_cc.sz(2:end));
            phs2cum=zeros(out1_cc.sz(2:end));
            out1_rm2=out1_rm;
            out1_cc2=out1_cc;
        end
        
        close all
    end
    %now combine the averages averages
    out1_av=op_averaging(out1_aa);
    if water
        outw_av=op_averaging(outw_aa);
    end
end

%now leftshift
out1_ls=op_leftshift(out1_av,out1_av.pointsToLeftshift);
if water
    outw_ls=op_leftshift(outw_av,outw_av.pointsToLeftshift);
end

%now choose global phase correction:
SpecTool(out1_ls,0.5,0,7);
disp('******************************************************************************');
disp('Use GUI Interface to adjust 0th order phase until NAA peak (2 ppm) is inverted');
disp('***NOTE If you are using the Siemens MEGA_PRESS WIP (WIP529), then you should');
disp('adjust the 0th order phase until the NAA peak is upright!***');
disp('************************************************************');
fprintf('\n');
ph0=input('Input desired phase shift (Hz): ');

switch alignSS    
    case 2
        figure('position',[0 50 560 420]);
        out1_ls_filt=op_addphase(op_filter(out1_ls,5),ph0);
        subSpecTool(out1_ls_filt,0,7);
        disp('***************************************************************************************');
        disp('Use GUI interface to align edit-ON and edit-OFF scans by adjusting Phase and Frequency.');
        disp('Try to minimize the residual water, residual Creatine, and residual Choline peaks!');
        disp('***NOTE If you are using the Siemens MEGA_PRESS WIP (WIP529), then you will');
        disp('have to add about 180 degrees of phase to the subspectrum!***');
        disp('*************************************************************');
        fprintf('\n');
        phshft1=input('Input Desired Phase Shift (Degrees) for first spectrum: ');
        frqshft1=input('Input Desired Frequncy Shift (Hz) for first spectrum: ');
        out1=op_freqshiftSubspec(op_addphaseSubspec(op_addphase(out1_ls,ph0),phshft1),frqshft1);
        close all;
        
    case 0
        out1=op_addphase(out1_ls,ph0);        
        
    otherwise
        error('ERROR: alignSS value not valid! ');
end


out_filt_diff=op_combinesubspecs(op_filter(out1,5),'diff');
close all;
figure('position',[0 50 560 420]);
plot(out_filt_diff.ppm,real(out_filt_diff.specs));
set(gca,'XDir','reverse');


%Make final fully processed data;
out1_diff=op_combinesubspecs(out1,'diff');
out1_sum=op_combinesubspecs(out1,'summ');



%Make final water unsuppressed data
if water
    if ~isempty(findstr(outw_ls.seq,'edit_529'));
        if outw_ls.dims.subSpecs
            outw=op_combinesubspecs(outw_ls,'diff');
        else
            outw=outw_ls;
        end
    else
        if outw_ls.dims.subSpecs
            outw=op_combinesubspecs(outw_ls,'summ');
        else
            outw=outw_ls;
        end
    end
    outw=op_addphase(outw,-phase(outw.fids(1))*180/pi,0,4.65,1);
else
    outw=0;
end


writ=input('Write? (y or n):  ','s');
if writ=='y' || writ=='Y'
    RF=io_writejmrui(out1_diff,[filestring '/MPdiff_full.txt']);
    RF=io_writejmrui(out1_sum,[filestring '/MPsum_full.txt']);
    if water
        RF=io_writejmrui(outw,[filestring '/MPw.txt']);
    end
end



