% run_megapressproc_auto.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% [out1_diff,out1_sum,out1,outw,coilcombos]=run_megapressproc_auto(filestring,coilcombos,avgAlignDomain,alignSS);
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
% outw               = Fully processed water unsuppressed spectrum. 

function [out1_diff,out1_sum,out1_off,out1_on,out1_cc,out1_rm,out1_aa,out1_noproc,nBadAvgTotal1,totalFreqDrift,totalPhaseDrift]=run_megapressproc_auto(filestring,coilcombos,avgAlignDomain,alignSS);

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

%make a new directory for the output report and figures:
mkdir([filestring '/report']);
mkdir([filestring '/report/figs']);

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

%generate unprocessed spectrum:
out1_noproc=op_combinesubspecs(op_averaging(out1_cc),'diff');


%Generate plots showing coil channels before and after phase alignment
%figure('position',[0 50 560 420]);
h=figure('visible','off');
subplot(1,2,1);
plot(raw1_av.ppm,real(raw1_av.specs(:,:,1,1)));xlim([1 5]);
set(gca,'FontSize',8);
set(gca,'XDir','reverse');
xlabel('Frequency (ppm)','FontSize',10);
ylabel('Amplitude (a.u.)','FontSize',10);
title('Before correction','FontSize',12);
box off;
subplot(1,2,2);
plot(raw1_av.ppm,real(spec1_av_pre(:,:,1,1)));xlim([1 5]);
set(gca,'FontSize',12);
set(gca,'XDir','reverse');
xlabel('Frequency (ppm)','FontSize',10);
ylabel('Amplitude(a.u.)','FontSize',10);
title('After correction','FontSize',12);
box off;
set(h,'PaperUnits','centimeters');
set(h,'PaperPosition',[0 0 20 10]);
saveas(h,[filestring '/report/figs/coilReconFig'],'jpg');
saveas(h,[filestring '/report/figs/coilReconFig'],'fig');
close(h);


%%%%%%%%OPTIONAL REMOVAL OF BAD AVERAGES FROM FIRST DATASET%%%%%%%%%%%%%%%%%%%%
close all;
out1_cc2=out1_cc;
nBadAvgTotal1=0;
nbadAverages1=1;
rmbadav1='y';
close all;
if rmbadav1=='n' || rmbadav1=='N'
    out1_rm=out1_cc;
    nsd1='N/A';
else
    sat1='n'
    while sat1=='n' || sat1=='N'
        nsd1=3.2; %Setting the number of standard deviations to 3.2;
        iter=1;
        nbadAverages1=1;
        nBadAvgTotal1=0;
        out1_cc2=out1_cc;
        while nbadAverages1>0
            [out1_rm,metric1{iter},badAverages1]=op_rmbadaverages(out1_cc2,nsd1,'t');
            badAverages1;
            nbadAverages1=length(badAverages1)*raw1.sz(raw1.dims.subSpecs);
            nBadAvgTotal1=nBadAvgTotal1+nbadAverages1;
            out1_cc2=out1_rm;
            iter=iter+1;
            disp([num2str(nbadAverages1) ' bad averages removed on this iteration.']);
            disp([num2str(nBadAvgTotal1) ' bad averages removed in total.']);
            close all;
        end
        %figure('position',[0 50 560 420]);
        %Make figure to show pre-post removal of averages
        h=figure('visible','off');
        subplot(2,2,1);
        plot(out1_cc.ppm,real(out1_cc.specs(:,:,1)));xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Edit-OFF Before','FontSize',12);
        box off;
        subplot(2,2,2);
        plot(out1_rm.ppm,real(out1_rm.specs(:,:,1)));xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Edit-OFF After','FontSize',12);
        box off;
        subplot(2,2,3);  
        plot(out1_cc.ppm,real(out1_cc.specs(:,:,2)));xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Edit-ON Before','FontSize',12);
        box off;
        subplot(2,2,4);
        plot(out1_rm.ppm,real(out1_rm.specs(:,:,2)));xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Edit-ON After','FontSize',12);
        box off;
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 20 15]);
        saveas(h,[filestring '/report/figs/rmBadAvg_prePostFig'],'jpg');
        saveas(h,[filestring '/report/figs/rmBadAvg_prePostFig'],'fig');
        close(h);
        
        %figure('position',[0 550 560 420]);
        h=figure('visible','off');
        plot([1:length(metric1{1})],metric1{1},'.r',[1:length(metric1{iter-1})],metric1{iter-1},'x','MarkerSize',16);
        set(gca,'FontSize',8);
        xlabel('Scan Number','FontSize',10);
        ylabel('Deviation Metric','FontSize',10);
        legend('Before rmBadAv','Before rmBadAv','After rmBadAv','After rmBadAv');
        legend boxoff;
        title('Deviation Metric','FontSize',12);
        box off;
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 20 10]);
        saveas(h,[filestring '/report/figs/rmBadAvg_scatterFig'],'jpg');
        saveas(h,[filestring '/report/figs/rmBadAvg_scatterFig'],'fig');
        close(h);
        
        %sat1=input('are you satisfied with the removal of bad averages? ','s');
        sat1='y';
                     
    end
end


%NOW ALIGN AVERAGES:  A.K.A. Frequency Drift Correction.
driftCorr='y';
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
            
            if driftCorr=='y' || driftCorr=='Y'
                out1_rm2=out1_aa;
            end
        end
        h=figure('visible','off');
        subplot(2,2,1);
        plot(out1_rm.ppm,real(out1_rm.specs(:,:,1)));xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Edit-OFF Before','FontSize',12);
        box off;
        subplot(2,2,2);
        plot(out1_aa.ppm,real(out1_aa.specs(:,:,1)));xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Edit-OFF After','FontSize',12);
        box off;
        subplot(2,2,3);
        plot(out1_rm.ppm,real(out1_rm.specs(:,:,2)));xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Edit-ON Before','FontSize',12);
        box off;
        subplot(2,2,4);
        plot(out1_aa.ppm,real(out1_aa.specs(:,:,2)));xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Edit-ON After','FontSize',12);
        box off;
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 20 15]);
        saveas(h,[filestring '/report/figs/alignAvgs_prePostFig'],'jpg');
        saveas(h,[filestring '/report/figs/alignAvgs_prePostFig'],'fig');
        close(h);
        
        h=figure('visible','off');
        plot([1:out1_aa.sz(out1_aa.dims.averages)],fs1cum,'.-','LineWidth',2);
        set(gca,'FontSize',8);
        xlabel('Scan Number','FontSize',10);
        ylabel('Frequency Drift [Hz]','FontSize',10);
        box off;
        legend('Edit-off scans','Edit-on scans','Location','SouthEast');
        legend boxoff;
        title('Estimated Freqeuncy Drift','FontSize',12);
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 10 10]);
        saveas(h,[filestring '/report/figs/freqDriftFig'],'jpg');
        saveas(h,[filestring '/report/figs/freqDriftFig'],'fig');
        close(h);
        
        h=figure('visible','off');
        plot([1:out1_aa.sz(out1_aa.dims.averages)],phs1cum,'.-','LineWidth',2);
        set(gca,'FontSize',8);
        xlabel('Scan Number','FontSize',10);
        ylabel('Phase Drift [Deg.]','FontSize',10);
        box off;
        legend('Edit-off scans','Edit-on scans');
        legend boxoff;
        title('Estimated Phase Drift','FontSize',12);
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 10 10]);
        saveas(h,[filestring '/report/figs/phaseDriftFig'],'jpg');
        saveas(h,[filestring '/report/figs/phaseDriftFig'],'fig');
        close(h);

        sat1='y';
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
        totalFreqDrift=mean(max(fs1cum)-min(fs1cum));
        totalPhaseDrift=mean(max(phs1cum)-min(phs1cum));
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

%now do automatic zero-order phase correction:
ph0=-phase(out1_ls.fids(1,1))*180/pi;
out1_ph=op_addphase(out1_ls,ph0+180);

%do same phase corection on unprocessed data
out1_noproc=op_addphase(out1_noproc,ph0+180);

%Now align subspecs if desired:
switch alignSS    
    
    case 2
        out1=op_alignMPSubspecs(out1_ph);
        
    case 0
        out1=out1_ph;        
  
    otherwise
        error('ERROR: alignSS value not valid! ');
end


out_filt_diff=op_combinesubspecs(op_filter(out1,5),'diff');

%Make final fully processed data;
out1_diff=op_combinesubspecs(out1,'diff');
out1_sum=op_combinesubspecs(out1,'summ');

out1_off=op_takesubspec(out1,1);
out1_off=op_addphase(out1_off,180);
out1_on=op_takesubspec(out1,2);


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

%Make figure to show the final spectrum:
h=figure('visible','off');
subplot(1,2,1);
plot(out1_off.ppm,out1_off.specs,out1_on.ppm,out1_on.specs,'linewidth',2);xlim([0.2 5.2]);
set(gca,'FontSize',8);
set(gca,'XDir','reverse');
xlabel('Frequency (ppm)','FontSize',10);
ylabel('Amplitude(a.u.)','FontSize',10);
legend('Edit-OFF','Edit-ON');
legend boxoff;
box off;
title('Result: Subspecs','FontSize',12);
subplot(1,2,2);
plot(out1_diff.ppm,out1_diff.specs,'linewidth',2);xlim([0.2 5.2]);
set(gca,'FontSize',8);
set(gca,'XDir','reverse');
xlabel('Frequency (ppm)','FontSize',10);
ylabel('Amplitude(a.u.)','FontSize',10);
legend('diff');
legend boxoff;
box off;
title('Result: Diff Spectrum','FontSize',12);
set(h,'PaperUnits','centimeters');
set(h,'PaperPosition',[0 0 20 10]);
saveas(h,[filestring '/report/figs/finalSpecFig'],'jpg');
saveas(h,[filestring '/report/figs/finalSpecFig'],'fig');

% writ=input('Write? (y or n):  ','s');
writ='y';
if writ=='y' || writ=='Y'
     RF=io_writelcm(out1_diff,[filestring '/' filestring '_diff_lcm'],out1_diff.te);
     RF=io_writelcm(out1_off,[filestring '/' filestring '_editOFF_lcm'],out1_off.te);
     RF=io_writelcm(out1_on,[filestring '/' filestring '_editON_lcm'],out1_on.te);
    if water
        RF=io_writejmrui(outw,[filestring '/MPw.txt']);
    end
end

close all;

%write an html report: 
fid=fopen([filestring '/report/report.html'],'w+');
fprintf(fid,'<!DOCTYPE html>');
fprintf(fid,'\n<html>');
fprintf(fid,'\n<img src= " /Users/jnear/Documents/data/studies/jn_1501_copy/FID-A/FID-A_Documentation/LOGO48.jpg " width="120" height="120"></body>');
fprintf(fid,'\n<h1>FID-A Processing Report</h1>');
fprintf(fid,'\n<h2>Processing pipeline applied to MEGA-PRESS data using run_megapressproc.m</h2>');
fprintf(fid,'\n<p>FILENAME: %s/%s/%s </p>',pwd,filestring,filename1);
fprintf(fid,'\n<p>DATE: %s </p>',date);
fprintf(fid,'\n\n<p> </p>');
fprintf(fid,'\n\n<h2>Results of multi-coil combination:</h2>');
fprintf(fid,'\n<img src= " %s/%s/report/figs/coilReconFig.jpg " width="800" height="400"></body>',pwd,filestring);
fprintf(fid,'\n\n<p> </p>');
fprintf(fid,'\n\n<h2>Results of removal of bad averages:</h2>');
fprintf(fid,'\n<p>Original number of averages: \t%5.6f </p>',raw1.sz(raw1.dims.averages)*raw1.sz(raw1.dims.subSpecs));
fprintf(fid,'\n<p>Number of bad Averages removed:  \t%5.6f </p>',nBadAvgTotal1);
fprintf(fid,'\n<p>Number of remaining averages in processed dataset:  \t%5.6f </p>',out1_rm.sz(out1_rm.dims.averages)*raw1.sz(raw1.dims.subSpecs));
fprintf(fid,'\n<p>Bad Averages Removal Threshold was:  \t%2.2f </p>',nsd1);
fprintf(fid,'\n<img src= " %s/%s/report/figs/rmBadAvg_prePostFig.jpg " width="800" height="600"><img src= " %s/%s/report/figs/rmBadAvg_scatterFig.jpg " width="800" height="400">',pwd,filestring,pwd,filestring);
fprintf(fid,'\n\n<p> </p>');
fprintf(fid,'\n\n<h2>Results of spectral registration:</h2>');
fprintf(fid,'\n<p>Total frequency drift was: \t%5.6f </p>',max(totalFreqDrift));
fprintf(fid,'\n<p>Total phase drift was: \t%5.6f </p>',max(totalPhaseDrift));
fprintf(fid,'\n<img src= " %s/%s/report/figs/alignAvgs_prePostFig.jpg " width="800" height="600">',pwd,filestring);
fprintf(fid,'\n<img src= " %s/%s/report/figs/freqDriftFig.jpg " width="400" height="400"><img src="%s/%s/report/figs/phaseDriftFig.jpg " width="400" height="400">',pwd,filestring,pwd,filestring);
fprintf(fid,'\n\n<p> </p>');
fprintf(fid,'\n\n<h2>Final Result:</h2>');
fprintf(fid,'\n<img src= " %s/%s/report/figs/finalSpecFig.jpg " width="800" height="400">',pwd,filestring);
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



