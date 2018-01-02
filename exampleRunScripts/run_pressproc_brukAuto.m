% run_pressproc_brukAuto.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% [out,out_w,out_noproc,out_w_noproc]=run_pressproc_brukAuto(filestring,filestringw,aaDomain,tmaxin,iterin);
% 
% DESCRIPTION:
% Fully automated processing script for Bruker PRESS MRS data in fid.raw 
% format (raw data).  Processing steps include removal of bad averages, 
% frequency drift correction, leftshifting and auto-phasing. This pipeline 
% requires no user interaction. This function automatically generates an 
% html report to describe the results of each processing step.
% 
% INPUTS:
% filestring    = String variable for the name of the directory containing
%                   the water suppressed fid.raw file.  
% filestringw   = String variable for the name of the directory containing
%                   the water unsuppressed fid.raw file.
% aaDomain      = (Optional) Perform the spectral registration (drift correction) using
%                   the full spectrum ('t'), or only a limited frequency
%                   range ('f').  Default is 'f'.
% tmaxin        = (Optional).  Duration (in sec.) of the time domain signal
%                   used in the spectral registration (drift correction).
%                   Default is 0.2 sec.
% iterin        = (Optional).  Maximum number of allowed iterations for the spectral
%                   registration to converge. Default is 20.
% 
% OUTPUTS:
% out           = Fully processed, water suppressed output spectrum.
% out_w         = Fully processed, water unsuppressed output spectrum.
% out_noproc    = Water suppressed output spectrum without pre-
%                   processing (No bad-averages removal, no frequency drift
%                   correction).
% out_w_noproc  = Water unsuppressed output spectrum without pre-
%                   processing.


function [out,outw,out_noproc,outw_noproc]=run_pressproc_brukAuto(filestring,filestringw,aaDomain,tmaxin,iterin);

if nargin<5
    iterin=20;
    if nargin<4
        tmaxin=0.2;
        if nargin<3
            aaDomain='t';
            if nargin<2
                water=false;
            else
                water=true;
            end
        end
    end
end

%make a new directory for the output report and figures:
mkdir([filestring '/report']);
mkdir([filestring '/report/figs']);


% %read in both datasets:
raw=io_loadspec_bruk(filestring);
if water
    raww=io_loadspec_bruk(filestringw);
else
    outw=struct();
    outw_noproc=struct();
end

%No need to combine coil channels with raw bruker data.
out_cc=raw;
if water
    outw_cc=raww;
end

%generate unprocessed spectrum:
out_noproc=op_averaging(out_cc);
if water
    outw_noproc=op_averaging(outw_cc);
end

%Generate plots showing coil channels before and after phase alignment
%figure('position',[0 50 560 420]);
%{
h=figure('visible','off');
subplot(1,2,1);
plot(raw_av.ppm,real(raw_av.specs(:,:,1)));xlim([1 5]);
set(gca,'FontSize',8);
set(gca,'XDir','reverse');
xlabel('Frequency (ppm)','FontSize',10);
ylabel('Amplitude (a.u.)','FontSize',10);
title('Before correction','FontSize',12);
box off;
subplot(1,2,2);
plot(raw_av.ppm,real(spec_av_pre(:,:,1)));xlim([1 5]);
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
%}


%%%%%%%%OPTIONAL REMOVAL OF BAD AVERAGES FROM FIRST DATASET%%%%%%%%%%%%%%%%%%%%
close all;
out_cc2=out_cc;
nBadAvgTotal=0;
nbadAverages=1;
rmbadav='y';
close all;
if rmbadav=='n' || rmbadav=='N'
    out_rm=out_cc;
    nsd='N/A';
else
    sat='n'
    while sat=='n' || sat=='N'
        nsd=4; %Setting the number of standard deviations;
        iter=1;
        nbadAverages=1;
        nBadAvgTotal=0;
        out_cc2=out_cc;
        while nbadAverages>0
            [out_rm,metric{iter},badAverages]=op_rmbadaverages(out_cc2,nsd,'t');
            badAverages;
            nbadAverages=length(badAverages);
            nBadAvgTotal=nBadAvgTotal+nbadAverages;
            out_cc2=out_rm;
            iter=iter+1;
            disp([num2str(nbadAverages) ' bad averages removed on this iteration.']);
            disp([num2str(nBadAvgTotal) ' bad averages removed in total.']);
            close all;
        end
        %figure('position',[0 50 560 420]);
        %Make figure to show pre-post removal of averages
        h=figure('visible','off');
        subplot(1,2,1);
        plot(out_cc.ppm,real(out_cc.specs(:,:)));xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Before','FontSize',12);
        box off;
        subplot(1,2,2);
        plot(out_rm.ppm,real(out_rm.specs(:,:)));xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('After','FontSize',12);
        box off;
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 20 15]);
        saveas(h,[filestring '/report/figs/rmBadAvg_prePostFig'],'jpg');
        saveas(h,[filestring '/report/figs/rmBadAvg_prePostFig'],'fig');
        close(h);
        
        %figure('position',[0 550 560 420]);
        h=figure('visible','off');
        plot([1:length(metric{1})],metric{1},'.r',[1:length(metric{iter-1})],metric{iter-1},'x','MarkerSize',16);
        set(gca,'FontSize',8);
        xlabel('Scan Number','FontSize',10);
        ylabel('Deviation Metric','FontSize',10);
        legend('Before rmBadAv','After rmBadAv');
        legend boxoff;
        title('Deviation Metric','FontSize',12);
        box off;
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 20 10]);
        saveas(h,[filestring '/report/figs/rmBadAvg_scatterFig'],'jpg');
        saveas(h,[filestring '/report/figs/rmBadAvg_scatterFig'],'fig');
        close(h);
        
        %sat1=input('are you satisfied with the removal of bad averages? ','s');
        sat='y';
                     
    end
end


%NOW ALIGN AVERAGES:  A.K.A. Frequency Drift Correction.
driftCorr='y';
if driftCorr=='n' || driftCorr=='N'
    out_av=op_averaging(out_rm);
    if water
        outw_av=op_averaging(outw_cc);
    end
    fs=0;
    phs=0;
else
    if water
        outw_aa=op_alignAverages(outw_cc,0.2,'n');
    end
    sat='n';
    out_rm2=out_rm;
    while sat=='n' || sat=='N'
        fsPoly=100;
        phsPoly=1000;
        fscum=zeros(out_rm2.sz(out_rm2.dims.averages),1);
        phscum=zeros(out_rm2.sz(out_rm2.dims.averages),1);
        iter=1;
        while (abs(fsPoly(1))>0.001 || abs(phsPoly(1))>0.01) && iter<iterin
            iter=iter+1
            close all
            tmax=0.05+0.01*randn(1);
            ppmmin=1.6+0.1*randn(1);
            ppmmaxarray=[3.5+0.1*randn(1,2),4+0.1*randn(1,3),5.5+0.1*randn(1,1)];
            ppmmax=ppmmaxarray(randi(6,1));
            switch aaDomain
                case 't'
                    [out_aa,fs,phs]=op_alignAverages(out_rm2);
                case 'f'
                    [out_aa,fs,phs]=op_alignAverages_fd(out_rm2,ppmmin,ppmmax,tmax,'y');
                otherwise
                    error('ERROR: avgAlignDomain not recognized!');
            end
            
            fsPoly=polyfit([1:out_aa.sz(out_aa.dims.averages)]',fs,1)
            phsPoly=polyfit([1:out_aa.sz(out_aa.dims.averages)]',phs,1)
            iter
            
            fscum=fscum+fs;
            phscum=phscum+phs;
            
            if driftCorr=='y' || driftCorr=='Y'
                out_rm2=out_aa;
            end
        end
        h=figure('visible','off');
        subplot(1,2,1);
        plot(out_rm.ppm,real(out_rm.specs(:,:)));xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Before','FontSize',12);
        box off;
        subplot(1,2,2);
        plot(out_aa.ppm,real(out_aa.specs(:,:)));xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('After','FontSize',12);
        box off;
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 20 15]);
        saveas(h,[filestring '/report/figs/alignAvgs_prePostFig'],'jpg');
        saveas(h,[filestring '/report/figs/alignAvgs_prePostFig'],'fig');
        close(h);
        
        h=figure('visible','off');
        plot([1:out_aa.sz(out_aa.dims.averages)],fscum,'.-','LineWidth',2);
        set(gca,'FontSize',8);
        xlabel('Scan Number','FontSize',10);
        ylabel('Frequency Drift [Hz]','FontSize',10);
        box off;
        legend('Frequency Drift','Location','SouthEast');
        legend boxoff;
        title('Estimated Freqeuncy Drift','FontSize',12);
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 10 10]);
        saveas(h,[filestring '/report/figs/freqDriftFig'],'jpg');
        saveas(h,[filestring '/report/figs/freqDriftFig'],'fig');
        close(h);
        
        h=figure('visible','off');
        plot([1:out_aa.sz(out_aa.dims.averages)],phscum,'.-','LineWidth',2);
        set(gca,'FontSize',8);
        xlabel('Scan Number','FontSize',10);
        ylabel('Phase Drift [Deg.]','FontSize',10);
        box off;
        legend('Phase Drift','Location','SouthEast');
        legend boxoff;
        title('Estimated Phase Drift','FontSize',12);
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 10 10]);
        saveas(h,[filestring '/report/figs/phaseDriftFig'],'jpg');
        saveas(h,[filestring '/report/figs/phaseDriftFig'],'fig');
        close(h);

        sat='y';
        if sat=='n'
            iter=0;
            p1=100;
            fscum=zeros(out_rm.sz(2:end));
            phscum=zeros(out_rm.sz(2:end));
            fs2cum=zeros(out_cc.sz(2:end));
            phs2cum=zeros(out_cc.sz(2:end));
            out_rm2=out_rm;
            out_cc2=out_cc;
        end
        totalFreqDrift=mean(max(fscum)-min(fscum));
        totalPhaseDrift=mean(max(phscum)-min(phscum));
        close all
    end
    %now combine the averages averages
    out_av=op_averaging(out_aa);
    if water
        outw_av=op_averaging(outw_aa);
    end
end

%now leftshift
out_ls=op_leftshift(out_av,out_av.pointsToLeftshift);
if water
    outw_ls=op_leftshift(outw_av,outw_av.pointsToLeftshift);
end

%now do automatic zero-order phase correction (Use Creatine Peak):
out_ls_zp=op_zeropad(out_ls,16);
index=find(abs(out_ls_zp.specs)==max(abs(out_ls_zp.specs(out_ls_zp.ppm>2.9 & out_ls_zp.ppm<3.1))));
ph0=-phase(out_ls_zp.specs(index))*180/pi;
out_ph=op_addphase(out_ls,ph0);
out_ls_zp=op_addphase(out_ls_zp,ph0);
%And now for water unsuppressed data (use water peak):
if water
    outw_ls_zp=op_zeropad(outw_ls,16);
    indexw=find(abs(outw_ls_zp.specs)==max(abs(outw_ls_zp.specs(outw_ls_zp.ppm>4 & outw_ls_zp.ppm<5.5))));
    ph0w=-phase(outw_ls_zp.specs(indexw))*180/pi;
    outw_ph=op_addphase(outw_ls,ph0w);
    outw_ls_zp=op_addphase(outw_ls_zp,ph0w);
end

%do same phase corection on unprocessed data
out_noproc=op_addphase(op_leftshift(out_noproc,out_noproc.pointsToLeftshift),ph0);
if water
    outw_noproc=op_addphase(op_leftshift(outw_noproc,outw_noproc.pointsToLeftshift),ph0w);
end

%Frequency shift all spectra so that Creatine appears at 3.027 ppm:
[~,frqShift]=op_ppmref(out_ls_zp,2.9,3.1,3.027);
out=op_freqshift(out_ph,frqShift);
out_noproc=op_freqshift(out_noproc,frqShift);
%And now for water unsuppressed data (user water peak and set to 4.65 ppm):
if water
    [~,frqShiftw]=op_ppmref(outw_ls_zp,4,5.5,4.65);
    outw=op_freqshift(outw_ph,frqShiftw);
    outw_noproc=op_freqshift(outw_noproc,frqShiftw);
end

%Make figure to show the final spectrum:
h=figure('visible','off');
plot(out.ppm,out.specs,'linewidth',2);xlim([0.2 5.2]);
set(gca,'FontSize',8);
set(gca,'XDir','reverse');
xlabel('Frequency (ppm)','FontSize',10);
ylabel('Amplitude(a.u.)','FontSize',10);
legend('diff');
legend boxoff;
box off;
title('Result: Final Spectrum','FontSize',12);
set(h,'PaperUnits','centimeters');
set(h,'PaperPosition',[0 0 20 10]);
saveas(h,[filestring '/report/figs/finalSpecFig'],'jpg');
saveas(h,[filestring '/report/figs/finalSpecFig'],'fig');

% writ=input('Write? (y or n):  ','s');
writ='y';
if writ=='y' || writ=='Y'
    mkdir([filestring '/press']);
    [~]=io_writelcm(out,[filestring '/press/press_lcm'],out.te);
    if water
        mkdir([filestring '/press_w']);
        RF=io_writelcm(outw,[filestring '/press_w/press_w_lcm'],outw.te);
    end
end

close all;

%write an html report: 
fid=fopen([filestring '/report/report.html'],'w+');
fprintf(fid,'<!DOCTYPE html>');
fprintf(fid,'\n<html>');
logoPath=which('FID-A_LOGO.jpg');
fprintf(fid,'\n<img src= " %s " width="120" height="120"></body>',logoPath);
fprintf(fid,'\n<h1>FID-A Processing Report</h1>');
fprintf(fid,'\n<h2>Processing pipeline applied to PRESS data using run_pressproc_brukAuto.m</h2>');
fprintf(fid,'\n<p>FILENAME: %s/%s/%s </p>',pwd,filestring);
fprintf(fid,'\n<p>DATE: %s </p>',date);
fprintf(fid,'\n\n<p> </p>');
%fprintf(fid,'\n\n<h2>Results of multi-coil combination:</h2>');
%fprintf(fid,'\n<img src= " %s/%s/report/figs/coilReconFig.jpg " width="800" height="400"></body>',pwd,filestring);
fprintf(fid,'\n\n<p> </p>');
fprintf(fid,'\n\n<h2>Results of removal of bad averages:</h2>');
fprintf(fid,'\n<p>Original number of averages: \t%5.6f </p>',raw.sz(raw.dims.averages));
fprintf(fid,'\n<p>Number of bad Averages removed:  \t%5.6f </p>',nBadAvgTotal);
fprintf(fid,'\n<p>Number of remaining averages in processed dataset:  \t%5.6f </p>',out_rm.sz(out_rm.dims.averages));
fprintf(fid,'\n<p>Bad Averages Removal Threshold was:  \t%2.2f </p>',nsd);
fprintf(fid,'\n<img src= " %s/%s/report/figs/rmBadAvg_prePostFig.jpg " width="800" height="600"><img src= " %s/%s/report/figs/rmBadAvg_scatterFig.jpg " width="800" height="400">',pwd,filestring,pwd,filestring);
fprintf(fid,'\n\n<p> </p>');
fprintf(fid,'\n\n<p> </p>');
fprintf(fid,'\n\n<h2>Results of spectral registration:</h2>');
fprintf(fid,'\n<p>Total frequency drift was: \t%5.6f </p>',max(totalFreqDrift));
fprintf(fid,'\n<p>Total phase drift was: \t%5.6f </p>',max(totalPhaseDrift));
fprintf(fid,'\n<img src= " %s/%s/report/figs/alignAvgs_prePostFig.jpg " width="800" height="600">',pwd,filestring);
fprintf(fid,'\n\n<p> </p>');
fprintf(fid,'\n<img src= " %s/%s/report/figs/freqDriftFig.jpg " width="400" height="400"><img src="%s/%s/report/figs/phaseDriftFig.jpg " width="400" height="400">',pwd,filestring,pwd,filestring);
fprintf(fid,'\n\n<p> </p>');
fprintf(fid,'\n\n<h2>Final Result:</h2>');
fprintf(fid,'\n<img src= " %s/%s/report/figs/finalSpecFig.jpg " width="800" height="400">',pwd,filestring);
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



