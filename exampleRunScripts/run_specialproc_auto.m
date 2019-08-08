% run_specialproc_auto.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% [out,out_w,out_noproc,out_w_noproc]=run_specialproc_auto(filestring,aaDomain,tmaxin,iterin);
% 
% DESCRIPTION:
% Processing script for Siemens SPECIAL MRS data in .dat format (twix raw 
% data).  Includes combination of reciever channels, removal of bad averages, 
% frequency drift correction, and leftshifting.
% 
% INPUTS:
% filestring   = String variable for the name of the directory containing
%                   the water suppressed .dat file.  Water unsuppressed
%                   .dat file should be contained in [filestring '_w/'];
% aaDomain     = (Optional) Perform the spectral registration (drift correction) using
%                   the full spectrum ('t'), or only a limited frequency
%                   range ('f').  Default is 'f'.
% tmaxin       = (Optional).  Duration (in sec.) of the time domain signal
%                   used in the spectral registration (drift correction).
%                   Default is 0.2 sec.
% iterin       = (Optional).  Maximum number of allowed iterations for the spectral
%                   registration to converge. Default is 20.
% 
% OUTPUTS:
% out          = Fully processed, water suppressed output spectrum.
% out_w        = Fully processed, water unsuppressed output spectrum.
% out_noproc   = Water suppressed output spectrum without pre-
%                   processing (No bad-averages removal, no frequency drift
%                   correction).
% out_w_noproc = Water unsuppressed output spectrum without pre-
%                   processing.

function [out,out_w,out_noproc,out_w_noproc]=run_specialproc_auto(filestring,aaDomain,tmaxin,iterin)

if nargin<4
    iterin=20;
    if nargin<3
        tmaxin=0.2;
        if nargin<2
            aaDomain='f';
        end
    end
end

%Ensure that filestring is an absolute path, otherwise you can get yourself
%in a tangle later inserting figures at the report stage.
if ~java.io.File(filestring).isAbsolute
    filestring = cd(cd(filestring)); 
end

%make a new directory for the output report and figures:
reportDir = fullfile(filestring,'report');
mkdir(reportDir);
reportFigsDir = fullfile(filestring,'report','figs');
mkdir(reportFigsDir);

%Find the filename of the SPECIAL dataset
close all
unixString1=fullfile(filestring,'*.dat');
[filename]=dir(unixString1);
filename=filename.name(1:end);

unixStringw=fullfile([filestring '_w'],'*.dat');
[filenamew]=dir(unixStringw);
if ~isempty(filenamew)
    filenamew=filenamew.name(1:end);
    water=true
else
    water=false
    out_w=struct();
    out_w_noproc=struct();
end

%read in the data:
out_raw=io_loadspec_twix(fullfile(filestring,filename));

%load water unsuppressed data and find the coil phases:
if water
    disp('***FOUND WATER UNSUPPRESSED DATA***');
    out_w_raw=io_loadspec_twix(fullfile([filestring '_w/'], filenamew));
    coilcombos=op_getcoilcombos(op_combinesubspecs(out_w_raw,'diff'),2);
else
    disp('***WATER UNSUPPRESSED DATA NOT FOUND***');
    coilcombos=op_getcoilcombos_specReg(op_combinesubspecs(op_averaging(out_raw),'diff'),0,0.01,2);
    out_w='';
    out_w_noproc='';
end

%now combine the coil channels:
[out_cc,fid_pre,spec_pre,ph,sig]=op_addrcvrs(out_raw,2,'w',coilcombos);
if water
    [out_w_cc,fid_w_pre,spec_w_pre,ph_w,sig_w]=op_addrcvrs(out_w_raw,2,'w',coilcombos);
end

%make the un-processed spectra, which may be optionally output from this function:
out_noproc=op_combinesubspecs(op_averaging(out_cc),'diff');
if water
    out_w_noproc=op_combinesubspecs(op_averaging(out_w_cc),'diff');
end

%plot the data before and after coil phasing:
%figure('position',[0 50 560 420]);
h=figure('visible','off');
subplot(1,2,1);
plot(out_raw.ppm,out_raw.specs(:,:,1,1));xlim([-1 7]);
set(gca,'FontSize',12);
set(gca,'XDir','reverse');
xlabel('Frequency (ppm)','FontSize',10);
ylabel('Amplitude (a.u.)','FontSize',10);
title('Before correction','FontSize',12);
box off;
subplot(1,2,2);
plot(out_raw.ppm,spec_pre(:,:,1,1));xlim([-1 7]);
set(gca,'FontSize',12);
set(gca,'XDir','reverse');
xlabel('Frequency (ppm)','FontSize',10);
ylabel('Amplitude(a.u.)','FontSize',10);
title('After correction','FontSize',12);
box off;
set(h,'PaperUnits','centimeters');
set(h,'PaperPosition',[0 0 20 10]);
saveas(h,fullfile(reportFigsDir,'coilReconFig'),'jpg');
saveas(h,fullfile(reportFigsDir,'coilReconFig'),'fig');
close(h);

if water
    %figure('position',[0 550 560 420]);
    h=figure('visible','off');
    subplot(1,2,1);
    plot(out_w_raw.ppm,out_w_raw.specs(:,:,1,1));xlim([4 5]);
    set(gca,'FontSize',12);
    set(gca,'XDir','reverse');
    xlabel('Frequency (ppm)','FontSize',10);
    ylabel('Amplitude (a.u.)','FontSize',10);
    title('Before correction','FontSize',12);
    box off;
    subplot(1,2,2);
    plot(out_w_raw.ppm,spec_w_pre(:,:,1,1));xlim([4 5]);
    set(gca,'FontSize',12);
    set(gca,'XDir','reverse');
    xlabel('Frequency (ppm)','FontSize',10);
    ylabel('Amplitude(a.u.)','FontSize',10);
    title('After correction','FontSize',12);
    box off;
    set(h,'PaperUnits','centimeters');
    set(h,'PaperPosition',[0 0 20 10]);
    saveas(h,fullfile(reportFigsDir,'coilReconFig_w'),'jpg');
    saveas(h,fullfile(reportFigsDir,'coilReconFig_w'),'fig');
    close(h);
end

%%%%%%%%%%%%%%%%%%%%%OPTIONAL ALIGNMENT OF SUBSPECTRA%%%%%%%%%%%%%%%%
close all;
fs_ai=[];
phs_ai=[];
alignISIS='y'
if strcmp(alignISIS,'y') || strcmp(alignISIS,'Y')
    %What we're actually doing is aligning the averages, then aligning the
    %subspectra, then aligning the averages again, and then aligning the
    %subspectra again.  
    [out_ai,fs_temp,phs_temp]=op_alignAverages(out_cc,0.4,'y');
    fs_ai=fs_temp;
    phs_ai=phs_temp;
    [out_ai,fs_temp,phs_temp]=op_alignISIS(out_ai,0.4);
    fs_ai(:,2)=fs_ai(:,2)+fs_temp;
    phs_ai(:,2)=phs_ai(:,2)+phs_temp;
    [out_ai,fs_temp,phs_temp]=op_alignAverages(out_ai,0.4,'y');
    fs_ai=fs_ai+fs_temp;
    phs_ai=phs_ai+phs_temp;
    [out_ai,fs_temp,phs_temp]=op_alignISIS(out_ai,0.4);
    fs_ai(:,2)=fs_ai(:,2)+fs_temp;
    phs_ai(:,2)=phs_ai(:,2)+phs_temp;
    
    %for fs_ai and phs_ai, take the average across both subspecs:
    fs_ai=mean(fs_ai,2);
    phs_ai=mean(phs_ai,2);
    
    if water
        %Now repeat above for water unsuppressed data:
        [out_w_ai,fs_w_temp,phs_w_temp]=op_alignAverages(out_w_cc,0.4,'y');
        fs_w_ai=fs_w_temp;
        phs_w_ai=phs_w_temp;
        [out_w_ai,fs_w_temp,phs_w_temp]=op_alignISIS(out_w_ai,0.4);
        fs_w_ai(:,2)=fs_w_ai(:,2)+fs_w_temp;
        phs_w_ai(:,2)=phs_w_ai(:,2)+phs_w_temp;
        [out_w_ai,fs_w_temp,phs_w_temp]=op_alignAverages(out_w_ai,0.4,'y');
        fs_w_ai=fs_w_ai+fs_w_temp;
        phs_w_ai=phs_w_ai+phs_w_temp;
        [out_w_ai,fs_w_temp,phs_w_temp]=op_alignISIS(out_w_ai,0.4);
        fs_w_ai(:,2)=fs_w_ai(:,2)+fs_w_temp;
        phs_w_ai(:,2)=phs_w_ai(:,2)+phs_w_temp;
        
        %for fs_w_ai and phs_w_ai, take the average across both subspecs:
        fs_w_ai=mean(fs_w_ai,2);
        phs_w_ai=mean(phs_w_ai,2);
    end
    
    %Now check the plots to make sure that they look okay:
    %First make combined subspecs plots:
    out_cc_temp=op_combinesubspecs(out_cc,'diff');
    out_ai_cc_temp=op_combinesubspecs(out_ai,'diff');
    if water
        out_w_cc_temp=op_combinesubspecs(out_w_cc,'diff');
        out_w_ai_cc_temp=op_combinesubspecs(out_w_ai,'diff');
    end
    
    %Now plot them
    close all
    %figure('position',[0 0 560 420]);
    h=figure('visible','off');
    subplot(1,2,1);
    plot(out_cc_temp.ppm,out_cc_temp.specs);
    set(gca,'FontSize',12);
    set(gca,'XDir','reverse');
    xlabel('Frequency (ppm)','FontSize',10);
    ylabel('Amplitude(a.u.)','FontSize',10);
    title('Subspecs not aligned: (all averages)','FontSize',12);
    xlim([0.2 5.2]);
    box off;
    subplot(1,2,2);
    plot(out_ai_cc_temp.ppm,out_ai_cc_temp.specs);
    set(gca,'FontSize',8);
    set(gca,'XDir','reverse');
    xlabel('Frequency (ppm)','FontSize',10);
    ylabel('Amplitude(a.u.)','FontSize',10);
    title('Subspecs aligned: (all averages)','FontSize',12);
    xlim([0.2 5.2]);
    box off;
    set(h,'PaperUnits','centimeters');
    set(h,'PaperPosition',[0 0 20 15]);
    saveas(h,[filestring '/report/figs/alignSubSpecs_prePostFig'],'jpg');
    saveas(h,[filestring '/report/figs/alignSubSpecs_prePostFig'],'fig');
    close(h);
    
    
    if water
        %figure('position',[0 550 560 420]);
        h=figure('visible','off');
        subplot(1,2,1);
        plot(out_w_cc_temp.ppm,out_w_cc_temp.specs);
        xlim([3.7 5.7]);
        set(gca,'FontSize',12);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Subspecs not aligned: (all averages)','FontSize',12);
        box off;
        subplot(1,2,2);
        plot(out_w_ai_cc_temp.ppm,out_w_ai_cc_temp.specs);
        xlim([3.7 5.7]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Subspecs aligned: (all averages)','FontSize',12);
        box off;
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 20 15]);
        saveas(h,fullfile(reportFigsDir,'alignSubSpecs_w_prePostFig'),'jpg');
        saveas(h,fullfile(reportFigsDir,'alignSubSpecs_w_prePostFig'),'fig');
        close(h);
    end
    
%     figure('position',[570 50 560 420]);
%     subplot(2,1,1);
%     plot([1:length(fs_ai)],fs_ai);
%     xlabel('Scan Number');
%     ylabel('Frequency Drift (Hz)');
%     title('Estimated Frequency Drift');
%     subplot(2,1,2);
%     plot([1:length(phs_ai)],phs_ai);
%     xlabel('Scan Number');
%     ylabel('Phase Drift (Degrees)');
%     title('Estimated Phase Drift');
    
    
    %sat=input('are you satisfied with alignment of subspecs? ','s');
    sat='y';
    if strcmp(sat,'n') || strcmp(sat,'N')
        out_ai=out_cc;
        if water
            out_w_ai=out_w_cc;
        end
    end
    
else
    out_ai=out_cc;
    out_w_ai=out_w_cc;
    fs_ai=zeros(size(out_cc.fids,out_cc.dims.averages),1);
    phs_ai=zeros(size(out_cc.fids,out_cc.dims.averages),1);
end
    

%Now combine the subspecs
out_cs=op_combinesubspecs(out_ai,'diff');
if water
    out_w_cs=op_combinesubspecs(out_w_ai,'diff');
end

%%%%%%%%%%%%%%%%%%%%%END OF ALIGNMENT OF SUBSPECTRA%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%OPTIONAL REMOVAL OF BAD AVERAGES%%%%%%%%%%%%%%%%
close all;
out_cs2=out_cs;
nBadAvgTotal=0;
nbadAverages=1;
allAveragesLeft=[1:out_cs.sz(out_cs.dims.averages)]';
allBadAverages=[];
rmbadav='y';
close all;
if rmbadav=='n' || rmbadav=='N'
    out_rm=out_cs;
    nsd='N/A';
else
    sat='n'
    while sat=='n'||sat=='N'
        nsd=3; %setting the number of standard deviations;
        iter=1;
        nbadAverages=1;
        nBadAvgTotal=0;
        allAveragesLeft=[1:out_cs.sz(out_cs.dims.averages)]';
        allBadAverages=[];
        out_cs2=out_cs;
        while nbadAverages>0;
            [out_rm,metric{iter},badAverages]=op_rmbadaverages(out_cs2,nsd,'t');
            badAverages;
            allBadAverages=[allBadAverages; allAveragesLeft(badAverages)];
            badavMask_temp=zeros(length(allAveragesLeft),1);
            badavMask_temp(badAverages)=1;
            allAveragesLeft=allAveragesLeft(~badavMask_temp);
            nbadAverages=numel(badAverages);
            nBadAvgTotal=nBadAvgTotal+nbadAverages;
            out_cs2=out_rm;
            iter=iter+1;
            disp([num2str(nbadAverages) ' bad averages removed on this iteration.']);
            disp([num2str(nBadAvgTotal) ' bad averages removed in total.']);
            close all;
        end
        
        %figure('position',[0 50 560 420]);
        h=figure('visible','off');
        subplot(1,2,1);
        plot(out_cs.ppm,out_cs.specs);xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Before removal of bad averages:','FontSize',12);
        box off;
        subplot(1,2,2);
        plot(out_rm.ppm,out_rm.specs);xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('After removal of bad averages:','FontSize',12);
        box off;
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 20 15]);
        saveas(h,fullfile(reportFigsDir,'rmBadAvg_prePostFig'),'jpg');
        saveas(h,fullfile(reportFigsDir,'rmBadAvg_prePostFig'),'fig');
        close(h);
        
        %figure('position',[0 550 560 420]);
        h=figure('visible','off');        
        plot(1:out_cs.sz(out_cs.dims.averages),metric{1},'.r','MarkerSize',16)
        hold on
        removedAvgs = find(~ismember(1:out_cs.sz(out_cs.dims.averages),allAveragesLeft));
        plot(removedAvgs,metric{1}(removedAvgs,1),'ko','MarkerSize',20)
        set(gca,'FontSize',8);
        xlabel('Scan Number','FontSize',10);
        ylabel('Deviation Metric','FontSize',10);
        legend('Original averages','Removed Avg')
        legend boxoff;
        title('Deviation Metric','FontSize',12);
        box off;
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 20 10]);
        saveas(h,fullfile(reportFigsDir,'rmBadAvg_scatterFig'),'jpg');
        saveas(h,fullfile(reportFigsDir,'rmBadAvg_scatterFig'),'fig');
        close(h);
        
        %sat=input('are you satisfied with bad averages removal? ','s');
        sat='y';
    end
end

close all;

%Now remove the entries from fs_ai and phs_ai that correspond to
%the bad-averages that were removed.
BadAvgMask=zeros(length(fs_ai),1);
BadAvgMask(allBadAverages)=1;
size(BadAvgMask)
size(fs_ai)
fs_ai=fs_ai(~BadAvgMask);
phs_ai=phs_ai(~BadAvgMask);

%%%%%%%%%%%%%%%%%%%%END OF BAD AVERAGES REMOVAL%%%%%%%%%%%%%%%%%%%%

%now align averages;
%driftCorr=input('Would you like to perform the frequency drift correction?  ','s');
driftCorr='y';
if driftCorr=='n'|| driftCorr=='N'
    out_aa=out_rm;
    if water
        out_w_aa=out_w_cs;
    end
else
    sat='n'
    while sat=='n' || sat=='N'
        out_rm2=out_rm;
        fsPoly=100;
        phsPoly=1000;
        fsCum=fs_ai;
        phsCum=phs_ai;
        iter=1;
        while (abs(fsPoly(1))>0.001 || abs(phsPoly(1))>0.01) && iter<iterin
            close all
            if aaDomain=='t' || aaDomain=='T'
                tmax=input('input tmax for drift correction: ');
                [out_aa,fs,phs]=op_alignAverages(out_rm2,tmax,'n');
            elseif aaDomain=='f' || aaDomain=='F'
                tmax=tmaxin+0.04*randn(1);
                fmin=1.8+0.1*randn(1);
                fmaxarray=[2.4,2.85,3.35,4.2,4.4,5.2];
                fmax=fmaxarray(randi(6,1))
                [out_aa,fs,phs]=op_alignAverages_fd(out_rm2,fmin,fmax,tmax,'n');
            end
            if water
                [out_w_aa,fs_w,phs_w]=op_alignAverages(out_w_cs,5*tmax,'n');
            end
            
            fsCum=fsCum+fs;
            phsCum=phsCum+phs;
            fsPoly=polyfit([1:out_aa.sz(out_aa.dims.averages)]',fs,1)
            phsPoly=polyfit([1:out_aa.sz(out_aa.dims.averages)]',phs,1)
            iter
            out_rm2=out_aa;
            if water;
                out_w_cs=out_w_aa;
            end
            iter=iter+1;
        end
        %Now plot the cumulative frequency drift correction:
        %figure('position',[0 550 1125 420]);
        h=figure('visible','off');
        subplot(2,1,1);
        plot(out_rm.ppm,out_rm.specs);xlim([0 6]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Before drift correction:');
        box off;
        subplot(2,1,2);
        plot(out_aa.ppm,out_aa.specs);xlim([0 6]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('After drift correction:');
        box off;
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 20 15]);
        saveas(h,fullfile(reportFigsDir,'alignAvgs_prePostFig'),'jpg');
        saveas(h,fullfile(reportFigsDir,'alignAvgs_prePostFig'),'fig');
        close(h);
        
        %figure('position',[0 50 560 420]);
        h=figure('visible','off');
        plot([1:out_aa.sz(out_aa.dims.averages)],phsCum);
        set(gca,'FontSize',8);
        xlabel('Scan Number','FontSize',10);
        ylabel('Phase Drift (deg.)','FontSize',10);
        title('Estimated Phase Drift','FontSize',12);
        box off;
        legend('Phase drift','Location','SouthEast');
        legend boxoff;
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 10 10]);
        saveas(h,fullfile(reportFigsDir,'phaseDriftFig'),'jpg');
        saveas(h,fullfile(reportFigsDir,'phaseDriftFig'),'fig');
        close(h);
        
        %figure('position',[570 50 560 420]);
        h=figure('visible','off');
        plot([1:out_aa.sz(out_aa.dims.averages)],fsCum);
        set(gca,'FontSize',8);
        xlabel('Scan Number','FontSize',10);
        ylabel('Frequency Drift [Hz]','FontSize',10);
        box off;
        legend('Frequency Drift','Location','SouthEast');
        legend boxoff;
        title('Estimated Frequency Drift','FontSize',12);
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 10 10]);
        saveas(h,fullfile(reportFigsDir,'freqDriftFig'),'jpg');
        saveas(h,fullfile(reportFigsDir,'freqDriftFig'),'fig');
        close(h);
        
        %sat=input('Are you satisfied with the drift correction? ','s');
        sat='y';
        totalFreqDrift=mean(max(fsCum)-min(fsCum));
        totalPhaseDrift=mean(max(phsCum)-min(phsCum));
    end
    
end



%now do the averaging and left shift to get rid of first order phase:
out_av=op_leftshift(op_averaging(out_aa),out_aa.pointsToLeftshift);
if water
    out_w_av=op_leftshift(op_averaging(out_w_aa),out_w_aa.pointsToLeftshift);
end

%Now do an automatic phase correction (Use Creatine peak):
out_av_zp=op_zeropad(out_av,16);
index=find(abs(out_av_zp.specs)==max(abs(out_av_zp.specs(out_av_zp.ppm>2.85 & out_av_zp.ppm<3.15,1))));
ph0=-phase(out_av_zp.specs(index,1))*180/pi;
out_ph=op_addphase(out_av,ph0);
out_noproc=op_addphase(op_leftshift(out_noproc,out_noproc.pointsToLeftshift),ph0);
%And now for water unsuppressed data (use water peak):
if water
    out_w_av_zp=op_zeropad(out_w_av,16);
    indexw=find(abs(out_w_av_zp.specs)==max(abs(out_w_av_zp.specs(out_w_av_zp.ppm>4 & out_w_av_zp.ppm<5.5))));
    ph0w=-phase(out_w_av_zp.specs(indexw))*180/pi;
    out_w_ph=op_addphase(out_w_av,ph0w);
    out_w_noproc=op_addphase(op_leftshift(out_w_noproc,out_w_noproc.pointsToLeftshift),ph0w);
    out_w_ph_zp=op_addphase(out_w_av_zp,ph0w);
end

%Frequency shift all spectra so that Creatine appears at 3.027 ppm:
[~,frqShift]=op_ppmref(out_av_zp,2.9,3.1,3.027);
out=op_freqshift(out_ph,frqShift);
out_noproc=op_freqshift(out_noproc,frqShift);
%And now for water unsuppressed data (user water peak and set to 4.65 ppm):
if water
    [~,frqShiftw]=op_ppmref(out_w_ph_zp,4,5.5,4.65);
    out_w=op_freqshift(out_w_ph,frqShiftw);
    out_w_noproc=op_freqshift(out_w_noproc,frqShiftw);
end

%Make figure to show the final spectrum:
h=figure('visible','off');
plot(out.ppm,real(out.specs),'linewidth',2);xlim([0.2 5.2]);
set(gca,'FontSize',8);
set(gca,'XDir','reverse');
xlabel('Frequency (ppm)','FontSize',10);
ylabel('Amplitude(a.u.)','FontSize',10);
box off;
title('Result: Final Spectrum','FontSize',12);
set(h,'PaperUnits','centimeters');
set(h,'PaperPosition',[0 0 20 10]);
saveas(h,fullfile(reportFigsDir,'finalSpecFig'),'jpg');
saveas(h,fullfile(reportFigsDir,'finalSpecFig'),'fig');


% wrt=input('write? ','s');
wrt='y';
if wrt=='y' || wrt=='Y'
    RF=io_writelcm(out,fullfile(filestring,'main_lcm'),out.te);
    RF=io_writelcm(out_noproc,fullfile(filestring,'unprocessed_lcm'),out_noproc.te);
    if water
        RF=io_writelcm(out_w,fullfile([filestring '_w'],'w_lcm'),out_w.te);
        RF=io_writelcm(out_w_noproc,fullfile([filestring '_w'],'w_unprocessed_lcm'),out_w_noproc.te);
    end
end

close all

%write an html report: 
fid=fopen(fullfile(reportDir,'report.html'),'w+');
fprintf(fid,'<!DOCTYPE html>');
fprintf(fid,'\n<html>');
logoPath=which('FID-A_LOGO.jpg');
fprintf(fid,'\n<img src= " %s " width="120" height="120"></body>',logoPath);
fprintf(fid,'\n<h1>FID-A Processing Report</h1>');
fprintf(fid,'\n<h2>Processing pipeline applied to SPECIAL data using run_specialproc_auto.m</h2>');
fprintf(fid,'\n<p>FILENAME: %s/%s/%s </p>',fullfile(filestring,filename));
fprintf(fid,'\n<p>DATE: %s </p>',date);
fprintf(fid,'\n\n<p> </p>');
fprintf(fid,'\n\n<h2>Results of multi-coil combination:</h2>');
fprintf(fid,'\n<img src= " %s " width="800" height="400"></body>',fullfile(reportFigsDir,'coilReconFig.jpg'));
fprintf(fid,'\n\n<p> </p>');
fprintf(fid,'\n\n<h2>Results of alignment of SPECIAL sub-spectra:</h2>');
fprintf(fid,'\n<img src= " %s " width="800" height="400"></body>',fullfile(reportFigsDir,'alignSubSpecs_prePostFig.jpg'));
fprintf(fid,'\n\n<p> </p>');
fprintf(fid,'\n\n<h2>Results of removal of bad averages:</h2>');
fprintf(fid,'\n<p>Original number of averages: \t%5.6f </p>',out_raw.sz(out_raw.dims.averages));
fprintf(fid,'\n<p>Number of bad Averages removed:  \t%5.6f </p>',nBadAvgTotal);
fprintf(fid,'\n<p>Number of remaining averages in processed dataset:  \t%5.6f </p>',out_rm.sz(out_rm.dims.averages));
fprintf(fid,'\n<p>Bad Averages Removal Threshold was:  \t%2.2f </p>',nsd);
fprintf(fid,'\n<img src= " %s " width="800" height="600"><img src= " %s " width="800" height="400">',fullfile(reportFigsDir,'rmBadAvg_prePostFig.jpg'),fullfile(reportFigsDir,'rmBadAvg_scatterFig.jpg'));
fprintf(fid,'\n\n<p> </p>');
fprintf(fid,'\n\n<h2>Results of spectral registration:</h2>');
fprintf(fid,'\n<p>Total frequency drift was: \t%5.6f </p>',max(totalFreqDrift));
fprintf(fid,'\n<p>Total phase drift was: \t%5.6f </p>',max(totalPhaseDrift));
fprintf(fid,'\n<img src= " %s " width="800" height="600">',fullfile(reportFigsDir,'alignAvgs_prePostFig.jpg'));
fprintf(fid,'\n\n<p> </p>');
fprintf(fid,'\n<img src= " %s " width="400" height="400"><img src=" %s " width="400" height="400">',fullfile(reportFigsDir,'freqDriftFig.jpg'),fullfile(reportFigsDir,'phaseDriftFig.jpg'));
fprintf(fid,'\n\n<p> </p>');
fprintf(fid,'\n\n<h2>Final Result:</h2>');
fprintf(fid,'\n<img src= " %s " width="800" height="400">',fullfile(reportFigsDir,'finalSpecFig.jpg'));
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








