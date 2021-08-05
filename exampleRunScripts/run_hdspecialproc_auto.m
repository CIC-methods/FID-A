% run_hd_specialproc_auto.m
% Jamie Near, Sunnybrook Research Institute 2021.
% 
% USAGE:
% [out1,out2,out1_w,out2_w]=run_hd_specialproc_auto(filestring,aaDomain,tmaxin,iterin);
% 
% DESCRIPTION:
% Processing script for Siemens Hadamard-encoded SPECIAL (HD-SPECIAL) MRS 
% data in .dat format (twix raw data).  Includes combination of reciever 
% channels, removal of bad averages, frequency drift correction, 
% leftshifting, and Hadamard Reconstruction.
% 
% INPUTS:
% filestring   = String variable for the name of the directory containing
%                   the water suppressed HD-SPECIAL .dat file.  Water unsuppressed
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
% out1         = Fully processed, water suppressed output spectrum from first voxel.
% out2         = Fully processed, water suppressed output spectrum from second voxel. 
% out1_w       = Fully processed, water unsuppressed output spectrum from first voxel.
% out2_w       = Fully processed, water unsuppressed output spectrum from second voxel. 


function [out1,out2,out1_w,out2_w]=run_hdspecialproc_auto(filestring,aaDomain,tmaxin,iterin)

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
    coilcombos1=op_getcoilcombos(op_combinesubspecs(op_fourStepCombine(out_w_raw,1),'diff'),2);
    coilcombos2=op_getcoilcombos(op_combinesubspecs(op_fourStepCombine(out_w_raw,3),'diff'),2);
else
    disp('***WATER UNSUPPRESSED DATA NOT FOUND***');
    coilcombos1=op_getcoilcombos_specReg(op_combinesubspecs(op_fourStepCombine(op_averaging(out_raw),1),'diff'),0,0.01,2);
    coilcombos2=op_getcoilcombos_specReg(op_combinesubspecs(op_fourStepCombine(op_averaging(out_raw),3),'diff'),0,0.01,2);
    out_w='';
    out_w_noproc='';
end

%now combine the coil channels:
[out1_cc,fid1_pre,spec1_pre,ph1,sig1]=op_addrcvrs(out_raw,2,'w',coilcombos1);
[out2_cc,fid2_pre,spec2_pre,ph2,sig2]=op_addrcvrs(out_raw,2,'w',coilcombos2);
if water
    [out1_w_cc,fid1_w_pre,spec1_w_pre,ph1_w,sig1_w]=op_addrcvrs(out_w_raw,2,'w',coilcombos1);
    [out2_w_cc,fid2_w_pre,spec2_w_pre,ph2_w,sig2_w]=op_addrcvrs(out_w_raw,2,'w',coilcombos2);
end

%Do the hadamard reconstruction:
out1_hd=op_fourStepCombine(out1_cc,1);
out2_hd=op_fourStepCombine(out2_cc,3);
if water
    out1_w_hd=op_fourStepCombine(out1_w_cc,1);
    out2_w_hd=op_fourStepCombine(out2_w_cc,3);
end        

%make the un-processed spectra, which may be optionally output from this function:
out1_noproc=op_combinesubspecs(op_averaging(out1_hd),'diff');
out2_noproc=op_combinesubspecs(op_averaging(out2_hd),'diff');
if water
    out1_w_noproc=op_combinesubspecs(op_averaging(out1_w_hd),'diff');
    out2_w_noproc=op_combinesubspecs(op_averaging(out2_w_hd),'diff');
end

%plot the data before and after coil phasing:
%figure('position',[0 50 560 420]);
h=figure('visible','off');
subplot(1,3,1);
plot(out_raw.ppm,out_raw.specs(:,:,1,1));xlim([-1 7]);
set(gca,'FontSize',12);
set(gca,'XDir','reverse');
xlabel('Frequency (ppm)','FontSize',10);
ylabel('Amplitude (a.u.)','FontSize',10);
title('Before correction','FontSize',12);
box off;
subplot(1,3,2);
plot(out_raw.ppm,spec1_pre(:,:,1,1));xlim([-1 7]);
set(gca,'FontSize',12);
set(gca,'XDir','reverse');
xlabel('Frequency (ppm)','FontSize',10);
ylabel('Amplitude(a.u.)','FontSize',10);
title('After correction','FontSize',12);
box off;
subplot(1,3,3);
plot(out_raw.ppm,spec2_pre(:,:,1,1));xlim([-1 7]);
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
    subplot(1,3,1);
    plot(out_w_raw.ppm,out_w_raw.specs(:,:,1,1));xlim([4 5]);
    set(gca,'FontSize',12);
    set(gca,'XDir','reverse');
    xlabel('Frequency (ppm)','FontSize',10);
    ylabel('Amplitude (a.u.)','FontSize',10);
    title('Before correction','FontSize',12);
    box off;
    subplot(1,3,2);
    plot(out_w_raw.ppm,spec1_w_pre(:,:,1,1));xlim([4 5]);
    set(gca,'FontSize',12);
    set(gca,'XDir','reverse');
    xlabel('Frequency (ppm)','FontSize',10);
    ylabel('Amplitude(a.u.)','FontSize',10);
    title('After correction','FontSize',12);
    box off;
    subplot(1,3,3);
    plot(out_w_raw.ppm,spec2_w_pre(:,:,1,1));xlim([4 5]);
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
fs1_ai=[];
phs1_ai=[];
fs2_ai=[];
phs2_ai=[];
alignISIS='y'
if strcmp(alignISIS,'y') || strcmp(alignISIS,'Y')
    %What we're actually doing is aligning the averages, then aligning the
    %subspectra, then aligning the averages again, and then aligning the
    %subspectra again.  
    [out1_ai,fs1_temp,phs1_temp]=op_alignAverages(out1_hd,0.4,'y');
    [out2_ai,fs2_temp,phs2_temp]=op_alignAverages(out2_hd,0.4,'y');
    fs1_ai=fs1_temp;
    fs2_ai=fs2_temp;
    phs1_ai=phs1_temp;
    phs2_ai=phs2_temp;
    
    [out1_ai,fs1_temp,phs1_temp]=op_alignISIS(out1_ai,0.4);
    [out2_ai,fs2_temp,phs2_temp]=op_alignISIS(out2_ai,0.4);
    fs1_ai(:,2)=fs1_ai(:,2)+fs1_temp;
    fs2_ai(:,2)=fs2_ai(:,2)+fs2_temp;
    phs1_ai(:,2)=phs1_ai(:,2)+phs1_temp;
    phs2_ai(:,2)=phs2_ai(:,2)+phs2_temp;
    
    [out1_ai,fs1_temp,phs1_temp]=op_alignAverages(out1_ai,0.4,'y');
    [out2_ai,fs2_temp,phs2_temp]=op_alignAverages(out2_ai,0.4,'y');
    fs1_ai=fs1_ai+fs1_temp;
    fs2_ai=fs2_ai+fs2_temp;
    phs1_ai=phs1_ai+phs1_temp;
    phs2_ai=phs2_ai+phs2_temp;
    
    [out1_ai,fs1_temp,phs1_temp]=op_alignISIS(out1_ai,0.4);
    [out2_ai,fs2_temp,phs2_temp]=op_alignISIS(out2_ai,0.4);
    fs1_ai(:,2)=fs1_ai(:,2)+fs1_temp;
    fs2_ai(:,2)=fs2_ai(:,2)+fs2_temp;
    phs1_ai(:,2)=phs1_ai(:,2)+phs1_temp;
    phs2_ai(:,2)=phs2_ai(:,2)+phs2_temp;
    
    %for fs_ai and phs_ai, take the average across both subspecs:
    fs1_ai=mean(fs1_ai,2);
    fs2_ai=mean(fs2_ai,2);
    phs1_ai=mean(phs1_ai,2);
    phs2_ai=mean(phs2_ai,2);
    
    if water
        %Now repeat above for water unsuppressed data:
        [out1_w_ai,fs1_w_temp,phs1_w_temp]=op_alignAverages(out1_w_hd,0.4,'y');
        [out2_w_ai,fs2_w_temp,phs2_w_temp]=op_alignAverages(out2_w_hd,0.4,'y');
        fs1_w_ai=fs1_w_temp;
        fs2_w_ai=fs2_w_temp;
        phs1_w_ai=phs1_w_temp;
        phs2_w_ai=phs2_w_temp;
        
        [out1_w_ai,fs1_w_temp,phs1_w_temp]=op_alignISIS(out1_w_ai,0.4);
        [out2_w_ai,fs2_w_temp,phs2_w_temp]=op_alignISIS(out2_w_ai,0.4);
        fs1_w_ai(:,2)=fs1_w_ai(:,2)+fs1_w_temp;
        fs2_w_ai(:,2)=fs2_w_ai(:,2)+fs2_w_temp;
        phs1_w_ai(:,2)=phs1_w_ai(:,2)+phs1_w_temp;
        phs2_w_ai(:,2)=phs2_w_ai(:,2)+phs2_w_temp;
        
        [out1_w_ai,fs1_w_temp,phs1_w_temp]=op_alignAverages(out1_w_ai,0.4,'y');
        [out2_w_ai,fs2_w_temp,phs2_w_temp]=op_alignAverages(out2_w_ai,0.4,'y');
        fs1_w_ai=fs1_w_ai+fs1_w_temp;
        fs2_w_ai=fs2_w_ai+fs2_w_temp;
        phs1_w_ai=phs1_w_ai+phs1_w_temp;
        phs2_w_ai=phs2_w_ai+phs2_w_temp;
        
        [out1_w_ai,fs1_w_temp,phs1_w_temp]=op_alignISIS(out1_w_ai,0.4);
        [out2_w_ai,fs2_w_temp,phs2_w_temp]=op_alignISIS(out2_w_ai,0.4);
        fs1_w_ai(:,2)=fs1_w_ai(:,2)+fs1_w_temp;
        fs2_w_ai(:,2)=fs2_w_ai(:,2)+fs2_w_temp;
        phs1_w_ai(:,2)=phs1_w_ai(:,2)+phs1_w_temp;
        phs2_w_ai(:,2)=phs2_w_ai(:,2)+phs2_w_temp;
        
        %for fs_w_ai and phs_w_ai, take the average across both subspecs:
        fs1_w_ai=mean(fs1_w_ai,2);
        fs2_w_ai=mean(fs2_w_ai,2);
        phs1_w_ai=mean(phs1_w_ai,2);
        phs2_w_ai=mean(phs2_w_ai,2);
    end
    
    %Now check the plots to make sure that they look okay:
    %First make combined subspecs plots:
    out1_hd_temp=op_combinesubspecs(out1_hd,'diff');
    out2_hd_temp=op_combinesubspecs(out2_hd,'diff');
    out1_ai_cc_temp=op_combinesubspecs(out1_ai,'diff');
    out2_ai_cc_temp=op_combinesubspecs(out2_ai,'diff');
    if water
        out1_w_hd_temp=op_combinesubspecs(out1_w_hd,'diff');
        out2_w_hd_temp=op_combinesubspecs(out2_w_hd,'diff');
        out1_w_ai_cc_temp=op_combinesubspecs(out1_w_ai,'diff');
        out2_w_ai_cc_temp=op_combinesubspecs(out2_w_ai,'diff');
    end
    
    %Now plot them
    close all
    %figure('position',[0 0 560 420]);
    h=figure('visible','off');
    subplot(2,2,1);
    plot(out1_hd_temp.ppm,out1_hd_temp.specs);
    set(gca,'FontSize',12);
    set(gca,'XDir','reverse');
    xlabel('Frequency (ppm)','FontSize',10);
    ylabel('Amplitude(a.u.)','FontSize',10);
    title('Subspecs not aligned: (all averages)','FontSize',12);
    xlim([0.2 5.2]);
    box off;
    subplot(2,2,3);
    plot(out2_hd_temp.ppm,out2_hd_temp.specs);
    set(gca,'FontSize',12);
    set(gca,'XDir','reverse');
    xlabel('Frequency (ppm)','FontSize',10);
    ylabel('Amplitude(a.u.)','FontSize',10);
    title('Subspecs not aligned: (all averages)','FontSize',12);
    xlim([0.2 5.2]);
    box off;
    subplot(2,2,2);
    plot(out1_ai_cc_temp.ppm,out1_ai_cc_temp.specs);
    set(gca,'FontSize',8);
    set(gca,'XDir','reverse');
    xlabel('Frequency (ppm)','FontSize',10);
    ylabel('Amplitude(a.u.)','FontSize',10);
    title('Subspecs aligned: (all averages)','FontSize',12);
    xlim([0.2 5.2]);
    box off;
    subplot(2,2,4);
    plot(out2_ai_cc_temp.ppm,out2_ai_cc_temp.specs);
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
        subplot(2,2,1);
        plot(out1_w_hd_temp.ppm,out1_w_hd_temp.specs);
        xlim([3.7 5.7]);
        set(gca,'FontSize',12);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Subspecs not aligned: (all averages)','FontSize',12);
        box off;
        subplot(2,2,3);
        plot(out2_w_hd_temp.ppm,out2_w_hd_temp.specs);
        xlim([3.7 5.7]);
        set(gca,'FontSize',12);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Subspecs not aligned: (all averages)','FontSize',12);
        box off;
        subplot(1,2,2);
        plot(out1_w_ai_cc_temp.ppm,out1_w_ai_cc_temp.specs);
        xlim([3.7 5.7]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Subspecs aligned: (all averages)','FontSize',12);
        box off;
        subplot(2,2,4);
        plot(out2_w_ai_cc_temp.ppm,out2_w_ai_cc_temp.specs);
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
        out1_ai=out1_hd;
        out2_ai=out2_hd;
        if water
            out1_w_ai=out1_w_hd;
            out2_w_ai=out2_w_hd;
        end
    end
    
else
    out1_ai=out1_hd;
    out2_ai=out2_hd;
    out1_w_ai=out1_w_hd;
    out2_w_ai=out2_w_hd;
    fs1_ai=zeros(size(out1_hd.fids,out1_hd.dims.averages),1);
    fs2_ai=zeros(size(out2_hd.fids,out2_hd.dims.averages),1);
    phs1_ai=zeros(size(out1_hd.fids,out1_hd.dims.averages),1);
    phs2_ai=zeros(size(out2_hd.fids,out2_hd.dims.averages),1);
end
    

%Now combine the subspecs
out1_cs=op_combinesubspecs(out1_ai,'diff');
out2_cs=op_combinesubspecs(out2_ai,'diff');
if water
    out1_w_cs=op_combinesubspecs(out1_w_ai,'diff');
    out2_w_cs=op_combinesubspecs(out2_w_ai,'diff');
end

%%%%%%%%%%%%%%%%%%%%%END OF ALIGNMENT OF SUBSPECTRA%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%OPTIONAL REMOVAL OF BAD AVERAGES%%%%%%%%%%%%%%%%
close all;
out1_cs2=out1_cs;
out2_cs2=out2_cs;
nBadAvgTotal1=0;
nBadAvgTotal2=0;
nbadAverages1=1;
nbadAverages2=1;
allAveragesLeft1=[1:out1_cs.sz(out1_cs.dims.averages)]';
allAveragesLeft2=[1:out2_cs.sz(out2_cs.dims.averages)]';
allBadAverages1=[];
allBadAverages2=[];
rmbadav='y';
close all;
if rmbadav=='n' || rmbadav=='N'
    out1_rm=out1_cs;
    out2_rm=out2_cs;
    nsd='N/A';
else
    sat='n'
    while sat=='n'||sat=='N'
        nsd=3; %setting the number of standard deviations;
        iter1=1;
        iter2=1;
        nbadAverages1=1;
        nbadAverages2=1;
        nBadAvgTotal1=0;
        nBadAvgTotal2=0;
        allAveragesLeft1=[1:out1_cs.sz(out1_cs.dims.averages)]';
        allAveragesLeft2=[1:out2_cs.sz(out2_cs.dims.averages)]';
        allBadAverages1=[];
        allBadAverages2=[];
        out1_cs2=out1_cs;
        out2_cs2=out2_cs;
        while nbadAverages1>0;
            [out1_rm,metric1{iter1},badAverages1]=op_rmbadaverages(out1_cs2,nsd,'t');
            allBadAverages1=[allBadAverages1; allAveragesLeft1(badAverages1)];
            badavMask1_temp=zeros(length(allAveragesLeft1),1);
            badavMask1_temp(badAverages1)=1;
            allAveragesLeft1=allAveragesLeft1(~badavMask1_temp);
            nbadAverages1=numel(badAverages1);
            nBadAvgTotal1=nBadAvgTotal1+nbadAverages1;
            out1_cs2=out1_rm;
            iter1=iter1+1;
            disp([num2str(nbadAverages1) ' bad averages removed from spectrum #1 on this iteration.']);
            disp([num2str(nBadAvgTotal1) ' bad averages removed from spectrum #1 in total.']);
            close all;
        end
        while nbadAverages2>0;
            [out2_rm,metric2{iter2},badAverages2]=op_rmbadaverages(out2_cs2,nsd,'t');
            allBadAverages2=[allBadAverages2; allAveragesLeft2(badAverages2)];
            badavMask2_temp=zeros(length(allAveragesLeft2),1);
            badavMask2_temp(badAverages2)=1;
            allAveragesLeft2=allAveragesLeft2(~badavMask2_temp);
            nbadAverages2=numel(badAverages2);
            nBadAvgTotal2=nBadAvgTotal2+nbadAverages2;
            out2_cs2=out2_rm;
            iter2=iter2+1;
            disp([num2str(nbadAverages2) ' bad averages removed from spectrum #2 on this iteration.']);
            disp([num2str(nBadAvgTotal2) ' bad averages removed from spectrum #2 in total.']);
            close all;
        end
        
        
        %figure('position',[0 50 560 420]);
        h=figure('visible','off');
        subplot(2,2,1);
        plot(out1_cs.ppm,out1_cs.specs);xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Before removal of bad averages:','FontSize',12);
        box off;
        subplot(2,2,3);
        plot(out2_cs.ppm,out2_cs.specs);xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Before removal of bad averages:','FontSize',12);
        box off;
        subplot(2,2,2);
        plot(out1_rm.ppm,out1_rm.specs);xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('After removal of bad averages:','FontSize',12);
        box off;
        subplot(2,2,4);
        plot(out2_rm.ppm,out2_rm.specs);xlim([1 5]);
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
        plot(1:out1_cs.sz(out1_cs.dims.averages),metric1{1},'.r','MarkerSize',16)
        hold on
        removedAvgs1 = find(~ismember(1:out1_cs.sz(out1_cs.dims.averages),allAveragesLeft1));
        plot(removedAvgs1,metric1{1}(removedAvgs1,1),'ko','MarkerSize',20)
        set(gca,'FontSize',8);
        xlabel('Scan Number','FontSize',10);
        ylabel('Deviation Metric','FontSize',10);
        legend('Original averages','Removed Avg')
        legend boxoff;
        title('Deviation Metric','FontSize',12);
        box off;
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 20 10]);
        saveas(h,fullfile(reportFigsDir,'rmBadAvg_vox1_scatterFig1'),'jpg');
        saveas(h,fullfile(reportFigsDir,'rmBadAvg_vox1_scatterFig1'),'fig');
        close(h);
        
        h=figure('visible','off');        
        plot(1:out2_cs.sz(out2_cs.dims.averages),metric2{1},'.r','MarkerSize',16)
        hold on
        removedAvgs2 = find(~ismember(1:out2_cs.sz(out2_cs.dims.averages),allAveragesLeft2));
        plot(removedAvgs2,metric2{1}(removedAvgs2,1),'ko','MarkerSize',20)
        set(gca,'FontSize',8);
        xlabel('Scan Number','FontSize',10);
        ylabel('Deviation Metric','FontSize',10);
        legend('Original averages','Removed Avg')
        legend boxoff;
        title('Deviation Metric','FontSize',12);
        box off;
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 20 10]);
        saveas(h,fullfile(reportFigsDir,'rmBadAvg_vox2_scatterFig'),'jpg');
        saveas(h,fullfile(reportFigsDir,'rmBadAvg_vox2_scatterFig'),'fig');
        close(h);
        
        %sat=input('are you satisfied with bad averages removal? ','s');
        sat='y';
    end
end

close all;

%Now remove the entries from fs_ai and phs_ai that correspond to
%the bad-averages that were removed.
BadAvgMask1=zeros(length(fs1_ai),1);
BadAvgMask2=zeros(length(fs2_ai),1);
BadAvgMask1(allBadAverages1)=1;
BadAvgMask2(allBadAverages2)=1;
size(BadAvgMask1)
size(BadAvgMask2)
size(fs1_ai)
size(fs2_ai)
fs1_ai=fs1_ai(~BadAvgMask1);
fs2_ai=fs2_ai(~BadAvgMask2);
phs1_ai=phs1_ai(~BadAvgMask1);
phs2_ai=phs2_ai(~BadAvgMask2);

%%%%%%%%%%%%%%%%%%%%END OF BAD AVERAGES REMOVAL%%%%%%%%%%%%%%%%%%%%

%now align averages;
%driftCorr=input('Would you like to perform the frequency drift correction?  ','s');
driftCorr='y';
if driftCorr=='n'|| driftCorr=='N'
    out1_aa=out1_rm;
    out2_aa=out2_rm;
    if water
        out1_w_aa=out1_w_cs;
        out2_w_aa=out2_w_cs;
    end
else
    sat='n'
    while sat=='n' || sat=='N'
        out1_rm2=out1_rm;
        out2_rm2=out2_rm;
        fsPoly1=100;
        fsPoly2=100;
        phsPoly1=1000;
        phsPoly2=1000;
        fsCum1=fs1_ai;
        fsCum2=fs2_ai;
        phsCum1=phs1_ai;
        phsCum2=phs2_ai;
        iter1=1;
        iter2=1;
        while (abs(fsPoly1(1))>0.001 || abs(phsPoly1(1))>0.01) && iter1<iterin
            close all
            if aaDomain=='t' || aaDomain=='T'
                tmax=input('input tmax for drift correction: ');
                [out1_aa,fs1,phs1]=op_alignAverages(out1_rm2,tmax,'n');
            elseif aaDomain=='f' || aaDomain=='F'
                tmax=tmaxin+0.04*randn(1);
                fmin=1.8+0.1*randn(1);
                fmaxarray=[2.4,2.85,3.35,4.2,4.4,5.2];
                fmax=fmaxarray(randi(6,1))
                [out1_aa,fs1,phs1]=op_alignAverages_fd(out1_rm2,fmin,fmax,tmax,'n');
            end
            if water
                [out1_w_aa,fs1_w,phs1_w]=op_alignAverages(out1_w_cs,5*tmax,'n');
            end
            
            fsCum1=fsCum1+fs1;
            phsCum1=phsCum1+phs1;
            fsPoly1=polyfit([1:out1_aa.sz(out1_aa.dims.averages)]',fs1,1)
            phsPoly1=polyfit([1:out1_aa.sz(out1_aa.dims.averages)]',phs1,1)
            iter1
            out1_rm2=out1_aa;
            if water;
                out1_w_cs=out1_w_aa;
            end
            iter1=iter1+1;
        end
        while (abs(fsPoly2(1))>0.001 || abs(phsPoly2(1))>0.01) && iter2<iterin
            close all
            if aaDomain=='t' || aaDomain=='T'
                tmax=input('input tmax for drift correction: ');
                [out2_aa,fs2,phs2]=op_alignAverages(out2_rm2,tmax,'n');
            elseif aaDomain=='f' || aaDomain=='F'
                tmax=tmaxin+0.04*randn(1);
                fmin=1.8+0.1*randn(1);
                fmaxarray=[2.4,2.85,3.35,4.2,4.4,5.2];
                fmax=fmaxarray(randi(6,1))
                [out2_aa,fs2,phs2]=op_alignAverages_fd(out2_rm2,fmin,fmax,tmax,'n');
            end
            if water
                [out2_w_aa,fs2_w,phs2_w]=op_alignAverages(out2_w_cs,5*tmax,'n');
            end
            
            fsCum2=fsCum2+fs2;
            phsCum2=phsCum2+phs2;
            fsPoly2=polyfit([1:out2_aa.sz(out2_aa.dims.averages)]',fs2,1)
            phsPoly2=polyfit([1:out2_aa.sz(out2_aa.dims.averages)]',phs2,1)
            iter2
            out2_rm2=out2_aa;
            if water;
                out2_w_cs=out2_w_aa;
            end
            iter2=iter2+1;
        end
        %Now plot the cumulative frequency drift correction:
        %figure('position',[0 550 1125 420]);
        h=figure('visible','off');
        subplot(2,2,1);
        plot(out1_rm.ppm,out1_rm.specs);xlim([0 6]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Before drift correction:');
        box off;
        subplot(2,2,2);
        plot(out2_rm.ppm,out2_rm.specs);xlim([0 6]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Before drift correction:');
        box off;
        subplot(2,2,3);
        plot(out1_aa.ppm,out1_aa.specs);xlim([0 6]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('After drift correction:');
        box off;
        subplot(2,2,4);
        plot(out2_aa.ppm,out2_aa.specs);xlim([0 6]);
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
        subplot(2,1,1);
        plot([1:out1_aa.sz(out1_aa.dims.averages)],phsCum1);
        set(gca,'FontSize',8);
        xlabel('Scan Number','FontSize',10);
        ylabel('Phase Drift (deg.)','FontSize',10);
        title('Estimated Phase Drift','FontSize',12);
        box off;
        subplot(2,1,2);
        plot([1:out2_aa.sz(out2_aa.dims.averages)],phsCum2);
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
        subplot(2,1,1);
        plot([1:out1_aa.sz(out1_aa.dims.averages)],fsCum1);
        set(gca,'FontSize',8);
        xlabel('Scan Number','FontSize',10);
        ylabel('Frequency Drift [Hz]','FontSize',10);
        box off;
        subplot(2,1,2);
        plot([1:out2_aa.sz(out2_aa.dims.averages)],fsCum2);
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
        totalFreqDrift1=mean(max(fsCum1)-min(fsCum1));
        totalFreqDrift2=mean(max(fsCum2)-min(fsCum2));
        totalPhaseDrift1=mean(max(phsCum1)-min(phsCum1));
        totalPhaseDrift2=mean(max(phsCum2)-min(phsCum2));
    end
    
end



%now do the averaging and left shift to get rid of first order phase:
out1_av=op_leftshift(op_averaging(out1_aa),out1_aa.pointsToLeftshift);
out2_av=op_leftshift(op_averaging(out2_aa),out2_aa.pointsToLeftshift);
if water
    out1_w_av=op_leftshift(op_averaging(out1_w_aa),out1_w_aa.pointsToLeftshift);
    out2_w_av=op_leftshift(op_averaging(out2_w_aa),out2_w_aa.pointsToLeftshift);
end

%Now do an automatic phase correction (Use Creatine peak):
out1_av_zp=op_zeropad(out1_av,16);
out2_av_zp=op_zeropad(out2_av,16);
index1=find(abs(out1_av_zp.specs)==max(abs(out1_av_zp.specs(out1_av_zp.ppm>2.85 & out1_av_zp.ppm<3.15,1))));
index2=find(abs(out2_av_zp.specs)==max(abs(out2_av_zp.specs(out2_av_zp.ppm>2.85 & out2_av_zp.ppm<3.15,1))));
ph01=-phase(out1_av_zp.specs(index1,1))*180/pi;
ph02=-phase(out2_av_zp.specs(index2,1))*180/pi;
out1_ph=op_addphase(out1_av,ph01);
out2_ph=op_addphase(out2_av,ph02);
out1_noproc=op_addphase(op_leftshift(out1_noproc,out1_noproc.pointsToLeftshift),ph01);
out2_noproc=op_addphase(op_leftshift(out2_noproc,out2_noproc.pointsToLeftshift),ph02);
%And now for water unsuppressed data (use water peak):
if water
    out1_w_av_zp=op_zeropad(out1_w_av,16);
    out2_w_av_zp=op_zeropad(out2_w_av,16);
    index1w=find(abs(out1_w_av_zp.specs)==max(abs(out1_w_av_zp.specs(out1_w_av_zp.ppm>4 & out1_w_av_zp.ppm<5.5))));
    index2w=find(abs(out2_w_av_zp.specs)==max(abs(out2_w_av_zp.specs(out2_w_av_zp.ppm>4 & out2_w_av_zp.ppm<5.5))));
    ph01w=-phase(out1_w_av_zp.specs(index1w))*180/pi;
    ph02w=-phase(out2_w_av_zp.specs(index2w))*180/pi;
    out1_w_ph=op_addphase(out1_w_av,ph01w);
    out2_w_ph=op_addphase(out2_w_av,ph02w);
    out1_w_noproc=op_addphase(op_leftshift(out1_w_noproc,out1_w_noproc.pointsToLeftshift),ph01w);
    out2_w_noproc=op_addphase(op_leftshift(out2_w_noproc,out2_w_noproc.pointsToLeftshift),ph02w);
    out1_w_ph_zp=op_addphase(out1_w_av_zp,ph01w);
    out2_w_ph_zp=op_addphase(out2_w_av_zp,ph02w);
end

%Frequency shift all spectra so that Creatine appears at 3.027 ppm:
[~,frqShift1]=op_ppmref(out1_av_zp,2.9,3.1,3.027);
[~,frqShift2]=op_ppmref(out2_av_zp,2.9,3.1,3.027);
out1=op_freqshift(out1_ph,frqShift1);
out2=op_freqshift(out2_ph,frqShift2);
out1_noproc=op_freqshift(out1_noproc,frqShift1);
out2_noproc=op_freqshift(out2_noproc,frqShift2);
%And now for water unsuppressed data (user water peak and set to 4.65 ppm):
if water
    [~,frqShift1w]=op_ppmref(out1_w_ph_zp,4,5.5,4.65);
    [~,frqShift2w]=op_ppmref(out2_w_ph_zp,4,5.5,4.65);
    out1_w=op_freqshift(out1_w_ph,frqShift1w);
    out2_w=op_freqshift(out2_w_ph,frqShift2w);
    out1_w_noproc=op_freqshift(out1_w_noproc,frqShift1w);
    out2_w_noproc=op_freqshift(out2_w_noproc,frqShift2w);
end

%Make figure to show the final spectrum:
h=figure('visible','off');
subplot(1,2,1);
plot(out1.ppm,real(out1.specs),'linewidth',2);xlim([0.2 5.2]);
set(gca,'FontSize',8);
set(gca,'XDir','reverse');
xlabel('Frequency (ppm)','FontSize',10);
ylabel('Amplitude(a.u.)','FontSize',10);
box off;
subplot(1,2,2);
plot(out2.ppm,real(out2.specs),'linewidth',2);xlim([0.2 5.2]);
set(gca,'FontSize',8);
set(gca,'XDir','reverse');
xlabel('Frequency (ppm)','FontSize',10);
ylabel('Amplitude(a.u.)','FontSize',10);
box off;
title('Result: Final Spectra','FontSize',12);
set(h,'PaperUnits','centimeters');
set(h,'PaperPosition',[0 0 20 10]);
saveas(h,fullfile(reportFigsDir,'finalSpecFig'),'jpg');
saveas(h,fullfile(reportFigsDir,'finalSpecFig'),'fig');


% wrt=input('write? ','s');
wrt='y';
if wrt=='y' || wrt=='Y'
    RF1=io_writelcm(out1,fullfile(filestring,'main1_lcm'),out1.te);
    RF2=io_writelcm(out2,fullfile(filestring,'main2_lcm'),out2.te);
    RF1=io_writelcm(out1_noproc,fullfile(filestring,'unprocessed1_lcm'),out1_noproc.te);
    RF2=io_writelcm(out2_noproc,fullfile(filestring,'unprocessed2_lcm'),out2_noproc.te);
    if water
        RF1=io_writelcm(out1_w,fullfile([filestring '_w'],'1w_lcm'),out1_w.te);
        RF2=io_writelcm(out2_w,fullfile([filestring '_w'],'2w_lcm'),out2_w.te);
        RF1=io_writelcm(out1_w_noproc,fullfile([filestring '_w'],'1w_unprocessed_lcm'),out1_w_noproc.te);
        RF2=io_writelcm(out2_w_noproc,fullfile([filestring '_w'],'2w_unprocessed_lcm'),out2_w_noproc.te);
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
fprintf(fid,'\n<h2>Processing pipeline applied to HD-SPECIAL data using run_hdspecialproc_auto.m</h2>');
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
fprintf(fid,'\n<p>Number of bad Averages removed from vox1:  \t%5.6f </p>',nBadAvgTotal1);
fprintf(fid,'\n<p>Number of bad Averages removed from vox2:  \t%5.6f </p>',nBadAvgTotal2);
fprintf(fid,'\n<p>Number of remaining averages in processed dataset (vox1):  \t%5.6f </p>',out1_rm.sz(out1_rm.dims.averages));
fprintf(fid,'\n<p>Number of remaining averages in processed dataset (vox2):  \t%5.6f </p>',out2_rm.sz(out2_rm.dims.averages));
fprintf(fid,'\n<p>Bad Averages Removal Threshold was:  \t%2.2f </p>',nsd);
fprintf(fid,'\n<img src= " %s " width="800" height="600"><img src= " %s " width="800" height="400"><img src= " %s " width="800" height="400">',fullfile(reportFigsDir,'rmBadAvg_prePostFig.jpg'),fullfile(reportFigsDir,'rmBadAvg_vox1_scatterFig.jpg'),fullfile(reportFigsDir,'rmBadAvg_vox2_scatterFig.jpg'));
fprintf(fid,'\n\n<p> </p>');
fprintf(fid,'\n\n<h2>Results of spectral registration:</h2>');
fprintf(fid,'\n<p>Total frequency drift from vox1 was: \t%5.6f </p>',max(totalFreqDrift1));
fprintf(fid,'\n<p>Total frequency drift from vox2 was: \t%5.6f </p>',max(totalFreqDrift2));
fprintf(fid,'\n<p>Total phase drift from vox1 was: \t%5.6f </p>',max(totalPhaseDrift1));
fprintf(fid,'\n<p>Total phase drift from vox2 was: \t%5.6f </p>',max(totalPhaseDrift2));
fprintf(fid,'\n<img src= " %s " width="800" height="600">',fullfile(reportFigsDir,'alignAvgs_prePostFig.jpg'));
fprintf(fid,'\n\n<p> </p>');
fprintf(fid,'\n<img src= " %s " width="400" height="400"><img src=" %s " width="400" height="400">',fullfile(reportFigsDir,'freqDriftFig.jpg'),fullfile(reportFigsDir,'phaseDriftFig.jpg'));
fprintf(fid,'\n\n<p> </p>');
fprintf(fid,'\n\n<h2>Final Result:</h2>');
fprintf(fid,'\n<img src= " %s " width="800" height="400">',fullfile(reportFigsDir,'finalSpecFig.jpg'));
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








