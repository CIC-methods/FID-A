% run_megapressproc_GEauto.m
% Jamie Near, McGill University 2017.
% 
% USAGE:
% [diffSpec,sumSpec,subSpec1,subSpec2,outw]=run_megapressproc_GEauto(filestring,coilcombos,avgAlignDomain,alignSS);
% 
% DESCRIPTION:
% Fully automated processing script for GE MEGA-PRESS MRS data in P-file 
% format (GE raw data).  Includes combination of reciever channels, 
% removal of bad averages, freqeuncy drift correction, manual alignment of 
% edit-on and edit-off spectra, and leftshifting.  This pipeline requires
% no user interaction. This function automatically 
% generates an html report to describe the results of each processing step.
% 
% INPUTS:
% filestring         = String variable for the name of the P file.
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
% diffSpec          = Fully processed difference spectrum.
% sumSpec           = Fully processed sum spectrum.
% subSpec1          = Fully processed MEGA-PRESS subspectrum #1.
% subSpec2          = Fully processed MEGA-PRESS subspectrum #2.
% outw              = Fully processed water unsuppressed spectrum. 
%

function [diffSpec,sumSpec,subSpec1,subSpec2,outw]=run_megapressproc_GEauto(filestring,coilcombos,avgAlignDomain,alignSS);

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
mkdir(['./report']);
mkdir(['./report/figs']);

% %read in both datasets:
%Assume spectral width = 5000 Hz;
%       txfrq          = 127.7 Hz/ppm;
%       nSubspecs      =2;
%       TE             =68 ms;
%       TR             =2000 ms;
[raw,raww]=io_loadspec_GE(filestring,2);
raw=op_complexConj(raw);

if isstruct(raww)
    water=true;
    raww=op_complexConj(raww);
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
        coilcombos=op_getcoilcombos(op_averaging(op_combinesubspecs(raw,'summ')),1);
    end
    
end
[out_cc,fid_pre,spec_pre,ph,sig]=op_addrcvrs(raw,1,'w',coilcombos);
[out_av_cc,fid_av_pre,spec_av_pre]=op_addrcvrs(op_averaging(raw),1,'w',coilcombos);
raw_av=op_averaging(raw);

%generate unprocessed spectrum:
out_noproc=op_combinesubspecs(op_averaging(out_cc),'diff');


%Generate plots showing coil channels before and after phase alignment
%figure('position',[0 50 560 420]);
h=figure('visible','off');
subplot(1,2,1);
plot(raw_av.ppm,real(raw_av.specs(:,:,1,1)));xlim([1 5]);
set(gca,'FontSize',12);
set(gca,'XDir','reverse');
xlabel('Frequency (ppm)','FontSize',10);
ylabel('Amplitude (a.u.)','FontSize',10);
title('Before correction','FontSize',12);
box off;
subplot(1,2,2);
plot(raw_av.ppm,real(spec_av_pre(:,:,1,1)));xlim([1 5]);
set(gca,'FontSize',12);
set(gca,'XDir','reverse');
xlabel('Frequency (ppm)','FontSize',10);
ylabel('Amplitude(a.u.)','FontSize',10);
title('After correction','FontSize',12);
box off;
set(h,'PaperUnits','centimeters');
set(h,'PaperPosition',[0 0 20 10]);
saveas(h,['./report/figs/coilReconFig'],'jpg');
saveas(h,['./report/figs/coilReconFig'],'fig');
close(h);


%%%%%%%%OPTIONAL REMOVAL OF BAD AVERAGES%%%%%%%%%%%%%%%%%%%%
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
            nbadAverages=length(badAverages)*raw.sz(raw.dims.subSpecs);
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
        subplot(2,2,1);
        plot(out_cc.ppm,real(out_cc.specs(:,:,1)));xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Edit-ON Before','FontSize',12);
        box off;
        subplot(2,2,2);
        plot(out_rm.ppm,real(out_rm.specs(:,:,1)));xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Edit-ON After','FontSize',12);
        box off;
        subplot(2,2,3);  
        plot(out_cc.ppm,real(out_cc.specs(:,:,2)));xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Edit-OFF Before','FontSize',12);
        box off;
        subplot(2,2,4);
        plot(out_rm.ppm,real(out_rm.specs(:,:,2)));xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Edit-OFF After','FontSize',12);
        box off;
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 20 15]);
        saveas(h,['./report/figs/rmBadAvg_prePostFig'],'jpg');
        saveas(h,['./report/figs/rmBadAvg_prePostFig'],'fig');
        close(h);
        
        %figure('position',[0 550 560 420]);
        h=figure('visible','off');
        plot([1:length(metric{1})],metric{1},'.r',[1:length(metric{iter-1})],metric{iter-1},'x','MarkerSize',16);
        set(gca,'FontSize',8);
        xlabel('Scan Number','FontSize',10);
        ylabel('Deviation Metric','FontSize',10);
        legend('Before rmBadAv','Before rmBadAv','After rmBadAv','After rmBadAv');
        legend boxoff;
        title('Deviation Metric','FontSize',12);
        box off;
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 20 10]);
        saveas(h,['./report/figs/rmBadAvg_scatterFig'],'jpg');
        saveas(h,['./report/figs/rmBadAvg_scatterFig'],'fig');
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
        iter=0;
        iterin=20;
        p=100;
        fscum=zeros(out_rm.sz(2:end));
        phscum=zeros(out_rm.sz(2:end));
        while (abs(p(1))>0.0003 && iter<iterin)
            iter=iter+1
            close all
            tmax=0.25+0.03*randn(1);
            ppmmin=1.6+0.1*randn(1);
            ppmmaxarray=[3.5+0.1*randn(1,2),4+0.1*randn(1,3),5.5+0.1*randn(1,1)];
            ppmmax=ppmmaxarray(randi(6,1));
            switch avgAlignDomain
                case 't'
                    [out_aa,fs,phs]=op_alignAverages(out_rm2,tmax,'y');
                case 'f'
                    [out_aa,fs,phs]=op_alignAverages_fd(out_rm2,ppmmin,ppmmax,tmax,'y');
                otherwise
                    error('ERROR: avgAlignDomain not recognized!');
            end
            
            x=repmat([1:size(fs,1)]',1,out_aa.sz(out_aa.dims.subSpecs));
            p=polyfit(x,fs,1)
            
            fscum=fscum+fs;
            phscum=phscum+phs;
            
            if driftCorr=='y' || driftCorr=='Y'
                out_rm2=out_aa;
            end
        end
        h=figure('visible','off');
        subplot(2,2,1);
        plot(out_rm.ppm,real(out_rm.specs(:,:,1)));xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Edit-ON Before','FontSize',12);
        box off;
        subplot(2,2,2);
        plot(out_aa.ppm,real(out_aa.specs(:,:,1)));xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Edit-ON After','FontSize',12);
        box off;
        subplot(2,2,3);
        plot(out_rm.ppm,real(out_rm.specs(:,:,2)));xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Edit-OFF Before','FontSize',12);
        box off;
        subplot(2,2,4);
        plot(out_aa.ppm,real(out_aa.specs(:,:,2)));xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Edit-OFF After','FontSize',12);
        box off;
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 20 15]);
        saveas(h,['./report/figs/alignAvgs_prePostFig'],'jpg');
        saveas(h,['./report/figs/alignAvgs_prePostFig'],'fig');
        close(h);
        
        h=figure('visible','off');
        plot([1:out_aa.sz(out_aa.dims.averages)],fscum,'.-','LineWidth',2);
        set(gca,'FontSize',8);
        xlabel('Scan Number','FontSize',10);
        ylabel('Frequency Drift [Hz]','FontSize',10);
        box off;
        legend('Edit-on scans','Edit-off scans','Location','SouthEast');
        legend boxoff;
        title('Estimated Freqeuncy Drift','FontSize',12);
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 10 10]);
        saveas(h,['./report/figs/freqDriftFig'],'jpg');
        saveas(h,['./report/figs/freqDriftFig'],'fig');
        close(h);
        
        h=figure('visible','off');
        plot([1:out_aa.sz(out_aa.dims.averages)],phscum,'.-','LineWidth',2);
        set(gca,'FontSize',8);
        xlabel('Scan Number','FontSize',10);
        ylabel('Phase Drift [Deg.]','FontSize',10);
        box off;
        legend('Edit-on scans','Edit-off scans');
        legend boxoff;
        title('Estimated Phase Drift','FontSize',12);
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 10 10]);
        saveas(h,['./report/figs/phaseDriftFig'],'jpg');
        saveas(h,['./report/figs/phaseDriftFig'],'fig');
        close(h);

        sat='y';
        if sat=='n'
            iter=0;
            p=100;
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

%now do automatic zero-order phase correction (Use water Peak):
out_ls_ss2=op_takesubspec(out_ls,2);
%SpecTool(out_ls_ss2,0.2,1,3.5);
[out_ls_ss2_ph,ph0]=op_autophase(out_ls_ss2,4,5.5);
%ph0=input('Input zero-order phase:');
out_ph=op_addphase(out_ls,ph0);
if water
    outw_ph=op_addphase(outw_ls,ph0);
end

%do same phase corection on unprocessed data
out_noproc=op_addphase(out_noproc,ph0);

%Now align subspecs if desired:
switch alignSS    
    
    case 2
        out=op_alignMPSubspecs(out_ph);
        
%         figure('position',[0 50 560 420]);
%         out_ph_filt=op_filter(out_ph,5);
%         subSpecTool(out_ph_filt,0,7);
%         disp('***************************************************************************************');
%         disp('Use GUI interface to align edit-ON and edit-OFF scans by adjusting Phase and Frequency.');
%         disp('Try to minimize the residual water, residual Creatine, and residual Choline peaks!');
%         disp('***NOTE If you are using the Siemens MEGA_PRESS WIP (WIP529), then you will');
%         disp('have to add about 180 degrees of phase to the subspectrum!***');
%         disp('*************************************************************');
%         fprintf('\n');
%         phshft1=input('Input Desired Phase Shift (Degrees) for first spectrum: ');
%         frqshft1=input('Input Desired Frequncy Shift (Hz) for first spectrum: ');
%         out=op_freqshiftSubspec(op_addphaseSubspec(out_ph,phshft1),frqshft1);
%         close all;
        
    case 0
        out=out_ph;        
  
    otherwise
        error('ERROR: alignSS value not valid! ');
end

%Make fully processed data;
diffSpec=op_combinesubspecs(out,'diff');
sumSpec=op_combinesubspecs(out,'summ');
subSpec1=op_takesubspec(out,1);
subSpec2=op_takesubspec(out,2);
subSpec2=op_addphase(subSpec2,180);

%Frequency shift all spectra so that Creatine appears at 3.027 ppm:
[subSpec1,frqShift]=op_ppmref(subSpec1,2.9,3.1,3.027);
diffSpec=op_freqshift(diffSpec,frqShift);
sumSpec=op_freqshift(sumSpec,frqShift);
subSpec2=op_freqshift(subSpec2,frqShift);


%Make final water unsuppressed data
if water
    if ~isempty(findstr(outw_ph.seq,'edit_529')) || ~isempty(findstr(outw_ph.seq,'jn_svs_special'))
        if outw_ph.dims.subSpecs
            outw=op_combinesubspecs(outw_ph,'diff');
        else
            outw=outw_ph;
        end
    else
        if outw_ph.dims.subSpecs
            outw=op_combinesubspecs(outw_ph,'summ');
        else
            outw=outw_ph;
        end
    end
    outw=op_addphase(outw,-phase(outw.fids(1))*180/pi,0,4.65,1);
else
    outw=0;
end

%Make figure to show the final spectrum:
h=figure('visible','off');
subplot(1,2,1);
plot(subSpec1.ppm,subSpec1.specs,subSpec2.ppm,subSpec2.specs,'linewidth',2);xlim([0.2 5.2]);
set(gca,'FontSize',8);
set(gca,'XDir','reverse');
xlabel('Frequency (ppm)','FontSize',10);
ylabel('Amplitude(a.u.)','FontSize',10);
legend('Edit-ON','Edit-OFF');
legend boxoff;
box off;
title('Result: Subspecs','FontSize',12);
subplot(1,2,2);
plot(diffSpec.ppm,diffSpec.specs,'linewidth',2);xlim([0.2 5.2]);
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
saveas(h,['./report/figs/finalSpecFig'],'jpg');
saveas(h,['./report/figs/finalSpecFig'],'fig');

% writ=input('Write? (y or n):  ','s');
writ='y';
if writ=='y' || writ=='Y'
     RF=io_writelcm(diffSpec,['./' filestring '_diff_lcm'],diffSpec.te);
     RF=io_writelcm(subSpec1,['./' filestring '_editON_lcm'],subSpec1.te);
     RF=io_writelcm(subSpec2,['./' filestring '_editOFF_lcm'],subSpec2.te);
    if water
        RF=io_writelcm(outw,['./' filestring '_w_lcm'],outw.te);
    end
end

close all;

%write an html report: 
fid=fopen(['./report/report.html'],'w+');
fprintf(fid,'<!DOCTYPE html>');
fprintf(fid,'\n<html>');
logoPath=which('FID-A_LOGO.jpg');
fprintf(fid,'\n<img src= " %s " width="120" height="120"></body>',logoPath);
fprintf(fid,'\n<h1>FID-A Processing Report</h1>');
fprintf(fid,'\n<h2>Processing pipeline applied to MEGA-PRESS data using run_megapressproc.m</h2>');
fprintf(fid,'\n<p>FILENAME: %s/%s </p>',pwd,filestring);
fprintf(fid,'\n<p>DATE: %s </p>',date);
fprintf(fid,'\n\n<p> </p>');
fprintf(fid,'\n\n<h2>Results of multi-coil combination:</h2>');
fprintf(fid,'\n<img src= " %s/report/figs/coilReconFig.jpg " width="800" height="400"></body>',pwd);
fprintf(fid,'\n\n<p> </p>');
fprintf(fid,'\n\n<h2>Results of removal of bad averages:</h2>');
fprintf(fid,'\n<p>Original number of averages: \t%5.6f </p>',raw.sz(raw.dims.averages)*raw.sz(raw.dims.subSpecs));
fprintf(fid,'\n<p>Number of bad Averages removed:  \t%5.6f </p>',nBadAvgTotal);
fprintf(fid,'\n<p>Number of remaining averages in processed dataset:  \t%5.6f </p>',out_rm.sz(out_rm.dims.averages)*raw.sz(raw.dims.subSpecs));
fprintf(fid,'\n<p>Bad Averages Removal Threshold was:  \t%2.2f </p>',nsd);
fprintf(fid,'\n<img src= " %s/report/figs/rmBadAvg_prePostFig.jpg " width="800" height="600"><img src= " %s/report/figs/rmBadAvg_scatterFig.jpg " width="800" height="400">',pwd,pwd);
fprintf(fid,'\n\n<p> </p>');
fprintf(fid,'\n\n<h2>Results of spectral registration:</h2>');
fprintf(fid,'\n<p>Total frequency drift was: \t%5.6f </p>',max(totalFreqDrift));
fprintf(fid,'\n<p>Total phase drift was: \t%5.6f </p>',max(totalPhaseDrift));
fprintf(fid,'\n<img src= " %s/report/figs/alignAvgs_prePostFig.jpg " width="800" height="600">',pwd);
fprintf(fid,'\n\n<p> </p>');
fprintf(fid,'\n<img src= " %s/report/figs/freqDriftFig.jpg " width="400" height="400"><img src="%s/report/figs/phaseDriftFig.jpg " width="400" height="400">',pwd,pwd);
fprintf(fid,'\n\n<p> </p>');
fprintf(fid,'\n\n<h2>Final Result:</h2>');
fprintf(fid,'\n<img src= " %s/report/figs/finalSpecFig.jpg " width="800" height="400">',pwd);
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



