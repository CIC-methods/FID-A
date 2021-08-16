% run_megapressproc_auto.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% [diffSpec,sumSpec,subSpec1,subSpec2,outw]=run_megapressproc_auto(filestring,coilcombos,avgAlignDomain,alignSS);
% 
% DESCRIPTION:
% Fully automated processing script for Siemens MEGA-PRESS MRS data in .dat 
% format (twix raw data).  Includes combination of reciever channels, 
% removal of bad averages, freqeuncy drift correction, manual alignment of 
% edit-on and edit-off spectra, and leftshifting.  This pipeline requires
% no user interaction. This function automatically 
% generates an html report to describe the results of each processing step.
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
% diffSpec           = Fully processed difference spectrum.
% sumSpec            = Fully processed sum spectrum.
% subSpec1           = Fully processed MEGA-PRESS subspectrum #1.
% subSpec2           = Fully processed MEGA-PRESS subspectrum #2.
% outw               = Fully processed water-unsuppressed spectrum.

function [diffSpec,sumSpec,subSpec1,subSpec2,outw]=run_megapressproc_auto(filestring,coilcombos,avgAlignDomain,alignSS)

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

%Find the filename of the MEGA_GABA dataset
close all
unixString=fullfile(filestring,'*.dat');
[filename]=dir(unixString);
filename=filename.name(1:end);

unixStringw=fullfile([filestring '_w'],'*.dat');
[filenamew]=dir(unixStringw);
if ~isempty(filenamew)
    filenamew=filenamew.name(1:end);
    water=true
else
    water=false
    outw=struct();
    outw_noproc=struct();
end

% %read in both datasets:
raw=io_loadspec_twix(fullfile(filestring,filename));
if water
    raww=io_loadspec_twix(fullfile([filestring '_w/'], filenamew));
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
saveas(h,fullfile(reportFigsDir,'coilReconFig'),'jpg');
saveas(h,fullfile(reportFigsDir,'coilReconFig'),'fig');
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
        trackRemovedAvg = 1:out_cc.sz(out_cc.dims.averages);
        while nbadAverages>0
            [out_rm,metric{iter},badAverages]=op_rmbadaverages(out_cc2,nsd,'t');
            trackRemovedAvg(badAverages) = [];            
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
        title('Edit-OFF Before','FontSize',12);
        box off;
        subplot(2,2,2);
        plot(out_rm.ppm,real(out_rm.specs(:,:,1)));xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Edit-OFF After','FontSize',12);
        box off;
        subplot(2,2,3);  
        plot(out_cc.ppm,real(out_cc.specs(:,:,2)));xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Edit-ON Before','FontSize',12);
        box off;
        subplot(2,2,4);
        plot(out_rm.ppm,real(out_rm.specs(:,:,2)));xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Edit-ON After','FontSize',12);
        box off;
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 20 15]);
        saveas(h,fullfile(reportFigsDir,'rmBadAvg_prePostFig'),'jpg');
        saveas(h,fullfile(reportFigsDir,'rmBadAvg_prePostFig'),'fig');
        close(h);
        
        %figure('position',[0 550 560 420]);
        h=figure('visible','off');
        plot(1:out_cc.sz(out_cc.dims.averages),metric{1}(:,1),'.r','MarkerSize',16)
        hold on
        plot(1:out_cc.sz(out_cc.dims.averages),metric{1}(:,2),'.b','MarkerSize',16)
        removedAvgs = find(~ismember(1:out_cc.sz(out_cc.dims.averages),trackRemovedAvg));
        plot(removedAvgs,metric{1}(removedAvgs,1),'ko','MarkerSize',20)
        plot(removedAvgs,metric{1}(removedAvgs,2),'ko','MarkerSize',20)
        set(gca,'FontSize',8);
        xlabel('Scan Number','FontSize',10);
        ylabel('Deviation Metric','FontSize',10);
        legend('Original subSpec 1','Original subSpec 2','Removed Avg');
        legend boxoff;
        title('Deviation Metric','FontSize',12);
        box off;
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 20 10]);
        saveas(h,fullfile(reportFigsDir,'rmBadAvg_scatterFig'),'jpg');
        saveas(h,fullfile(reportFigsDir,'rmBadAvg_scatterFig'),'fig');
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
            iter=iter+1;
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
            p=polyfit(x,fs,1);
            
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
        title('Edit-OFF Before','FontSize',12);
        box off;
        subplot(2,2,2);
        plot(out_aa.ppm,real(out_aa.specs(:,:,1)));xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Edit-OFF After','FontSize',12);
        box off;
        subplot(2,2,3);
        plot(out_rm.ppm,real(out_rm.specs(:,:,2)));xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Edit-ON Before','FontSize',12);
        box off;
        subplot(2,2,4);
        plot(out_aa.ppm,real(out_aa.specs(:,:,2)));xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Edit-ON After','FontSize',12);
        box off;
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 20 15]);
        saveas(h,fullfile(reportFigsDir,'alignAvgs_prePostFig'),'jpg');
        saveas(h,fullfile(reportFigsDir,'alignAvgs_prePostFig'),'fig');
        close(h);
        
        h=figure('visible','off');
        plot([1:out_aa.sz(out_aa.dims.averages)],fscum,'.-','LineWidth',2);
        set(gca,'FontSize',8);
        xlabel('Scan Number','FontSize',10);
        ylabel('Frequency Drift [Hz]','FontSize',10);
        box off;
        legend('Edit-off scans','Edit-on scans','Location','SouthEast');
        legend boxoff;
        title('Estimated Freqeuncy Drift','FontSize',12);
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 10 10]);
        saveas(h,fullfile(reportFigsDir,'freqDriftFig'),'jpg');
        saveas(h,fullfile(reportFigsDir,'freqDriftFig'),'fig');
        close(h);
        
        h=figure('visible','off');
        plot([1:out_aa.sz(out_aa.dims.averages)],phscum,'.-','LineWidth',2);
        set(gca,'FontSize',8);
        xlabel('Scan Number','FontSize',10);
        ylabel('Phase Drift [Deg.]','FontSize',10);
        box off;
        legend('Edit-off scans','Edit-on scans');
        legend boxoff;
        title('Estimated Phase Drift','FontSize',12);
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 10 10]);
        saveas(h,fullfile(reportFigsDir,'phaseDriftFig'),'jpg');
        saveas(h,fullfile(reportFigsDir,'phaseDriftFig'),'fig');
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

%now do automatic zero-order phase correction (Use Creatine Peak):
out_ls_ss2=op_takesubspec(out_ls,2);
%SpecTool(out_ls_ss2,0.2,1,3.5);
[out_ls_ss2_ph,ph0]=op_autophase(out_ls_ss2,2.9,3.1);
%ph0=input('Input zero-order phase:');
out_ph=op_addphase(out_ls,ph0);

%do same phase corection on unprocessed data
out_noproc=op_addphase(out_noproc,ph0);

%Now align subspecs:
if out_ph.dims.subSpecs
    %Make water suppressed subspecs 'out of phase'.
    out=op_alignMPSubspecs(out_ph);
else
    out=out_ph;    
end
if water
    if outw_ls.dims.subSpecs
        outw_as=op_alignMPSubspecs_fd(outw_ls,3,6.5,'i');
    else
        outw_as=outw_ls;
    end
end
        
%Make fully processed data;
diffSpec=op_combinesubspecs(out,'diff');
sumSpec=op_combinesubspecs(out,'summ');
subSpec1=op_takesubspec(out,1);
subSpec1=op_addphase(subSpec1,180);
subSpec2=op_takesubspec(out,2);

%Frequency shift all spectra so that Creatine appears at 3.027 ppm:
[subSpec1,frqShift]=op_ppmref(subSpec1,2.9,3.1,3.027);
diffSpec=op_freqshift(diffSpec,frqShift);
sumSpec=op_freqshift(sumSpec,frqShift);
subSpec2=op_freqshift(subSpec2,frqShift);

%Make final water unsuppressed data
if water
    if outw_ls.dims.subSpecs
        outw=op_combinesubspecs(outw_as,'diff');
    else
        outw=outw_ls;
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
legend('Edit-OFF','Edit-ON');
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
saveas(h,fullfile(reportFigsDir,'finalSpecFig'),'jpg');
saveas(h,fullfile(reportFigsDir,'finalSpecFig'),'fig');

% writ=input('Write? (y or n):  ','s');
writ='y';
if writ=='y' || writ=='Y'
     RF=io_writelcm(diffSpec,fullfile(filestring,'diff_lcm'),diffSpec.te);
     RF=io_writelcm(subSpec1,fullfile(filestring,'editOFF_lcm'),subSpec1.te);
     RF=io_writelcm(subSpec2,fullfile(filestring,'editON_lcm'),subSpec2.te);
    if water
        RF=io_writelcm(outw,fullfile([filestring '_w/'],'w_lcm'),outw.te);
    end
end

close all;

%write an html report: 
fid=fopen(fullfile(reportDir,'report.html'),'w+');
fprintf(fid,'<!DOCTYPE html>');
fprintf(fid,'\n<html>');
logoPath=which('FID-A_LOGO.jpg');
fprintf(fid,'\n<img src= " %s " width="120" height="120"></body>',logoPath);
fprintf(fid,'\n<h1>FID-A Processing Report</h1>');
fprintf(fid,'\n<h2>Processing pipeline applied to MEGA-PRESS data using run_megapressproc_auto.m</h2>');
fprintf(fid,'\n<p>FILENAME: %s </p>',fullfile(filestring,filename));
fprintf(fid,'\n<p>DATE: %s </p>',date);
fprintf(fid,'\n\n<p> </p>');
fprintf(fid,'\n\n<h2>Results of multi-coil combination:</h2>');
fprintf(fid,'\n<img src= " %s " width="800" height="400"></body>',fullfile(reportFigsDir,'coilReconFig.jpg'));
fprintf(fid,'\n\n<p> </p>');
fprintf(fid,'\n\n<h2>Results of removal of bad averages:</h2>');
fprintf(fid,'\n<p>Original number of averages: \t%5.6f </p>',raw.sz(raw.dims.averages)*raw.sz(raw.dims.subSpecs));
fprintf(fid,'\n<p>Number of bad Averages removed:  \t%5.6f </p>',nBadAvgTotal);
fprintf(fid,'\n<p>Number of remaining averages in processed dataset:  \t%5.6f </p>',out_rm.sz(out_rm.dims.averages)*raw.sz(raw.dims.subSpecs));
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



