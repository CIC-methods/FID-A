% run_specialproc_fmrs.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% [out_stimOFF,out_stimON,out_w]=run_specialproc_fmrs(filestring,blockDesign,leadingAvgsToRmv);
% 
% DESCRIPTION:
% Processing script for functional MRS data acquired using the SPECIAL MRS 
% sequence.  This script accepts data in Siemens .dat format (twix raw data).  
% Processing steps include combination of reciever channels, removal of bad 
% averages, freqeuncy drift correction, and leftshifting.  Also includes
% prior knowledge of the stimulation paradigm (given by the 'blockDesign'
% vector) and returns the averaged stimulus OFF and simulus ON spectra.
% 
% INPUTS:
% filestring       = String variable for the name of the directory containing
%                       the water suppressed .dat file.  Water unsuppressed
%                       .dat file should be contained in [filestring '_w/'];
% blockDesign      = This is a vector of positive and negative even integers that
%                       make up the ON/OFF block design.  Each integer 
%                       represents the number of sequential averages in a 
%                       block.  Positive integers refer to ON blocks, and 
%                       negative integers refer to OFF blocks.  For
%                       example, for a block design consisting of 30 OFF
%                       averages followed by 20 ON averages followed by 10
%                       OFF averages, the blockDesign vector would be: 
%                       [-30 20 -10];
% leadingAvgsToRmv = The number of averages to omit from the beginning of
%                       each block.  This is done to account for a lag in the 
%                       neurochemical response to stimulus.  Must be an 
%                       even integer.  (optional. Default=0);
% 
% OUTPUTS:
% out_stimOFF      = Fully processed water suppressed spectrum from the sum 
%                       of the stimulus OFF periods.
% out_stimON       = Fully processed water suppressed spectrum from the sum
%                       of the stimulus ON periods.
% out_w            = Fully processed, water unsuppressed output spectrum.

function [out_stimOFF,out_stimON,out_w]=run_specialproc_fmrs(filestring,blockDesign,leadingAvgsToRmv);

if nargin<3
    leadingAvgsToRmv=0;
end

%Check to make sure the elements of blockDesign are even integers:
if sum(mod(blockDesign,2))>0
    error('ERROR:  blockDesign must consist of only even integers!  ');
end
%Check to make sure that leadingAvgsToRmv is even:
if mod(leadingAvgsToRmv,2)>0
    error('ERROR: leadingAvgsToRmv must be an even integer!  ');
end

%initialize some parameters.
iterin=20;
tmaxin=0.2;
aaDomain='f';
      
%Find the filename of the first MEGA_GABA dataset
close all
unixString=['ls ' filestring '/*.dat'];
[status, filename1]=unix(unixString);
filename1=filename1(1:end-1);

unixStringw=['ls ' filestring '_w/*.dat'];
[status,filename2]=unix(unixStringw);
filename2=filename2(1:end-1);

%read in the data:
out_raw=io_loadspec_twix(filename1);

%load water unsuppressed data and find the coil phases:
if exist(filename2)
    disp('***FOUND WATER UNSUPPRESSED DATA***');
    out_w_raw=io_loadspec_twix(filename2);
    coilcombos=op_getcoilcombos(op_combinesubspecs(out_w_raw,'diff'),2);
else
    disp('***WATER UNSUPPRESSED DATA NOT FOUND***');
    coilcombos=op_getcoilcombos_specReg(op_combinesubspecs(op_averaging(out_raw),'diff'),0,0.01,2);
end

%now combine the coil channels:
[out_cc,fid_pre,spec_pre,ph,sig]=op_addrcvrs(out_raw,2,'w',coilcombos);
if exist(filename2)
    [out_w_cc,fid_w_pre,spec_w_pre,ph_w,sig_w]=op_addrcvrs(out_w_raw,2,'w',coilcombos);
end

%make the un-processed spectra, which may be optionally output from this function:
out_noproc=op_combinesubspecs(op_averaging(out_cc),'diff');
if exist(filename2)
    out_w_noproc=op_combinesubspecs(op_averaging(out_w_cc),'diff');
end

%plot the data before and after coil phasing:
figure('position',[0 50 560 420]);
subplot(2,1,1);
plot(out_raw.ppm,out_raw.specs(:,:,1,1));xlim([-1 7]);
xlabel('Frequency (ppm)');
ylabel('Amplitude (a.u.)');
title('Multi-channel (water supp.) data before phase correction');
subplot(2,1,2);
plot(out_raw.ppm,spec_pre(:,:,1,1));xlim([-1 7]);
xlabel('Frequency (ppm)');
ylabel('Amplitude (a.u.)');
title('Multi-channel (water supp.) data after phase correction');

if exist(filename2)
    figure('position',[0 550 560 420]);
    subplot(2,1,1);
    plot(out_w_raw.ppm,out_w_raw.specs(:,:,1,1));xlim([4 5]);
    xlabel('Frequency (ppm)');
    ylabel('Amplitude (a.u.)');
    title('Multi-channel (water unsupp.) data before phase correction');
    subplot(2,1,2);
    plot(out_w_raw.ppm,spec_w_pre(:,:,1,1));xlim([4 5]);
    xlabel('Frequency (ppm)');
    ylabel('Amplitude (a.u.)');
    title('Multi-channel (water unsupp.) data after phase correction');
end
pause;
close all;

%Now combine the subspecs
out_cs=op_combinesubspecs(out_cc,'diff');
if exist(filename2)
    out_w_cs=op_combinesubspecs(out_w_cc,'diff');
end

%%%%%%%%%%%%%%%%%%%%%OPTIONAL REMOVAL OF BAD AVERAGES%%%%%%%%%%%%%%%%
close all
figure('position',[0 50 560 420]);
plot(out_cs.ppm,out_cs.specs);
xlabel('Frequency (ppm)');
ylabel('Amplitude (a.u.)');
title('Water suppressed spectra (all averages)');

out_cs2=out_cs;
nBadAvgTotal=0;
nbadAverages=1;
allAveragesLeft=[1:out_cs.sz(out_cs.dims.averages)]';
allBadAverages=[];
rmbadav=input('would you like to remove bad averages?  ','s');
close all;
if rmbadav=='n' || rmbadav=='N'
    out_rm=out_cs;
else
    sat='n'
    while sat=='n'||sat=='N'
        nsd=input('input number of standard deviations.  ');
        iter=1;
        nbadAverages=1;
        nBadAvgTotal=0;
        allAveragesLeft=[1:out_cs.sz(out_cs.dims.averages)]';
        allBadAverages=[];
        out_cs2=out_cs;
        while nbadAverages>0;
            [out_rm,metric{iter},badAverages]=op_rmbadaverages(out_cs2,nsd,'t');
            allBadAverages=[allBadAverages; allAveragesLeft(badAverages)];
            badavMask_temp=zeros(length(allAveragesLeft),1);
            badavMask_temp(badAverages)=1;
            allAveragesLeft=allAveragesLeft(~badavMask_temp);
            nbadAverages=length(badAverages)*out_raw.sz(out_raw.dims.subSpecs);
            nBadAvgTotal=nBadAvgTotal+nbadAverages;
            out_cs2=out_rm;
            iter=iter+1;
            disp([num2str(nbadAverages) ' bad averages removed on this iteration.']);
            disp([num2str(nBadAvgTotal) ' bad averages removed in total.']);
            disp('Press any key to continue...');
            pause
            close all;
        end
        figure('position',[0 50 560 420]);
        subplot(1,2,1);
        plot(out_cs.ppm,out_cs.specs);xlim([1 5]);
        xlabel('Frequency (ppm)');
        ylabel('Amplitude (a.u.)');
        title('Before removal of bad averages:');
        subplot(1,2,2);
        plot(out_rm.ppm,out_rm.specs);xlim([1 5]);
        xlabel('Frequency (ppm)');
        ylabel('Amplitude (a.u.)');
        title('After removal of bad averages:');
        figure('position',[0 550 560 420]);
        plot([1:length(metric{1})],metric{1},'.r',[1:length(metric{iter-1})],metric{iter-1},'x');
        xlabel('Scan Number');
        ylabel('Unlikeness metric (a.u.)');
        title('Unlikeness metrics before and after bad averages removal');
        legend('before','after');
        legend boxoff;
        sat=input('are you satisfied with bad averages removal? ','s');
    end
end

%write a readme file to record the number of dropped avgs
fid=fopen([filestring '/readme.txt'],'w+');
fprintf(fid,'Original number of averages: \t%5.6f',out_raw.sz(out_raw.dims.averages)*2);
disp(['Original number of averages:  ' num2str(out_raw.sz(out_raw.dims.averages)*2)]);
fprintf(fid,'\nNumber of bad Averages removed:  \t%5.6f',nBadAvgTotal);
disp(['Number of bad averages removed:  ' num2str(nBadAvgTotal)]);
fprintf(fid,'\nNumber of remaining averages in processed dataset:  \t%5.6f',out_rm.sz(out_rm.dims.averages)*2);
disp(['Number of remaining averages in processed dataset:  ' num2str(out_rm.sz(out_rm.dims.averages)*2)]);
fclose(fid);

close all;

%%%%%%%%%%%%%%%%%%%%END OF BAD AVERAGES REMOVAL%%%%%%%%%%%%%%%%%%%%

%now align averages;
driftCorr=input('Would you like to perform the frequency drift correction?  ','s');
if driftCorr=='n'|| driftCorr=='N'
    out_aa=out_rm;
    out_w_aa=out_w_cs;
else
    sat='n'
    while sat=='n' || sat=='N'
        out_rm2=out_rm;
        fsPoly=100;
        phsPoly=1000;
        fsCum=zeros(out_rm2.sz(out_rm2.dims.averages),1);
        phsCum=zeros(out_rm2.sz(out_rm2.dims.averages),1);
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
            if exist(filename2)
                [out_w_aa,fs_w,phs_w]=op_alignAverages(out_w_cs,5*tmax,'n');
            end
            
            fsCum=fsCum+fs;
            phsCum=phsCum+phs;
            fsPoly=polyfit([1:out_aa.sz(out_aa.dims.averages)]',fs,1)
            phsPoly=polyfit([1:out_aa.sz(out_aa.dims.averages)]',phs,1)
            iter
            out_rm2=out_aa;
            if exist(filename2);
                out_w_cs=out_w_aa;
            end
            iter=iter+1;
        end
        %Now plot the cumulative frequency drift correction:
        figure('position',[0 550 1125 420]);
        subplot(2,1,1);
        plot(out_rm.ppm,out_rm.specs);xlim([0 6]);
        xlabel('Frequnecy (ppm)');
        ylabel('Amplitude (a.u.)');
        title('Before drift correction:');
        subplot(2,1,2);
        plot(out_aa.ppm,out_aa.specs);xlim([0 6]);
        xlabel('Frequency (ppm)');
        ylabel('Amplitude (a.u.)');
        title('After drift correction:');
        
        figure('position',[0 50 560 420]);
        plot([1:out_aa.sz(out_aa.dims.averages)],phsCum);
        xlabel('Scan Number');
        ylabel('Phase Drift (deg.)');
        title('Estimated Phase Drift');
        
        figure('position',[570 50 560 420]);
        plot([1:out_aa.sz(out_aa.dims.averages)],fsCum);
        xlabel('Scan Number');
        ylabel('Frequency Drift (Hz)');
        title('Estimated Frequency Drift');
        
        sat=input('Are you satisfied with the drift correction? ','s');
    end
    
end

%Now separate the spectrum into the "Stimulus ON" and "Stimulus OFF"
%blocks.  Remove the correct number of leading averages from each block.

%First check to make sure that the length of the block design is equal to
%the number of averages:
if sum(abs(blockDesign))==out_aa.rawAverages;
    blockDesign=blockDesign/2;
elseif sum(abs(blockDesign))==out_aa.rawAverages/2;
    %do nothing
else
    error('ERROR:  Length of blockDesign must be equal to the number of averages!  ');
end

%Initialize the ONaverages and OFFaverages vectors.  These will be vectors
%of zeros and ones corresponding to the averages that you would like to
%take to make up the ON and OFF averages.
ONaverages=[];
OFFaverages=[];

%If the block design starts with a negative number, then the experiment
%started with an OFF period.  The ONaverages and OFFaverages are filled as
%follows:
if blockDesign(1)<0
    blockDesignTrimmedON=blockDesign-leadingAvgsToRmv/2;
    if isrow(blockDesign)
        blockDesignTrimmedOFF=blockDesign+[0 (leadingAvgsToRmv/2)*ones(1,length(blockDesign)-1)];
    elseif iscolumn(blockDesign)
        blockDesignTrimmedOFF=blockDesign+[0;(leadingAvgsToRmv/2)*ones(length(blockDesign)-1,1)];
    end
    for n=1:length(blockDesignTrimmedON);
        ONaverages=[ONaverages;zeros(abs(blockDesignTrimmedON(n)),1)+mod(n+1,2)];
    end
    for n=1:length(blockDesignTrimmedOFF);
        OFFaverages=[OFFaverages;zeros(abs(blockDesignTrimmedOFF(n)),1)+mod(n,2)];
    end
    ONaverages=ONaverages(1:sum(abs(blockDesign)));
    OFFaverages=OFFaverages(1:sum(abs(blockDesign)));
   
%If the block design starts with a positive number, then the experiment 
%started with an ON period.  The ON averages and OFF averages are filled as
%follows:
elseif blockDesign(1)>0
    blockDesignTrimmedOFF=blockDesign-leadingAvgsToRmv/2;
    if isrow(blockDesign)
        blockDesignTrimmedON=blockDesign+[0 (leadingAvgsToRmv/2)*ones(1,length(blockDesign)-1)];
    elseif iscolumn(blockDesign)
        blockDesignTrimmedON=blockDesign+[0;(leadingAvgsToRmv/2)*ones(length(blockDesign)-1,1)];
    end
    for n=1:length(blockDesignTrimmedON);
        ONaverages=[ONaverages;zeros(abs(blockDesignTrimmedON(n)),1)+mod(n,2)];
    end
    for n=1:length(blockDesignTrimmedOFF);
        OFFaverages=[OFFaverages;zeros(abs(blockDesignTrimmedOFF(n)),1)+mod(n+1,2)];
    end
    ONaverages=ONaverages(1:sum(abs(blockDesign)));
    OFFaverages=OFFaverages(1:sum(abs(blockDesign)));
end

%Now remove the entries from ONaverages and OFFaverages that correspond to
%the bad-averages that were removed.
BadAvgMask=zeros(length(ONaverages),1);
BadAvgMask(allBadAverages)=1;
ONaverages=ONaverages(~BadAvgMask);
OFFaverages=OFFaverages(~BadAvgMask);

%Now use op_takeaverages to make the stimulation ON and stimulation OFF
%spectra.
ONaverages=ONaverages>0; %ONaverages must be a logical;
OFFaverages=OFFaverages>0; %OFFaverages must be a logical;
outON_aa=op_takeaverages(out_aa,ONaverages);
outOFF_aa=op_takeaverages(out_aa,OFFaverages);


%now do the averaging and left shift to get rid of first order phase:
outON_av=op_leftshift(op_averaging(outON_aa),outON_aa.pointsToLeftshift);
outOFF_av=op_leftshift(op_averaging(outOFF_aa),outOFF_aa.pointsToLeftshift);
if exist(filename2)
    out_w_av=op_leftshift(op_averaging(out_w_aa),out_w_aa.pointsToLeftshift);
end

%Do a manual phase correction:
SpecTool(outON_av,0.05,-2,7);
ph0=input('input 0 order phase correction: ');
ph1=input('input 1st order phase correction: ');

out_stimON=op_addphase(outON_av,ph0,ph1);
out_stimOFF=op_addphase(outOFF_av,ph0,ph1);

%Now do a manual phase correction on the water unsuppressed data:
if exist(filename2)
    SpecTool(out_w_av,0.05,-2,7);
    ph0=input('input 0 order phase correction: ');
    ph1=input('input 1st order phase correction: ');
    
    out_w=op_addphase(out_w_av,ph0,ph1);
end

wrt=input('write? ','s');
if wrt=='y' || wrt=='Y'
    RF=io_writelcm(out_stimON,[filestring '/' filestring '_stimON_lcm'],8.5);
    RF=io_writelcm(out_stimOFF,[filestring '/' filestring '_stimOFF_lcm'],8.5);
    if exist(filename2)
        RF=io_writelcm(out_w,[filestring '_w/' filestring '_w_lcm'],8.5);
    end
end











