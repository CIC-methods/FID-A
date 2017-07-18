% run_specialproc_fmrs_slidingWindow.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% [out1,out_w]=run_specialproc_fmrs_slidingWindow(filestring,windowSize);
% 
% DESCRIPTION:
% Processing script for functional MRS data acquired using the SPECIAL MRS 
% sequence.  This script accepts data in Siemens .dat format (twix raw data).  
% Processing steps include combination of reciever channels, removal of bad 
% averages, freqeuncy drift correction, and leftshifting.  This code
% generates a 'sliding window timecourse' of MR spectra by combining the
% averages within a small window given by the windowSize argument, and then
% sliding the window by 1 average and combining again.  Each summed window
% is output as an LCModel text file to be analyzed in LCModel.  As a
% result, this function generates many text output files.  
% 
% INPUTS:
% filestring       = String variable for the name of the directory containing
%                       the water suppressed .dat file.  Water unsuppressed
%                       .dat file should be contained in [filestring '_w/'];
% windowSize       = This is an integer that specifies the number of averages
%                       that are stored within the sliding window.  It is
%                       recommended to choose a window size that is
%                       divisible by the number of phase cycles so that the
%                       window does not contain any partial phase cycles.
% 
% OUTPUTS:
% out1             = The first sliding window spectrum (the others are 
%                    written to LCModel format).  
% out_w            = Fully processed, water unsuppressed output spectrum.

function [out1,out_w]=run_specialproc_fmrs_slidingWindow(filestring,windowSize);

if nargin<3
    leadingAvgsToRmv=0;
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
        out_cs2=out_cs;
        while nbadAverages>0;
            [out_rm,metric{iter},badAverages]=op_rmbadaverages(out_cs2,nsd,'t');
            badAverages;
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

%Do a phase correction
SpecTool(op_leftshift(op_averaging(out_aa),out_aa.pointsToLeftshift),0.05,-2,7);
ph0=input('input 0 order phase correction: ');
ph1=input('input 1st order phase correction: ');
out_aa=op_addphase(out_aa,ph0,ph1);


%Now use op_takeaverages to perform the sliding window:
windowIndices=zeros(out_aa.sz(out_aa.dims.averages),1);
slidingWindow=[1:windowSize];
windowIndices(slidingWindow)=ones(length(slidingWindow),1);
windowIndices=windowIndices>0; %make logical
close all;
wrt=input('write? ','s');
for n=1:out_aa.sz(out_aa.dims.averages)-(length(slidingWindow)-1)
    %take sliding window:
    eval(['out' num2str(n) '_aa=op_takeaverages(out_aa,windowIndices);']);
    
    %Averaging and left shifting
    eval(['out' num2str(n) '=op_leftshift(op_averaging(out' num2str(n) '_aa),out' num2str(n) '_aa.pointsToLeftshift);']);
    
    if wrt=='y' || wrt=='Y'
        %Write to LCModel text file
        eval(['RF=io_writelcm(out' num2str(n) ',[filestring ''/'' filestring ''_out' num2str(n) '_lcm''],8.5);']);
    end
    
    windowIndices(slidingWindow)=zeros(length(slidingWindow),1);
    slidingWindow=slidingWindow+1;
    windowIndices(slidingWindow)=ones(length(slidingWindow),1);
end

%now do the averaging and left shift for water unsuppressed data to get rid of first order phase:
if exist(filename2)
    out_w_av=op_leftshift(op_averaging(out_w_aa),out_w_aa.pointsToLeftshift);
end

%Now do a manual phase correction on the water unsuppressed data:
if exist(filename2)
    SpecTool(out_w_av,0.05,-2,7);
    ph0=input('input 0 order phase correction: ');
    ph1=input('input 1st order phase correction: ');
    
    out_w=op_addphase(out_w_av,ph0,ph1);
end

if wrt=='y' || wrt=='Y'
    if exist(filename2)
            RF=io_writelcm(out_w,[filestring '_w/' filestring '_w_lcm'],8.5);
    end
end










