% run_pressproc.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% [out,out_w,out_noproc,out_w_noproc]=run_pressproc(filestring,aaDomain,tmaxin,iterin);
% 
% DESCRIPTION:
% Processing script for Siemens PRESS MRS data in .dat format (twix raw data).  
% Includes combination of reciever channels, removal of bad averages, 
% freqeuncy drift correction, and leftshifting.
% 
% INPUTS:
% filestring    = String variable for the name of the directory containing
%                   the water suppressed .dat file.  Water unsuppressed
%                   .dat file should be contained in [filestring '_w/'];
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

function [out,out_w,out_noproc,out_w_noproc]=run_pressproc(filestring,aaDomain,tmaxin,iterin);

if nargin<4
    iterin=20;
    if nargin<3
        tmaxin=0.2;
        if nargin<2
            aaDomain='f';
        end
    end
end


%Find the filename of the first MEGA_GABA dataset
close all
unixString=['ls ' filestring '/*.dat'];
[status, filename1]=unix(unixString);
filename1=filename1(1:end-1);

unixStringw=['ls ' filestring '_w/*.dat'];
[status,filename2]=unix(unixStringw);
filename2=filename2(1:end-1);

%read in the data.
out_raw=io_loadspec_twix(filename1);

%load water unsuppressed data and find the coil phases:
if exist(filename2)
    out_w_raw=io_loadspec_twix(filename2);
    coilcombos=op_getcoilcombos(out_w_raw,out_w_raw.pointsToLeftshift+1);
else
    coilcombos=op_getcoilcombos(op_averaging(out_raw),out_raw.pointsToLeftshift+1);
end

%now combine the coil channels:
[out_cc,fid_pre,spec_pre,ph,sig]=op_addrcvrs(out_raw,out_raw.pointsToLeftshift+1,'w',coilcombos);
if exist(filename2)
    [out_w_cc,fid_w_pre,spec_w_pre,ph_w,sig_w]=op_addrcvrs(out_w_raw,out_w_raw.pointsToLeftshift+1,'w',coilcombos);
end
%make the un-processed spectra:
out_noproc=op_leftshift(op_averaging(out_cc),out_cc.pointsToLeftshift);
if exist(filename2)
    out_w_noproc=op_leftshift(op_averaging(out_w_cc),out_w_cc.pointsToLeftshift);
end

%plot the data before and after coil phasing:
figure('position',[0 0 560 420]);
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
    figure('position',[0 500 560 420]);
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

%%%%%%%%%%%%%%%%%%%%%OPTIONAL REMOVAL OF BAD AVERAGES%%%%%%%%%%%%%%%%
close all
figure('position',[0 0 560 420]);
plot(out_cc.ppm,out_cc.specs);
xlabel('Frequency (ppm)');
ylabel('Amplitude (a.u.)');
title('Water suppressed spectra (all averages)');
       
out_cc2=out_cc;
nBadAvgTotal=0;
nbadAverages=1;
rmbadav=input('would you like to remove bad averages?  ','s');
close all;
if rmbadav=='n' || rmbadav=='N'
    out_rm=out_cc;
else
    sat='n'
    while sat=='n'||sat=='N'
        nsd=input('input number of standard deviations.  ');
        iter=1;
        nbadAverages=1;
        nBadAvgTotal=0;
        out_cc2=out_cc;
        while nbadAverages>0;
            [out_rm,metric{iter},badAverages]=op_rmbadaverages(out_cc2,nsd,'t');
            badAverages
            nbadAverages=length(badAverages);
            nBadAvgTotal=nBadAvgTotal+nbadAverages
            out_cc2=out_rm;
            iter=iter+1;
            pause
            close all;
        end
        figure('position',[0 0 560 420]);
        subplot(1,2,1);
        plot(out_cc.ppm,out_cc.specs);xlim([1 5]);
        xlabel('Frequency (ppm)');
        ylabel('Amplitude (a.u.)');
        title('Before removal of bad averages:');
        subplot(1,2,2);
        plot(out_rm.ppm,out_rm.specs);xlim([1 5]);
        xlabel('Frequency (ppm)');
        ylabel('Amplitude (a.u.)');
        title('After removal of bad averages:');
        figure('position',[0 500 560 420]);
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
fprintf(fid,'Original number of averages: \t%5.6f',out_raw.sz(out_raw.dims.averages));
disp(['Original number of averages:  ' num2str(out_raw.sz(out_raw.dims.averages))]);
fprintf(fid,'\nNumber of bad Averages removed:  \t%5.6f',nBadAvgTotal);
disp(['Number of bad averages removed:  ' num2str(nBadAvgTotal)]);
fprintf(fid,'\nNumber of remaining averages in processed dataset:  \t%5.6f',out_rm.sz(out_rm.dims.averages));
disp(['Number of remaining averages in processed dataset:  ' num2str(out_rm.sz(out_rm.dims.averages))]);
fclose(fid);
%%%%%%%%%%%%%%%%%%%%END OF BAD AVERAGES REMOVAL%%%%%%%%%%%%%%%%%%%%



%now align averages;
driftCorr=input('Would you like to perform the frequency drift correction?  ','s');
if driftCorr=='n'|| driftCorr=='N'
    out_aa=out_rm;
    if exist(filename2)
        out_w_aa=out_w_cc;
    end
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
                %tmax=input('input tmax for drift correction: ');
                tmax=tmaxin+0.04*randn(1);
                %fmin=input('input fmin for drift correction: ');
                fmin=1.8+0.04*randn(1);
                %fmax=input('input fmax for drift correction: ');
                fmaxarray=[3.5,3.5,4.2,4.2,7.5,4.2];
                fmax=fmaxarray(randi(6,1))
                [out_aa,fs,phs]=op_alignAverages_fd(out_rm2,fmin,fmax,tmax,'n');
            end
            if exist(filename2)
                [out_w_aa,fs_w,phs_w]=op_alignAverages(out_w_cc,5*tmax,'n');
                %[out_w_aa,fs_w,phs_w]=op_alignAverages(out_w_aa,0.5,'n');
            end
            
            fsCum=fsCum+fs;
            phsCum=phsCum+phs;
            
            
            fsPoly=polyfit([1:out_aa.sz(out_aa.dims.averages)]',fs,1)
            phsPoly=polyfit([1:out_aa.sz(out_aa.dims.averages)]',phs,1)
            iter
            
            
            out_rm2=out_aa;
            if exist(filename2);
                out_w_cc=out_w_aa;
            end
            
            
            iter=iter+1;
        end
        %Now plot the cumulative frequency drift correction:
        figure('position',[0 500 1125 420]);
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
        
        figure('position',[0 0 560 420]);
        plot([1:out_aa.sz(out_aa.dims.averages)],phsCum);
        xlabel('Scan Number');
        ylabel('Phase Drift (deg.)');
        title('Estimated Phase Drift');
        
        figure('position',[570 0 560 420]);
        plot([1:out_aa.sz(out_aa.dims.averages)],fsCum);
        xlabel('Scan Number');
        ylabel('Frequency Drift (Hz)');
        title('Estimated Frequency Drift');
        
        sat=input('Are you satisfied with the frequency drift correction? ','s');
    end
    
end
close all;
%now do the averaging and leftshift to get rid of first order phase:
out=op_leftshift(op_averaging(out_aa),out_aa.pointsToLeftshift);
if exist(filename2)
    out_w=op_leftshift(op_averaging(out_w_aa),out_w_aa.pointsToLeftshift);
end

SpecTool(out,0.05,-2,7);
ph0=input('input 0 order phase correction: ');
ph1=input('input 1st order phase correction: ');

out=op_addphase(out,ph0,ph1);
out_noproc=op_addphase(out_noproc,ph0,ph1);

if exist(filename2)
    SpecTool(out_w,0.05,-2,7);
    ph0=input('input 0 order phase correction: ');
    ph1=input('input 1st order phase correction: ');

    out_w=op_addphase(out_w,ph0,ph1);
    out_w_noproc=op_addphase(out_w_noproc,ph0,ph1);
end

wrt=input('write? ','s');
if wrt=='y' || wrt=='Y'
    RF=io_writelcm(out,[filestring '/' filestring '_lcm'],out.te);
    RF=io_writelcm(out_noproc,[filestring '/' filestring '_lcm_unprocessed'],out_noproc.te);
    if exist(filename2)
        RF=io_writelcm(out_w,[filestring '_w/' filestring '_w_lcm'],out_w.te);
        RF=io_writelcm(out_w_noproc,[filestring '_w/' filestring '_w_lcm_unprocessed'],out_w_noproc.te);
    end
end




