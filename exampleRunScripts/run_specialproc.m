% run_specialproc.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% [out,out_w,out_noproc,out_w_noproc]=run_specialproc(filestring,aaDomain,tmaxin,iterin);
% 
% DESCRIPTION:
% Processing script for Siemens SPECIAL MRS data in .dat format (twix raw 
% data).  Includes combination of reciever channels, removal of bad averages, 
% freqeuncy drift correction, and leftshifting.
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

function [out,out_w,out_noproc,out_w_noproc]=run_specialproc(filestring,aaDomain,tmaxin,iterin);

if nargin<4
    iterin=20;
    if nargin<3
        tmaxin=0.2;
        if nargin<2
            aaDomain='f';
        end
    end
end

%Find the filename of the SPECIAL dataset
close all
unixString1=[filestring '/*.dat'];
[filename]=dir(unixString1);
filename=filename.name(1:end);

unixStringw=[filestring '_w/*.dat'];
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
out_raw=io_loadspec_twix([filestring '/' filename]);

%load water unsuppressed data and find the coil phases:
if water
    disp('***FOUND WATER UNSUPPRESSED DATA***');
    out_w_raw=io_loadspec_twix([filestring '_w/' filenamew]);
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

if water
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

%%%%%%%%%%%%%%%%%%%%%OPTIONAL ALIGNMENT OF SUBSPECTRA%%%%%%%%%%%%%%%%
fs_ai=[];
phs_ai=[];
alignISIS=input('would you like to align subspectra?  ','s');
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
    figure('position',[0 0 560 420]);
    subplot(1,2,1);
    plot(out_cc_temp.ppm,out_cc_temp.specs);
    xlim([0 5]);
    xlabel('Frequency (ppm)');
    ylabel('Amplitude (a.u.)');
    title('Subspecs not aligned: (all averages)');
    subplot(1,2,2);
    plot(out_ai_cc_temp.ppm,out_ai_cc_temp.specs);
    xlim([0 5]);
    xlabel('Frequency (ppm)');
    ylabel('Amplitude (a.u.)');
    title('Subspecs aligned: (all averages)');
    
    if water
        figure('position',[0 550 560 420]);
        subplot(1,2,1);
        plot(out_w_cc_temp.ppm,out_w_cc_temp.specs);
        xlim([3.7 5.7]);
        xlabel('Frequency (ppm)');
        ylabel('Amplitude (a.u.)');
        title('Subspecs not aligned: (all averages)');
        subplot(1,2,2);
        plot(out_w_ai_cc_temp.ppm,out_w_ai_cc_temp.specs);
        xlim([3.7 5.7]);
        xlabel('Frequency (ppm)');
        ylabel('Amplitude (a.u.)');
        title('Subspecs aligned: (all averages)');
    end
    
    figure('position',[570 50 560 420]);
    subplot(2,1,1);
    plot([1:length(fs_ai)],fs_ai);
    xlabel('Scan Number');
    ylabel('Frequency Drift (Hz)');
    title('Estimated Frequency Drift');
    subplot(2,1,2);
    plot([1:length(phs_ai)],phs_ai);
    xlabel('Scan Number');
    ylabel('Phase Drift (Degrees)');
    title('Estimated Phase Drift');
    
    
    sat=input('are you satisfied with alignment of subspecs? ','s');
    if strcmp(sat,'n') || strcmp(sat,'N')
        out_ai=out_cc;
        if water
            out_w_ai=out_w_cc;
        end
    end
    
else
    out_ai=out_cc;
    if water
        out_w_ai=out_w_cc;
    end
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
            badAverages;
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
driftCorr=input('Would you like to perform the frequency drift correction?  ','s');
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



%now do the averaging and left shift to get rid of first order phase:
out_av=op_leftshift(op_averaging(out_aa),out_aa.pointsToLeftshift);
if water
    out_w_av=op_leftshift(op_averaging(out_w_aa),out_w_aa.pointsToLeftshift);
end

%Do a manual phase correction:
SpecTool(out_av,0.05,-2,7);
ph0=input('input 0 order phase correction: ');
ph1=input('input 1st order phase correction: ');

out=op_addphase(out_av,ph0,ph1);
out_noproc=op_addphase(out_noproc,ph0,ph1);


%Now do a manual phase correction on the water unsuppressed data:
if water
    SpecTool(out_w_av,0.05,-2,7);
    ph0=input('input 0 order phase correction: ');
    ph1=input('input 1st order phase correction: ');
    
    out_w=op_addphase(out_w_av,ph0,ph1);
    out_w_noproc=op_addphase(out_w_noproc,ph0,ph1);
    
end

wrt=input('write? ','s');
if wrt=='y' || wrt=='Y'
    RF=io_writelcm(out,[filestring '/' filestring '_lcm'],out.te);
    RF=io_writelcm(out_noproc,[filestring '/' filestring '_lcm_unprocessed'],out_noproc.te);
    if water
        RF=io_writelcm(out_w,[filestring '_w/' filestring '_w_lcm'],out_w.te);
        RF=io_writelcm(out_w_noproc,[filestring '_w/' filestring '_w_unprocessed'],out_w_noproc.te);
    end
end











