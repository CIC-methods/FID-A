% op_alignAverages_fd.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% [out,fs,phs]=op_alignAverages_fd(in,minppm,maxppm,tmax,med,ref);
% 
% DESCRIPTION:
% Perform time-domain spectral registration using a limited range of
% frequencies to correct frequency and phase drifts.  As described in Near
% et al.  Frequency and phase drift correction of magnetic resonance 
% spectroscopy data by spectral registration in the time domain. Magn Reson 
% Med 2015; 73(1):44-50.
% 
% INPUTS:
% in        = Input data structure.
% minppm	= Minimum of frequency range (ppm).
% maxppm	= Maximum of frequnecy range (ppm).
% tmax      = Maximum time (s) in time domain to use for alignment.
% med       = Align averages to the median of the averages? ('y','n', 'a' or 
%             'r').  If you select 'n', all averages will be aligned to a 
%             single average.  The average chosen as the reference 
%             average will be the one with the lowest 'unlikeness' metric 
%             (see 'op_rmbadaverages.m').  If you select 'y', all
%             averages will be aligned to the median of the averages.  If
%             you select 'a', all averages will be aligned to the average
%             of the averages.  If you select 'r', all averages will be 
%             aligned to an externally provided reference spectrum.
% ref       = An externally provided reference spectrum that you would like
%             to align everything to (Required only if med = 'r').  
%
% OUTPUTS:
% out       = Output following alignment of averages.  
% fs        = Vector of frequency shifts (in Hz) used for alignment.
% phs       = Vector of phase shifts (in degrees) used for alignment.

function [out,fs,phs]=op_alignAverages_fd(in,minppm,maxppm,tmax,med,ref)

if ~in.flags.addedrcvrs
    error('ERROR:  I think it only makes sense to do this after you have combined the channels using op_addrcvrs.  ABORTING!!');
end

if (strcmp(med,'r') || strcmp(med,'R'))
    if nargin<6
        error('ERROR:  If using the ''r'' option for input variable ''med'', then a 6th input argument must be provided');
    end
else
    if nargin<6
        ref=struct();
    end
end

parsFit=[0,0];

if in.dims.subSpecs==0
    B=1;
else
    B=in.sz(in.dims.subSpecs);
end

fs=zeros(in.sz(in.dims.averages),B);
phs=zeros(in.sz(in.dims.averages),B);
fids=zeros(in.sz(in.dims.t),1,B);
for m=1:B
    if med=='y' || med=='Y'
        disp('Aligning all averages to the Average of the averages.');
        base=op_averaging(in);
        base=op_freqrange(base,minppm,maxppm);
        base=[real(base.fids(base.t>=0 & base.t<tmax,m));imag(base.fids(base.t>=0 & base.t<tmax,m))];
        ind_min=0;
    elseif med=='n' || med=='N'
        %First find the average that is most similar to the total average:
        inavg=op_median(in);
        for k=1:in.sz(in.dims.averages)
            for l=1:B
                metric(k,l)=sum((real(in.fids(in.t>=0 & in.t<=tmax,k,l))-(real(inavg.fids(inavg.t>=0 & inavg.t<=tmax,l)))).^2);
            end
        end
        [temp,ind_min]=min(metric(:,m));
        
        %Now set the base function using the index of the most similar
        %average:
        disp(['Aligning all averages to average number ' num2str(ind_min) '.']);
        base=op_freqrange(in,minppm,maxppm);
        base=[real(base.fids(base.t>=0 & base.t<tmax,ind_min,m));imag(base.fids(base.t>=0 & base.t<tmax,ind_min,m))];
        fids(:,ind_min,m)=in.fids(:,ind_min,m);
    elseif med=='r' || med=='R'
        disp('Aligning all averages to an externally provided reference spectrum.');
        base=ref;
        base=op_freqrange(base,minppm,maxppm);
        base=[real(base.fids(base.t>=0 & base.t<tmax,m));imag(base.fids(base.t>=0 & base.t<tmax,m))];
        ind_min=0;
    end
    for n=1:in.sz(in.dims.averages)
        if n~=ind_min
            parsGuess=parsFit;
            %parsGuess(1)=parsGuess(1);
            %disp(['fitting subspec number ' num2str(m) ' and average number ' num2str(n)]);
            datarange=op_freqrange(in,minppm,maxppm);
            start=datarange.fids(datarange.t>=0 & datarange.t<tmax,n,m);
            parsFit=nlinfit(start,base,@op_freqPhaseShiftComplexRangeNest,parsGuess);
            fids(:,n,m)=op_freqPhaseShiftNest(parsFit,in.fids(:,n,m));
            fs(n,m)=parsFit(1);
            phs(n,m)=parsFit(2);
            %plot(in.ppm,fftshift(ifft(fids(:,1,m))),in.ppm,fftshift(ifft(fids(:,n,m))));
        end
    end
end


%re-calculate Specs using fft
specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);


%FILLING IN DATA STRUCTURE
out=in;
out.fids=fids;
out.specs=specs;

%FILLING IN THE FLAGS
out.flags=in.flags;
out.flags.writtentostruct=1;
out.flags.freqcorrected=1;

    function y=op_freqPhaseShiftComplexRangeNest(pars,input)
        f=pars(1);     %Frequency Shift [Hz]
        p=pars(2);     %Phase Shift [deg]
        
        
        dwelltime=datarange.dwelltime;
        t=0:dwelltime:(length(input)-1)*dwelltime;
        fid=input(:);
        
        shifted=addphase(fid.*exp(1i*t'*f*2*pi),p);
        
        y=[real(shifted);imag(shifted)];
        %y=real(fid.*exp(-1i*t'*f*2*pi));
        
    end

    function y=op_freqPhaseShiftNest(pars,input)
        f=pars(1);     %Frequency Shift [Hz]
        p=pars(2);     %Phase Shift [deg]
        
        
        dwelltime=in.dwelltime;
        t=0:dwelltime:(length(input)-1)*dwelltime;
        fid=input(:);
        
        y=addphase(fid.*exp(1i*t'*f*2*pi),p);
        %y=real(fid.*exp(-1i*t'*f*2*pi));
        
    end
end
