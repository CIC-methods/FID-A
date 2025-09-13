% op_alignAverages_robust.m
% Jamie Near, McGill University 2014.
% Mark Mikkelsen, Johns Hopkins University 2020.
%
% USAGE:
% [out,fs,phs]=op_alignAverages_robust(in);
%
% DESCRIPTION:
% Perform robust spectral registration in the time domain to correct
% frequency and phase drifts. As described in Mikkelsen et al.
% Correcting frequency and phase offsets in MRS data using robust
% spectral registration. NMR Biomed 2020;e4368.
%
% INPUTS:
% in        = Input data structure.
%
% OUTPUTS:
% out       = Output following alignment of averages.
% fs        = Vector of frequency shifts (in Hz) used for alignment.
% phs       = Vector of phase shifts (in degrees) used for alignment.


function [out,fs,phs]=op_alignAverages_robust(in)

if ~in.flags.addedrcvrs
    error('ERROR:  I think it only makes sense to do this after you have combined the channels using op_addrcvrs.  ABORTING!!');
end

if in.dims.averages==0
    %DO NOTHING
    disp('WARNING:  No averages found.  Returning input without modification!');
    out=in;
    fs=0;
    phs=0;
    
else
    
    parsFit=[0,0];
    
    if in.dims.subSpecs==0
        B=1;
    else
        B=in.sz(in.dims.subSpecs);
    end
    
    fs=zeros(in.sz(in.dims.averages),B);
    phs=zeros(in.sz(in.dims.averages),B);
    fids=zeros(in.sz(in.dims.t),in.sz(in.dims.averages),B);
    MSEfun=@(a,b) sum((a-b).^2)/length(a);
    for m=1:B
        %Use first n points of time-domain data, where n is the last point where abs(diff(mean(SNR))) > 0.5
        signal=abs(in.fids(:,:,m));
        noise=2*std(signal(ceil(0.75*size(signal,1)):end,:));
        SNR=signal./repmat(noise, [in.sz(in.dims.t) 1]);
        SNR=abs(diff(mean(SNR,2)));
        SNR=SNR(in.t<=0.2); %Use no more than 200 ms of data
        tmax=find(SNR>0.5,1,'last');
        if isempty(tmax) || tmax<find(in.t<=0.1,1,'last') %Use at least 100 ms of data
            tmax=find(in.t<=0.1,1,'last');
        end
        tmax=in.t(tmax);
        %Determine optimal alignment order by calculating a similarity metric (mean squared error)
        D=zeros(in.sz(in.dims.averages));
        for ii=1:in.sz(in.dims.averages)
            D(ii,:)=MSEfun(real(in.fids(in.t>=0 & in.t<tmax,ii,m)),real(in.fids(in.t>=0 & in.t<tmax,:,m)));
        end
        D(~D)=NaN;
        [~,alignOrd]=sort(median(D,'omitnan'));
        disp('Aligning all averages to a weighted reference spectrum.');
        base=[real(in.fids(in.t>=0 & in.t<tmax,alignOrd(1),m));imag(in.fids(in.t>=0 & in.t<tmax,alignOrd(1),m))];
        for n=alignOrd
            parsGuess=parsFit;
            %disp(['fitting subspec number ' num2str(m) ' and average number ' num2str(n)]);
            parsFit=nlinfit(in.fids(in.t>=0 & in.t<tmax,n,m),base,@op_freqPhaseShiftComplexNest,parsGuess);
            fids(:,n,m)=op_freqPhaseShiftNest(parsFit,in.fids(:,n,m));
            fs(n,m)=parsFit(1);
            phs(n,m)=parsFit(2);
            %plot(in.ppm,real(fftshift(ifft(fids(:,1,m)))),in.ppm,real(fftshift(ifft(fids(:,n,m)))));
            %Update reference
            w=0.5*corr(base,[real(fids(in.t>=0 & in.t<tmax,n,m));imag(fids(in.t>=0 & in.t<tmax,n,m))]).^2;
            base=(1-w)*base+w*[real(fids(in.t>=0 & in.t<tmax,n,m));imag(fids(in.t>=0 & in.t<tmax,n,m))];
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
    
end


    function y=op_freqPhaseShiftComplexNest(pars,input)
        f=pars(1);     %Frequency Shift [Hz]
        p=pars(2);     %Phase Shift [deg]
        
        
        dwelltime=in.dwelltime;
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
