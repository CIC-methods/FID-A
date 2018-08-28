% op_alignAverages.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% [out,fs,phs]=op_alignAverages(in,tmax,med,ref);
% 
% DESCRIPTION:
% Perform spectral registration in the time domain to correct frequency and
% phase drifts.  As described in Near et al.  Frequency and phase drift 
% correction of magnetic resonance spectroscopy data by spectral 
% registration in the time domain. Magn Reson Med 2015; 73(1):44-50.
%
% June 15th 2017:  Made the tmax and med arguments optional.  If tmax is 
% not specified, the value is determined automatically by finding the time
% at which the SNR of the FID drops permanently below 5.  This idea
% was suggested by Mark Mikkelsen.  Thanks Mark!!
% 
% INPUTS:
% in        = Input data structure.
% tmax      = Maximum time (s) in time domain to use for alignment.
%             (Optional.  Default is the time at which SNR drops below 5)
% med       = Align averages to the median of the averages? ('y','n', 'a', or 
%             'r').  (Optional.  Default = 'n').  If you select 'n', all 
%             averages will be aligned to a single average.  The average 
%             chosen as the reference average will be the one with the 
%             lowest 'unlikeness' metric (see 'op_rmbadaverages.m').  If 
%             select 'y', all averages will be alinged to the median of
%             the averages.  If you select 'a', all averages will be 
%             aligned to the average of the averages.  If you select 'r', 
%             all averages will be aligned to an externally provided 
%             reference spectrum.
% ref       = An externally provided reference spectrum that you would like
%             to align everything to (Required only if med = 'r').  
%
% OUTPUTS:
% out       = Output following alignment of averages.  
% fs        = Vector of frequency shifts (in Hz) used for alignment.
% phs       = Vector of phase shifts (in degrees) used for alignment.


function [out,fs,phs]=op_alignAverages(in,tmax,med,initPars)

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
    
    if nargin<3
        med='n'
        if nargin<2
            %if tmax is not specified, find the time at which the SNR
            %drops below 5
            disp('tmax not supplied.  Calculating tmax....');
            sig=abs(in.fids);
            noise=std(real(in.fids(ceil(0.75*end):end,:,:)),[]);
            noise=mean(mean(mean(noise,2),3),4);
            snr=sig/noise;
            
            for n=1:(numel(snr)/size(snr,1))
                N=find(snr(:,n)>5);
                tmax_est(n)=in.t(N(end));
            end
            tmax=median(tmax_est);
            disp(['tmax = ' num2str(tmax*1000) 'ms.']);
        end
    end
    
    if (strcmp(med,'r') || strcmp(med,'R'))
        if nargin<4
            error('ERROR:  If using the ''r'' option for input variable ''med'', then a 4th input argument must be provided');
        end
    else
        if nargin<4
            ref=struct();
        end
    end
    
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
            disp('Aligning all averages to the median of the averages.');
            base=op_median(in);
            base=[real(base.fids( in.t>=0 & in.t<tmax ,m));imag(base.fids( in.t>=0 & in.t<tmax ,m))];
            ind_min=0;
        elseif med=='a' || med=='a'
            disp('Aligning all averages to the average of the averages.');
            base=op_averaging(in);
            base=[real(base.fids( in.t>=0 & in.t<tmax ,m));imag(base.fids( in.t>=0 & in.t<tmax ,m))];
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
            base=[real(in.fids(in.t>=0 & in.t<tmax,ind_min,m));imag(in.fids(in.t>=0 & in.t<tmax,ind_min,m))];
            fids(:,ind_min,m)=in.fids(:,ind_min,m);
        elseif med=='r' || med=='R'
            disp('Aligning all averages to an externally provided reference spectrum.');
            base=ref;
            base=[real(base.fids( in.t>=0 & in.t<tmax ,m));imag(base.fids( in.t>=0 & in.t<tmax ,m))];
            ind_min=0;
        end
        for n=1:in.sz(in.dims.averages)
            if n~=ind_min
                parsGuess=parsFit;
                %disp(['fitting subspec number ' num2str(m) ' and average number ' num2str(n)]);
                parsFit=nlinfit(in.fids(in.t>=0 & in.t<tmax,n,m),base,@op_freqPhaseShiftComplexNest,parsGuess);
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
