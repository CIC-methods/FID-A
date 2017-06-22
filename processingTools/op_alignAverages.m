% op_alignAverages.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% [out,fs,phs]=op_alignAverages(in,tmax,avg,initPars);
% 
% DESCRIPTION:
% Perform spectral registration in the time domain to correct frequency and
% phase drifts.  As described in Near et al.  Frequency and phase drift 
% correction of magnetic resonance spectroscopy data by spectral 
% registration in the time domain. Magn Reson Med. 2014 Jan 16. 
% doi: 10.1002/mrm.25094. [Epub ahead of print]
%
% June 15th 2017:  Made the tmax and avg arguments optional.  If tmax is 
% not specified, the value is determined automatically by finding the time
% at which the SNR of the FID drops permanently below 5.  This idea
% was suggested by Mark Mikkelsen.  Thanks Mark!!
% 
% INPUTS:
% in        = Input data structure.
% tmax      = Maximum time (s) in time domain to use for alignment.
%             (Optional.  Default is the time at which SNR drops below 3)
% avg       = Align averages to the average of the averages? (y or n)
%             (Optional.  Default = 'n')
% initPars	= (Optional) Initial fit parameters [freq(Hz), phase(degrees)]. Default=[0,0];

function [out,fs,phs]=op_alignAverages(in,tmax,avg,initPars)

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

    if nargin<4
        parsFit=[0,0];
        if nargin<3
            avg='n'
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
    else
        parsFit=initPars;
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
        if avg=='y' || avg=='Y'
            disp('aligning all averages to the Average of the averages');
            base=op_averaging(in);
            base=[real(base.fids( in.t>=0 & in.t<tmax ,m));imag(base.fids( in.t>=0 & in.t<tmax ,m))];
            begin=1;
        else
            disp('aligning all averages to the first average');
            base=[real(in.fids(in.t>=0 & in.t<tmax,1,m));imag(in.fids(in.t>=0 & in.t<tmax,1,m))];
            begin=2;
            fids(:,1,m)=in.fids(:,1,m);
        end
        for n=begin:in.sz(in.dims.averages)
            parsGuess=parsFit;
            %disp(['fitting subspec number ' num2str(m) ' and average number ' num2str(n)]);
            parsFit=nlinfit(in.fids(in.t>=0 & in.t<tmax,n,m),base,@op_freqPhaseShiftComplexNest,parsGuess);
            fids(:,n,m)=op_freqPhaseShiftNest(parsFit,in.fids(:,n,m));
            fs(n,m)=parsFit(1);
            phs(n,m)=parsFit(2);
            %plot(in.ppm,fftshift(ifft(fids(:,1,m))),in.ppm,fftshift(ifft(fids(:,n,m))));
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
