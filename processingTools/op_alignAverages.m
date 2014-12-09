%op_alignAverages.m
%Jamie Near, McGill University 2014.
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
% INPUTS:
% in        = Input data structure.
% tmax      = Maximum time (s) in time domain to use for alignment.
% avg       = Align averages to the average of the averages? (y or n)
% initPars	= (Optional) Initial fit parameters [freq(Hz), phase(degrees)]. Default=[0,0];

function [out,fs,phs]=op_alignAverages(in,tmax,avg,initPars)

if ~in.flags.addedrcvrs
    error('ERROR:  I think it only makes sense to do this after you have combined the channels using op_addrcvrs.  ABORTING!!');
end

if nargin<4
    parsFit=[0,0];
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
        base=[real(base.fids( in.t>=0 & in.t<tmax ,m));imag(base.fids( in.t>=0 & in.t<tmax ,m))]/in.sz(in.dims.averages);
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
