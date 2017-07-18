% op_freqAlignAverages.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% [out,fs]=op_freqAlignAverages(in,tmax,avg,initPars);
% 
% DESCRIPTION:
% Perform spectral registration in the time domain using only frequency 
% adjustment (no phase adjustment).
% 
% INPUTS:
% in         = Input data structure.
% tmax       = Maximum time (s) in time domain to use for alignment.
% avg        = Align averages to the average of the averages? (y or n)
% initPars	 = (Optional) Initial fit parameters [freq(Hz), phase(degrees)]. Default=[0,0];
%
% OUTPUTS:
% out        = Output following alignment of averages.  
% fs         = Vector of frequencies (in Hz) used for alignment.

function [out,fs]=op_freqAlignAverages(in,tmax,avg,initPars);

if ~in.flags.addedrcvrs
    error('ERROR:  I think it only makes sense to do this after you have combined the channels using op_addrcvrs.  ABORTING!!');
end

if nargin<4
    parsFit=0;
else
    parsFit=initPars;
end

if in.dims.subSpecs==0
    B=1;
else
    B=in.sz(in.dims.subSpecs);
end

fs=zeros(in.sz(in.dims.averages),B);
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
        parsFit=nlinfit(in.fids(in.t>=0 & in.t<tmax,n,m),base,@op_freqShiftComplexNest,parsGuess);
        fids(:,n,m)=op_freqShiftNest(parsFit,in.fids(:,n,m));
        fs(n,m)=parsFit(1);
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


    function y=op_freqShiftComplexNest(pars,input);
        f=pars(1);     %Frequency Shift [Hz]
        
        dwelltime=in.dwelltime;
        t=[0:dwelltime:(length(input)-1)*dwelltime];
        fid=input(:);
        
        shifted=fid.*exp(1i*t'*f*2*pi);
        
        y=[real(shifted);imag(shifted)];
        %y=real(fid.*exp(1i*t'*f*2*pi));
        
    end

    function y=op_freqShiftNest(pars,input);
        f=pars(1);     %Frequency Shift [Hz]
        
        
        dwelltime=in.dwelltime;
        t=[0:dwelltime:(length(input)-1)*dwelltime];
        fid=input(:);
        
        y=fid.*exp(1i*t'*f*2*pi);
        %y=real(fid.*exp(1i*t'*f*2*pi));
        
    end

end
