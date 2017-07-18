% op_alignScans_fd.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% [out,ph,frq]=op_alignScans_fd(in,in1,fmin,fmax,tmax,mode);
% 
% DESCRIPTION:
% Use spectral registration on a limited frequnecy range to align two 
% separate MRS datasets (align in to in1);
% 
% INPUTS:
% in	= input (spectrum to be registered to the base)
% in1	= base (spectrum that the input is to be registered to).
% fmin	= Minimum frequency for spectral alignment (ppm).
% fmax	= Maximum frequency for spectral alignment (ppm).
% tmax	= Maximum time (s) in time domain to use for registration.
% mode  = (optional)'f' - Frequency align only
%                  'p' - Phase align only
%                  'fp or pf' - Frequency and phase align (default)
%
% OUTPUTS:
% out   = Output following alignment of input (in1) to base.  
% ph    = Phase shift (in degrees) used for alignment.
% frq   = Frequency shift (in Hz) used for alignment.

function [out,ph,frq]=op_alignScans_fd(in,in1,fmin,fmax,tmax,mode);

if ~in1.flags.addedrcvrs || ~in.flags.addedrcvrs
    error('ERROR:  I think it only makes sense to do this after you have combined the channels using op_addrcvrs.  ABORTING!!');
end

if ~in1.flags.averaged || ~in.flags.averaged
    error('ERROR:  I think it only makes sense to do this after you have combined the averages using op_averaging.  ABORTING!!');
end

if in1.flags.isISIS || in.flags.isISIS
    error('ERROR:  I think it only makes sense to do this after you have performed op_ISIScombine.  ABORTING!!');
end

% if ~in1.flags.subtracted  || ~in.flags.subtracted
%     error('ERROR:  I think it only makes sense to do this after you have combined the subspectra using op_combinesubspecs.  ABORTING!!');
% end

if nargin<6
    mode='fp';
end

switch mode
    case 'f'
        parsFit=0;
    case 'p'
        parsFit=0;
    case 'fp'
        parsFit=[0 0];
    case 'pf'
        mode='fp'
        parsFit=[0 0];
    otherwise
        error('ERROR:mode unrecognized.  Use "1" or "2".  ABORTING!');
end


baserange=op_freqrange(in1,fmin,fmax);
base=[real(baserange.fids(baserange.t>=0 & baserange.t<tmax));imag(baserange.fids(baserange.t>=0 & baserange.t<tmax))];
%plot(in1.t,in1.fids,in.t,in.fids);

parsGuess=parsFit;
%disp(['fitting subspec number ' num2str(m) ' and average number ' num2str(n)]);
switch mode
    case 'f'
        startrange=op_freqrange(in,fmin,fmax);
        start=startrange.fids(startrange.t>=0 & startrange.t<tmax);
        parsFit=nlinfit(start,base,@op_freqShiftComplexRangeNest,parsGuess);
        fids=op_freqShiftNest(parsFit,in.fids);
    case 'p'
        startrange=op_freqrange(in,fmin,fmax);
        start=startrange.fids(startrange.t>=0 & startrange.t<tmax);
        parsFit=nlinfit(start,base,@op_phaseShiftComplexRangeNest,parsGuess);
        fids=op_phaseShiftNest(parsFit,in.fids);
    case 'fp'
        startrange=op_freqrange(in,fmin,fmax);
        start=startrange.fids(startrange.t>=0 & startrange.t<tmax);
        parsFit=nlinfit(start,base,@op_freqPhaseShiftComplexRangeNest,parsGuess);
        fids=op_freqPhaseShiftNest(parsFit,in.fids);
end
ph=parsFit(2);
frq=parsFit(1);
%figure
%plot(in1.t,in1.fids,in.t,fids);


%re-calculate Specs using fft
specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);

%plot(in1.ppm,combinedSpecs);

%FILLING IN DATA STRUCTURES
out=in;
out.fids=fids;
out.specs=specs;

%FILLING IN THE FLAGS
out.flags=in1.flags;


    function y=op_freqPhaseShiftComplexRangeNest(pars,input);
        f=pars(1);     %Frequency Shift [Hz]
        p=pars(2);     %Phase Shift [deg]
        
        
        dwelltime=startrange.dwelltime;
        t=[0:dwelltime:(length(input)-1)*dwelltime];
        fid=input(:);
        
        shifted=addphase(fid.*exp(-1i*t'*f*2*pi),p);
        
        y=[real(shifted);imag(shifted)];
        %y=real(fid.*exp(-1i*t'*f*2*pi));
        
    end

    function y=op_freqPhaseShiftNest(pars,input);
        f=pars(1);     %Frequency Shift [Hz]
        p=pars(2);     %Phase Shift [deg]
        
        
        dwelltime=in.dwelltime;
        t=[0:dwelltime:(length(input)-1)*dwelltime];
        fid=input(:);
        
        y=addphase(fid.*exp(-1i*t'*f*2*pi),p);
        %y=real(fid.*exp(-1i*t'*f*2*pi));
    end

    function y=op_freqShiftComplexRangeNest(pars,input);
        f=pars(1);     %Frequency Shift [Hz]
                
        dwelltime=startrange.dwelltime;
        t=[0:dwelltime:(length(input)-1)*dwelltime];
        fid=input(:);
        
        shifted=fid.*exp(-1i*t'*f*2*pi);
        
        y=[real(shifted);imag(shifted)];
        %y=real(fid.*exp(-1i*t'*f*2*pi));
        
    end

    function y=op_freqShiftNest(pars,input);
        f=pars(1);     %Frequency Shift [Hz]
               
        dwelltime=in.dwelltime;
        t=[0:dwelltime:(length(input)-1)*dwelltime];
        fid=input(:);
        
        y=fid.*exp(-1i*t'*f*2*pi);
        %y=real(fid.*exp(-1i*t'*f*2*pi));
    end
    function y=op_phaseShiftComplexRangeNest(pars,input);
        p=pars(1);     %Phase Shift [deg]
        
        
        dwelltime=startrange.dwelltime;
        t=[0:dwelltime:(length(input)-1)*dwelltime];
        fid=input(:);
        
        shifted=addphase(fid,p);
        
        y=[real(shifted);imag(shifted)];
        %y=real(fid.*exp(-1i*t'*f*2*pi));
        
    end

    function y=op_phaseShiftNest(pars,input);
        p=pars(1);     %Phase Shift [deg]
        
        
        dwelltime=in.dwelltime;
        t=[0:dwelltime:(length(input)-1)*dwelltime];
        fid=input(:);
        
        y=addphase(fid,p);
        %y=real(fid.*exp(-1i*t'*f*2*pi));
    end
    
end
