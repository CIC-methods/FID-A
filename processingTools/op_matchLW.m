% op_matchLW.m
% Jamie Near, Floris van Vugt, McGill University 2018.
% 
% USAGE:
% [out,lb,k]=op_matchLW(in,ref,ppmmin,ppmmax,tmax,mode);
% 
% DESCRIPTION:
% Line-broaden the input spectrum (in) to match the linewidth of the
% reference spectrum (ref).  A spectral similarity metric, similar to the
% one used in 'op_rmbadaverages.m' is used as the stopping criterion.  The
% script allows the choice of a spectral range (ppmmin and ppmmax) over
% which to compute the similarity measure (it is done in the time
% domain).  Finally, a 'mode' parameter allows the choice of different
% apodization functions (lorentzian, gaussian or voigt).  
% 
% INPUTS: 
% in	 = input (spectrum to be registered to the reference)
% ref	 = reference (spectrum that the input is to be broadened to match).
% ppmmin = Minimum frequency for spectral alignment (ppm).
% ppmmax = Maximum frequency for spectral alignment (ppm).
% tmax   = Maximum time in the time domain over which to compute
%          similarity (s).
% mode   = Apodization function to use (optional)
%                   'l' - lorentzian (default)
%                   'g' - gaussian
%                   'lg' - mixture of lorentizan and gaussian
%
% OUTPUTS:
% out   = Output following broadening of input (in1) to reference.  
% lb    = Line broadening factor (in Hz) that was used.
% k     = Proportion of lorentzian to gaussian broadening that was used. 

function [out,lb,k]=op_matchLW(in,ref,ppmmin,ppmmax,tmax,mode);

if ~in.flags.addedrcvrs || ~ref.flags.addedrcvrs
    error('ERROR:  I think it only makes sense to do this after you have combined the channels using op_addrcvrs.  ABORTING!!');
end

if ~in.flags.averaged || ~ref.flags.averaged
    error('ERROR:  I think it only makes sense to do this after you have combined the averages using op_averaging.  ABORTING!!');
end

if in.flags.isISIS || ref.flags.isISIS
    error('ERROR:  I think it only makes sense to do this after you have performed op_ISIScombine.  ABORTING!!');
end

% if ~in1.flags.subtracted  || ~in.flags.subtracted
%     error('ERROR:  I think it only makes sense to do this after you have combined the subspectra using op_combinesubspecs.  ABORTING!!');
% end

if nargin<6
    mode='l';
end

switch mode
    case 'l'
        parsFit=0.5;
    case 'g'
        parsFit=0.5;
    case 'lg'
        parsFit=[0.5 0.5];
    case 'gl'
        mode='lg';
        parsFit=[0.5 0.5];
    otherwise
        error('ERROR:mode unrecognized.  Use "1" or "2".  ABORTING!');
end


baserange=op_freqrange(ref,ppmmin,ppmmax);
base=[real(baserange.fids(baserange.t>=0 & baserange.t<tmax));imag(baserange.fids(baserange.t>=0 & baserange.t<tmax))];
%plot(in1.t,in1.fids,in.t,in.fids);

parsGuess=parsFit;
%disp(['fitting subspec number ' num2str(m) ' and average number ' num2str(n)]);
switch mode
    case 'l'
        startrange=op_freqrange(in,ppmmin,ppmmax);
        start=startrange.fids(startrange.t>=0 & startrange.t<tmax);
        parsFit=nlinfit(start,base,@op_lorFilterComplexRangeNest,parsGuess);
        fids=op_lorFilterNest(parsFit,in.fids);
    case 'g'
        startrange=op_freqrange(in,ppmmin,ppmmax);
        start=startrange.fids(startrange.t>=0 & startrange.t<tmax);
        parsFit=nlinfit(start,base,@op_gaussFilterComplexRangeNest,parsGuess);
        fids=op_gaussFilterNest(parsFit,in.fids);
    case 'lg'
        startrange=op_freqrange(in,ppmmin,ppmmax);
        start=startrange.fids(startrange.t>=0 & startrange.t<tmax);
        parsFit=nlinfit(start,base,@op_lorGaussFilterComplexRangeNest,parsGuess);
        fids=op_lorGaussFilterNest(parsFit,in.fids);
        k=parsFit(2)
end
lb=parsFit(1);
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
out.flags=in.flags;


    function y=op_lorFilterComplexRangeNest(pars,input);
        lb=pars(1);     %Line broadening [Hz]
        %k=pars(2);     %proportion of Lor/(Lor+Gaus)
        
        
        dwelltime=startrange.dwelltime;
        t=[0:dwelltime:(length(input)-1)*dwelltime];
        fid=input(:);
        t2=1/(pi*lb);
        lor=exp(-t/t2);
        
        broadened=fid.*lor';
        
        y=[real(broadened);imag(broadened)];
        %y=real(fid.*exp(-1i*t'*f*2*pi));
        
    end

    function y=op_lorFilterNest(pars,input);
        lb=pars(1);     %line broadening [Hz]
        %k=pars(2);     %proportion of Lor/(Lor+Gaus)
        
        
        dwelltime=in.dwelltime;
        t=[0:dwelltime:(length(input)-1)*dwelltime];
        fid=input(:);
        t2=1/(pi*lb);
        lorFilt=exp(-t/t2);
        
        y=fid.*lorFilt';
    end

    function y=op_gaussFilterComplexRangeNest(pars,input);
        lb=pars(1);     %Line broadening [Hz]
        %k=pars(2);     %proportion of Lor/(Lor+Gaus)
        
        
        dwelltime=startrange.dwelltime;
        t=[0:dwelltime:(length(input)-1)*dwelltime];
        fid=input(:);
        thalf=log(0.5)/(pi*0.5*lb);
        sigma=sqrt((thalf^2)/(-2*log(0.5)));
        gausFilt=exp(-(t.^2)/(2*(sigma^2)));
        
        broadened=fid.*gausFilt';
        
        y=[real(broadened);imag(broadened)];
        %y=real(fid.*exp(-1i*t'*f*2*pi));
        
    end

    function y=op_gaussFilterNest(pars,input);
        lb=pars(1);     %line broadening [Hz]
        %k=pars(2);     %proportion of Lor/(Lor+Gaus)
        
        
        dwelltime=in.dwelltime;
        t=[0:dwelltime:(length(input)-1)*dwelltime];
        fid=input(:);
        thalf=log(0.5)/(pi*0.5*lb);
        sigma=sqrt((thalf^2)/(-2*log(0.5)));
        gausFilt=exp(-(t.^2)/(2*(sigma^2)));
        
        y=fid.*gausFilt';
    end

    function y=op_lorGaussFilterComplexRangeNest(pars,input);
        lb=pars(1);    %Line broadening [Hz]
        k=pars(2);     %proportion of Lor/(Lor+Gaus)
        
        
        dwelltime=startrange.dwelltime;
        t=[0:dwelltime:(length(input)-1)*dwelltime];
        fid=input(:);
        t2=1/(pi*lb);
        lorFilt=exp(-t/t2);
        thalf=log(0.5)/(pi*0.5*lb);
        sigma=sqrt((thalf^2)/(-2*log(0.5)));
        gausFilt=exp(-(t.^2)/(2*(sigma^2)));
        lorGausFilt=(k.*lorFilt) + ((1-k)*gausFilt);
        
        broadened=fid.*lorGausFilt';
        
        y=[real(broadened);imag(broadened)];
        %y=real(fid.*exp(-1i*t'*f*2*pi));
        
    end

    function y=op_lorGaussFilterNest(pars,input);
        lb=pars(1);     %line broadening [Hz]
        k=pars(2);      %proportion of Lor/(Lor+Gaus)
        
        
        dwelltime=in.dwelltime;
        t=[0:dwelltime:(length(input)-1)*dwelltime];
        fid=input(:);
        t2=1/(pi*lb);
        lorFilt=exp(-t/t2);
        thalf=log(0.5)/(pi*0.5*lb);
        sigma=sqrt((thalf^2)/(-2*log(0.5)));
        gausFilt=exp(-(t.^2)/(2*(sigma^2)));
        lorGausFilt=(k.*lorFilt) + ((1-k)*gausFilt);
        
        y=fid.*lorGausFilt';
    end

    
end
