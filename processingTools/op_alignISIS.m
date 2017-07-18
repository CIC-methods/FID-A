% op_alignISIS.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% [out,fs,phs]=op_alignISIS(in,tmax,initPars);
% 
% DESCRIPTION:
% Apply spectral registration to align ISIS subspectra prior to 
% subtraction.  This is intended to be used prior to averaging, so that the
% alignment can be performed independently for each average.  
% 
% INPUTS:
% in        = Input data structure.
% tmax      = Maximum time (s) in time domain to use for alignment.
% initPars	= (Optional) Initial fit parameters [freq(Hz), phase(degrees)]. Default=[0,0];
%
% OUTPUTS:
% out       = Output following alignment of ISIS subspectra.  
% fs        = Vector of frequency shifts (in Hz) used for alignment.
% phs       = Vector of phase shifts (in degrees) used for alignment.

function [out,fs,phs]=op_alignISIS(in,tmax,initPars)

if ~in.flags.addedrcvrs
    error('ERROR:  I think it only makes sense to do this after you have combined the channels using op_addrcvrs.  ABORTING!!');
end

if in.dims.subSpecs==0
    error('ERROR:  Must have multiple subspectra.  ABORTING!!');
end

if nargin<3
    parsGuess=[0,0];
else
    parsGuess=initPars;
end

if in.dims.averages
    fs=zeros(in.sz(in.dims.averages),1);
    phs=zeros(in.sz(in.dims.averages),1);
else
    fs=0;
    phs=0;
end
fids=zeros(in.sz(in.dims.t),1);


disp('aligning all averages to the Average ISIS subtracted spectrum');
if in.dims.averages
    base0=op_averaging(op_combinesubspecs(in,'diff'));
else
    base0=op_combinesubspecs(in,'diff');
end
base=[real(base0.fids( base0.t>=0 & base0.t<tmax ));imag(base0.fids( base0.t>=0 & base0.t<tmax ))];
begin=1;
if in.dims.averages
    for n=begin:in.sz(in.dims.averages)
        %disp(['fitting subspec number ' num2str(m) ' and average number ' num2str(n)]);
        parsFit=nlinfit(squeeze(in.fids(in.t>=0 & in.t<tmax,n,:)),base,@op_freqPhaseShiftComplexNest,parsGuess);
        A=op_freqPhaseShiftNest(parsFit,in.fids(:,n,:));
        size(A);
        size(fids);
        fids(:,n,1)=A(:,1);
        fids(:,n,2)=A(:,2);
        fs(n)=parsFit(1);
        phs(n)=parsFit(2);
        %plot(in.ppm,fftshift(ifft(fids(:,1,m))),in.ppm,fftshift(ifft(fids(:,n,m))));   
    end
else
    n=1;
    %disp(['fitting subspec number ' num2str(m) ' and average number ' num2str(n)]);
    parsFit=nlinfit(squeeze(in.fids(in.t>=0 & in.t<tmax,:)),base,@op_freqPhaseShiftComplexNest,parsGuess);
    A=op_freqPhaseShiftNest(parsFit,in.fids(:,:));
    size(A);
    size(fids);
    fids(:,1)=A(:,1);
    fids(:,2)=A(:,2);
    fs(n)=parsFit(1);
    phs(n)=parsFit(2);
    %plot(in.ppm,fftshift(ifft(fids(:,1,m))),in.ppm,fftshift(ifft(fids(:,n,m))));
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
        fid=input(:,2);
        
        shifted=addphase(fid.*exp(-1i*t'*f*2*pi),p);
        subtracted=input(:,1)+shifted;
        subtracted=subtracted/2;
        
        y=[real(subtracted);imag(subtracted)];
        %y=real(fid.*exp(-1i*t'*f*2*pi));
        
    end

    function y=op_freqPhaseShiftNest(pars,input)
        f=pars(1);     %Frequency Shift [Hz]
        p=pars(2);     %Phase Shift [deg]
        
        
        dwelltime=in.dwelltime;
        t=0:dwelltime:(length(input)-1)*dwelltime;
        fid=input(:,2);
        
        shifted=addphase(fid.*exp(-1i*t'*f*2*pi),p);
        y(:,1)=input(:,1);
        y(:,2)=shifted;
        %y=real(fid.*exp(-1i*t'*f*2*pi));
        
    end

end
