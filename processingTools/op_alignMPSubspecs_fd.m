% op_alignMPSubspecs_fd.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% [out,fs,phs]=op_alignMPSubspecs_fd(in,initPars);
% 
% DESCRIPTION:
% Apply spectral registration to align MEGA-PRESS subspectra prior to 
% subtraction.  This function is designed to minimize subtraction artefacts 
% from choline, and residual water.  This is intended to be used after averaging.
% 
% INPUTS:
% in        = Input data structure.
% initPars	= (Optional) Initial fit parameters [freq(Hz), phase(degrees)]. Default=[0,0];

function [out,fs,phs]=op_alignMPSubspecs_fd(in,initPars)

if ~in.flags.addedrcvrs
    error('ERROR:  I think it only makes sense to do this after you have combined the channels using op_addrcvrs.  ABORTING!!');
end

if in.dims.subSpecs==0
    error('ERROR:  Must have multiple subspectra.  ABORTING!!');
end

if in.dims.averages
   error('ERROR:  Signal averaging must be performed before this step.  ABORTING!!');
end

if nargin<3
    parsGuess=[0,0];
else
    parsGuess=initPars;
end

fs=0;
phs=0;
fids=zeros(in.sz(in.dims.t),1);


disp('aligning the MEGA-PRESS edit-ON sub-spectrum to the edit-OFF sub-spectrum');

%NOTE:  Alignment is not performed in the regions of the spectrum where
%the difference spectrum is large (ie.  near NAA, 3ppm GABA, or 3.7ppm
%Glx).
minfrq1=2.46;
maxfrq1=2.9;
minfrq2=3.1;
maxfrq2=3.6;
minfrq3=3.85;
maxfrq3=4.1;
minfrq4=4.6;
maxfrq4=5.2;

%freqWindow1=in.ppm>minfrq1 & in.ppm<maxfrq1; 
%freqWindow2=in.ppm>minfrq2 & in.ppm<maxfrq2; 
%freqWindow3=in.ppm>minfrq3 & in.ppm<maxfrq3; 
%freqWindow4=in.ppm>minfrq4 & in.ppm<maxfrq4;

%freqWindows=freqWindow1+freqWindow2+freqWindow3+freqWindow4;
%freqWindows=logical(freqWindows);

base0=op_freqrange(op_takesubspec(in,1),minfrq3,maxfrq4);
start=op_freqrange(op_takesubspec(in,2),minfrq3,maxfrq4);

%base=[real(base0.specs( freqWindows ));imag(base0.specs( freqWindows ))];
base=[real(base0.fids);imag(base0.fids)];
begin=1;

%size(in.ppm(freqWindows))
%size(base)
%plot([in.ppm(freqWindows) in.ppm(freqWindows)],base,'.');

%DO FITTING
n=1;
%disp(['fitting subspec number ' num2str(m) ' and average number ' num2str(n)]);
parsFit=nlinfit(start.fids,base,@op_freqPhaseShiftComplexNest,parsGuess);
A=op_freqPhaseShiftNest(parsFit,in.fids(:,2));
size(A);
size(fids);
fids(:,1)=in.fids(:,1);
fids(:,2)=A;
fs=parsFit(1);
phs=parsFit(2)+180;
%plot(in.ppm,fftshift(ifft(fids(:,1,m))),in.ppm,fftshift(ifft(fids(:,n,m))));

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
        
        
        dwelltime=start.dwelltime;
        t=0:dwelltime:(length(input)-1)*dwelltime;
        fid=input;
        
        shiftedFids=addphase(fid.*exp(-1i*t'*f*2*pi),p);
%         shiftedSpecs=fftshift(ifft(shiftedFids,[],1),1);
%         shiftedSpecsWindow=shiftedSpecs(ppm>minfrq3 & ppm<maxfrq4);
        %y=real(shiftedSpecs);
        y=[real(shiftedFids);imag(shiftedFids)];
        %y=real(fid.*exp(-1i*t'*f*2*pi));
        
    end

    function y=op_freqPhaseShiftNest(pars,input)
        f=pars(1);     %Frequency Shift [Hz]
        p=pars(2);     %Phase Shift [deg]
        
        
        dwelltime=start.dwelltime;
        t=0:dwelltime:(length(input)-1)*dwelltime;
        fid=input;
        p=p+180;
        
        shiftedFids=addphase(fid.*exp(-1i*t'*f*2*pi),p);
        y=shiftedFids;
        %y=real(fid.*exp(-1i*t'*f*2*pi));
        
    end

end
