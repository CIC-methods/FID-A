% op_alignMPSubspecs.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% [out,fs,phs]=op_alignMPSubspecs(in,initPars);
% 
% DESCRIPTION:
% Apply spectral registration to align MEGA-PRESS subspectra prior to 
% subtraction.  This function is designed to minimize subtraction artefacts 
% from choline, and residual water.  This is intended to be used after averaging.
% 
% INPUTS:
% in        = Input data structure.
% initPars	= (Optional) Initial fit parameters [freq(Hz), phase(degrees)]. Default=[0,0];
%
% OUTPUTS:
% out       = Output following alignment of MEGA-PRESS subspectra.  
% fs        = Vector of frequency shifts (in Hz) used for alignment.
% phs       = Vector of phase shifts (in degrees) used for alignment.

function [out,fs,phs]=op_alignMPSubspecs(in,initPars)

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

base0=op_takesubspec(in,1);
base=[real(base0.specs);imag(base0.specs)];
begin=1;

%DO FITTING
n=1;
%disp(['fitting subspec number ' num2str(m) ' and average number ' num2str(n)]);
parsFit=nlinfit(in.fids(:,2),base,@op_freqPhaseShiftComplexNest,parsGuess);
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
        
        
        t=in.t;
        fid=input;
        
        shiftedFids=addphase(fid.*exp(-1i*t'*f*2*pi),p);
        shiftedSpecs=fftshift(ifft(shiftedFids,[],1),1);
        %shiftedSpecsWindow=shiftedSpecs(freqWindows);
        y=real(shiftedSpecs);
        y=[real(shiftedSpecs);imag(shiftedSpecs)];
        %y=real(fid.*exp(-1i*t'*f*2*pi));
        
    end

    function y=op_freqPhaseShiftNest(pars,input)
        f=pars(1);     %Frequency Shift [Hz]
        p=pars(2);     %Phase Shift [deg]
        
        
        t=in.t;
        fid=input;
        p=p+180;
        
        shiftedFids=addphase(fid.*exp(-1i*t'*f*2*pi),p);
        y=shiftedFids;
        %y=real(fid.*exp(-1i*t'*f*2*pi));
        
    end

end
