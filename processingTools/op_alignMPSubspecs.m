% op_alignMPSubspecs.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% [out,fs,phs]=op_alignMPSubspecs(in,mode,initPars,ppmWeights);
% 
% DESCRIPTION:
% Apply spectral registration to align MEGA-PRESS subspectra prior to 
% subtraction.  By default this function is designed to minimize subtraction artefacts 
% accross the entrie spectrum. This can be adjust by suppling weights to focus on
% specific parts of the spectrum.
% 
% This is intended to be used after averaging.
%
% August 2021:  In the default operation of this function, the edit-on and
% edit-off sub-spectra will end up 180 degrees out of phase.  However, I
% have now added an option to return the output with the sub-spectra
% in-phase.  This choice is controlled by the new input parameter (mode).  
% 
% INPUTS:
% in        = Input data structure.
% mode      = (Optional) Mode 'o' will return the output such that the 
%             sub-spectra are 180 degrees out of phase.  Mode 'i' will
%             return the output such that the sub-spectra are in phase.
%             Default = 'o'.  
% initPars	= (Optional) Initial fit parameters [freq(Hz), phase(degrees)].
%             Default=[0,0];
% ppmWeights= (Optional) Vector (size(in.ppm)) specifing positive weights
%             for nlinfit.
%             Default=ones(size(in.ppm));
%
% OUTPUTS:
% out       = Output following alignment of MEGA-PRESS subspectra.  
% fs        = Vector of frequency shifts (in Hz) used for alignment.
% phs       = Vector of phase shifts (in degrees) used for alignment.

function [out,fs,phs]=op_alignMPSubspecs(in,mode,initPars,ppmWeights)

if ~in.flags.addedrcvrs
    error('ERROR:  I think it only makes sense to do this after you have combined the channels using op_addrcvrs.  ABORTING!!');
end

if in.dims.subSpecs==0
    error('ERROR:  Must have multiple subspectra.  ABORTING!!');
end

if in.dims.averages
   error('ERROR:  Signal averaging must be performed before this step.  ABORTING!!');
end

if nargin<4
   ppmWeights = ones(size(in.ppm));
   if nargin<3
       parsGuess=[0,0];
       if nargin<2
           mode='o';
       end
   end
else
   parsGuess=initPars;
end

% Normalize weights 
ppmWeights = ppmWeights / sum(ppmWeights);

% do basic error Check on ppmWeights
if (~isequal(size(ppmWeights), size(in.ppm)) || ...
    min(ppmWeights) <= 0)
    error('ppmWeights must be A vector of real positive weights the same size as in.ppm');
end


%Set the phase shift depending on whether in-phase or out-of-phase
%subs-spectra are requested.  
if strcmp(mode,'o') || strcmp(mode,'O')
    phShift=180;
elseif strcmp(mode,'i') || strcmp(mode,'I')
    phShift=0;
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
parsFit=nlinfit(in.fids(:,2),base,@op_freqPhaseShiftComplexNest,parsGuess,...
                'weights', [ppmWeights ppmWeights]');
A=op_freqPhaseShiftNest(parsFit,in.fids(:,2),phShift);
size(A);
size(fids);
fids(:,1)=in.fids(:,1);
fids(:,2)=A;
fs=parsFit(1);
phs=parsFit(2)+phShift;
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

    function y=op_freqPhaseShiftNest(pars,input,phadd)
        f=pars(1);     %Frequency Shift [Hz]
        p=pars(2);     %Phase Shift [deg]
        
        
        t=in.t;
        fid=input;
        p=p+phadd;
        
        shiftedFids=addphase(fid.*exp(-1i*t'*f*2*pi),p);
        y=shiftedFids;
        %y=real(fid.*exp(-1i*t'*f*2*pi));
        
    end

end
