%getLW.m
%Jamie Near, McGill University 2014.
%
%USAGE:
%[FWHM]=getLW(in,zpfactor,Refppmmin,Refppmmax);
%
%DESCRIPTION:
%Estimates the linewidth of a reference peak in the spectrum.  Two methods are
%used to estimate the linewidth:  1.  FWHM is measured by simply taking the
%full width at half max of the reference peak.  2.  The FWHM is measured by
%fitting the reference peak to a lorentzian lineshape and determine the FWHM of the
%best fit.  The output FWHM is given by the average of these two measures.
%
%INPUTS:
%in         = input spectrum in structure format.
%zpfactor   = zero-padding factor (used for method 1.)
%Refppmmin  = Min of frequency range (ppm) in which to search for reference peak.
%Refppmmax  = Max of frequency range to (ppm) in which search for reference peak

function [FWHM]=getLW(in,zpfactor,Refppmmin,Refppmmax);

%in=op_readlcmraw(filestring,'dat');
in=op_zeropad(in,zpfactor);

%FIRST FIND FWHM USING TWO METHODS:

%METHOD 1:  ACTUALLY MEAUSURE FWHM OF WATER PEAK
Refwindow=in.specs(in.ppm>Refppmmin & in.ppm<Refppmmax);
ppmwindow=in.ppm(in.ppm>Refppmmin & in.ppm<Refppmmax);

maxRef_index=find(abs(real(Refwindow))==max(abs(real((Refwindow)))));
maxRef=real(Refwindow(maxRef_index));

plot(ppmwindow,abs(real(Refwindow)),'.');

gtHalfMax=find(abs(real(Refwindow)) >= 0.5*abs(maxRef));

FWHM1=ppmwindow(gtHalfMax(1)) - ppmwindow(gtHalfMax(end));
FWHM1=FWHM1*(42.577*in.Bo);  %Assumes proton.


%METHOD 2:  FIT WATER PEAK TO DETERMINE FWHM PARAM
sat='n'
waterFreq=ppmwindow(maxRef_index);
while sat=='n'
    parsGuess=zeros(1,5);
    parsGuess(1)=maxRef; %AMPLITUDE
    parsGuess(2)=5/(42.577*in.Bo); %FWHM.  Assumes Proton
    parsGuess(3)=waterFreq; %FREQUENCY
    parsGuess(4)=0; %Baseline Offset
    parsGuess(5)=0; %Phase
    
    yGuess=op_lorentz(parsGuess,ppmwindow);
    parsFit=nlinfit(ppmwindow,real(Refwindow'),@op_lorentz,parsGuess);
    yFit=op_lorentz(parsFit,ppmwindow);
    
    figure;
    plot(ppmwindow,Refwindow,'.',ppmwindow,yGuess,':',ppmwindow,yFit);
    legend('data','guess','fit');
    
    sat=input('are you satisfied with fit? y/n [y] ','s');
    if isempty(sat)
        sat='y';
    end
    if sat=='n';
        waterFreq=input('input new water frequency guess: ');
    end

end


FWHM2=abs(parsFit(2));
FWHM2=FWHM2*(42.577*in.Bo);  %Assumes Proton.

FWHM=mean([FWHM1 FWHM2]);  