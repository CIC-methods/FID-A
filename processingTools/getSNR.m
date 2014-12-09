%getSNR.m
%Jamie Near, McGill University 2014.
%
%USAGE:
%[SNR]=getSNR(in,NAAppmmin,NAAppmmax,noiseppmmin,noiseppmmax);
%
%DESCRIPTION:
%Find the SNR of the NAA peak in a spectrum.
%
%INPUTS:
%in             = input data in matlab structure format
%NAAppmmin      = min of frequncy range in which to search for NAA peak.
%NAAppmmax      = max of frequncy range in which to search for NAA peak.
%noiseppmmin    = min of frequency range in which to measure noise.
%noiseppmmax    = max of frequency range in which to measure noise.

function [SNR]=getSNR(in,NAAppmmin,NAAppmmax,noiseppmmin,noiseppmmax);


%FIRST FIND THE NAA PEAK HEIGHT:
NAAwindow=in.specs(in.ppm>NAAppmmin & in.ppm<NAAppmmax);
ppmwindow=in.ppm(in.ppm>NAAppmmin & in.ppm<NAAppmmax);

maxNAA_index=find(abs(real(NAAwindow))==max(abs(real((NAAwindow)))));
maxNAA=real(NAAwindow(maxNAA_index))

figure;
plot(ppmwindow,abs(real(NAAwindow)));

figure
plot(in.ppm,in.specs);
%noiseppmmin=input('input lower ppm limit for noise: ');
%noiseppmmax=input('input upper ppm limit for noise: ');

%NOW FIND THE STANDARD DEVIATION OF THE NOISE:
noisewindow=in.specs(in.ppm>noiseppmmin & in.ppm<noiseppmmax);
ppmwindow2=in.ppm(in.ppm>noiseppmmin & in.ppm<noiseppmmax)';

P=polyfit(ppmwindow2,noisewindow,2);
noise=noisewindow-polyval(P,ppmwindow2);
figure
plot(ppmwindow2,noisewindow,...
    ppmwindow2,polyval(P,ppmwindow2),...
    ppmwindow2,noise);

signal=(maxNAA-mean(real(noisewindow))) %Removes DC offset

noisesd=std(real(noise))

%SNR=maxNAA/noisesd
SNR=signal/noisesd;