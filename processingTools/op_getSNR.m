% op_getSNR.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% [SNR]=op_getSNR(in,NAAppmmin,NAAppmmax,noiseppmmin,noiseppmmax);
% 
% DESCRIPTION:
% Find the SNR of the NAA peak in a spectrum.
% 
% INPUTS:
% in             = input data in matlab structure format
% NAAppmmin      = min of frequncy range in which to search for NAA peak.
%                  (Optional.  Default = 1.8 ppm);
% NAAppmmax      = max of frequncy range in which to search for NAA peak.
%                  (Optional.  Default = 2.2 ppm);
% noiseppmmin    = min of frequency range in which to measure noise.
%                  (Optional.  Default = -2 ppm);
% noiseppmmax    = max of frequency range in which to measure noise.
%                  (Optional.  Default = 0 ppm);
%
% OUTPUTS:
% SNR            = Estimated SNR of the input spectrum.
% signal         = The measured raw signal intensity
% noisesd        = The measured noise standard deviation

function [SNR,signal,noisesd]=op_getSNR(in,NAAppmmin,NAAppmmax,noiseppmmin,noiseppmmax);


if nargin<5
    noiseppmmax=0;
    if nargin<4
        noiseppmmin=-2;
        if nargin<3
            NAAppmmax=2.2;
            if nargin<2
                NAAppmmin=1.8;
            end
        end
    end
end


%FIRST FIND THE NAA SIGNAL INTENSITY.  USE THE MAX PEAK HEIGHT OF THE 
%MAGNITUDE SPECTRUM INSIDE THE DESIRED SPECTRAL RANGE:
NAAwindow=in.specs(in.ppm>NAAppmmin & in.ppm<NAAppmmax);
ppmwindow=in.ppm(in.ppm>NAAppmmin & in.ppm<NAAppmmax);

maxNAA_index=find(abs(NAAwindow)==max(abs((NAAwindow))));
maxNAA=abs(NAAwindow(maxNAA_index))

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