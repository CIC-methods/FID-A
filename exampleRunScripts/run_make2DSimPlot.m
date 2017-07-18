% run_make2DSimPlot.m
% Jamie Near, McGill University 2015.
% 
% USAGE:
% []=run_make2DSimPlot(in,ppmmin,ppmmax,plotDiff)
% 
% DESCRIPTION:
% This function takes the output of a spatially resolved simulation, and
% plots the array of spectra on a single figure.  The input should be a
% cell array where the grid of elements of the cell array are simulated 
% spectra from a corresponding grid of spatial positions in the spatially 
% resolved simulation.  Each element of the cell array is also in FID-A 
% data struture format.  By including the optional input arguement ppmmin 
% and ppmmax, only a the corresponding range of each spectrum will be 
% plotted.
% 
% INPUTS:
% in          = input cell array of simulated spectra from a spatially resolved simulation
% ppmmin      = lower limit of ppm range to plot [ppm]
% ppmmax      = upper limit of ppm range to plot [ppm]

function []=run_make2DSimPlot(in,ppmmin,ppmmax)

if nargin<3
    if nargin<2
        ppmmin=min(in{1}{1}.ppm);
        ppmmax=max(in{1}{1}.ppm);
    else
        ppmmax=max(in{1}{1}.ppm);
    end
end
   
%initialize the figure and axes.  The axes are in arbitrary units, and will 
%span from -1 to the size of the input cell array + 1;
X=size(in,1);
Y=size(in,2);
figure;
hold;
xlim([-1 X+1]);
ylim([-1 Y+1]);

%Find out the current range of the ppm scale so that each ppm scale can be
%transformed down to a range of 0.8 (ie.  the first column will be plotted
%from -0.4 to 0.4 on the x-axis).  Here, we assume that the ppm range in
%each one of the cell array elements is the same.  
ppmrange=ppmmax-ppmmin;
scalefactorX=0.8/ppmrange;

%Estimate out the current range of amplitudes of the data so that the data
%can be plotted all on the same plotting array.
if X<2 || Y<2
    temp1=op_freqrange(in{ceil(X/2)}{ceil(Y/2)},ppmmin,ppmmax);
    temp2=op_freqrange(in{ceil(X/2)}{ceil(Y/2)},ppmmin,ppmmax);
else
    temp1=op_freqrange(in{floor(X/2)}{floor(Y/2)},ppmmin,ppmmax);
    temp2=op_freqrange(in{ceil(X/2)}{ceil(Y/2)},ppmmin,ppmmax);
end
tempppm=temp1.ppm;
yrange1=max(real(temp1.specs))-min(real(temp1.specs));
yrange2=max(real(temp2.specs))-min(real(temp2.specs));
yrange=max([yrange1 yrange2]);
scalefactorY=0.8/yrange;

for x=1:size(in,1)
    for y=1:size(in,2)
        %first scale the ppm scale so that the range is correct;
        ppm=tempppm*scalefactorX;
        %now zero it;
        ppm=ppm-min(ppm);
        %now shift it to the correct x-position;
        ppm=ppm+(x-1.4);
        
        %now generate the spectrum across the appropriate range
        tempSpec=op_freqrange(in{x}{y},ppmmin,ppmmax);
        
        %now do amplitude scaling on the spectrum
        tempSpec=op_ampScale(tempSpec,scalefactorY);
        
        %Now start plotting
        plot(ppm,tempSpec.specs+(y-1));
    end
end
