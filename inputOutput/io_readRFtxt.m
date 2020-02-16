%io_readRFtxt.m
%Jamie Near, McGill University 2016.
%
% USAGE:
% [rf,info]=io_readRFtxt(filename)
% 
% DESCRIPTION:
% Read an RF pulse in basic .txt format into matlab.  The text file should
% contain two, three or four columns of data, with the first column 
% specifying the magnitude (arbitrary units) and the second column 
% specifying the phase (in degrees) of the RF waveform.  If a third column 
% exists, it will be the timestep waveform (arbitrary units).  If a fourth 
% column exists, it will be the gradient strength (in G/cm).  The resulting 
% RF matrix will have the same number of columns specifying phase, 
% magnitude, timestep, and gradient strength (if supplied).
% 
% INPUTS:
% filename   = filename of the .txt file to read in.  
%
% OUTPUTS:
% rf         = Input rf pulse waveform saved as a matlab array with 3 or 4
%               columns (phase, magnitude, duration, gradient (if supplied)).
% info       = Empty. Not required.  

function [rf,info]=io_readRFtxt(filename)

% For these basic text files, the data can be read simply using the load
% function:

RF=load(filename);
   
rf(:,1)=RF(:,2);
rf(:,2)=RF(:,1);
if size(RF,2)<3
    rf(:,3)=ones(length(RF(:,1)),1);
else
    rf(:,3)=RF(:,3);
    if size(RF,2)==4
        rf(:,4)=RF(:,4);
    end
end


info=[];


   