%io_readRFtxt.m
%Jamie Near, McGill University 2016.
%
% USAGE:
% [rf,info]=io_readRFtxt(filename,colOrder)
% 
% DESCRIPTION:
% Read an RF pulse in basic .txt format into matlab.  The text file should
% contain two, three or four columns of data, specifying the magnitude 
% (arbitrary units) and phase (in degrees) of the RF waveform, and 
% optionally, the timestep waveform (arbitrary units) and gradient strength 
% (in G/cm).  The resulting RF matrix will have the same number of columns 
% specifying phase, magnitude, timestep, and gradient strength (if 
% supplied).  This function allows you to specify the column numbers of
% each waveform component in the input data, to account for different text
% file conventions.  
% 
% INPUTS:
% filename   = filename of the .txt file to read in.
% colOrder   = Vector specifying the order of the RF amplitude, phase, 
%              time-step and gradient columns of the input RF file.  If a 
%              column does not  exist, it should be labelled as a zero. For
%              example, if the input file has columns in the order [phase,
%              amplitude, gradient], then colOrder should be [2,1,0,3].  If 
%              the input has columns in the order [Amplitude, phase, time],
%              then colOrder should be [1,2,3,0] or [1,2,3].
%
% OUTPUTS:
% rf         = Input rf pulse waveform saved as a matlab array with 3 or 4
%               columns (phase, magnitude, duration, gradient (if supplied)).
% info       = Empty. Not required.  

function [rf,info]=io_readRFtxt(filename,colOrder)

% For these basic text files, the data can be read simply using the load
% function:

RF=load(filename);

%Figure out the column order:
if nargin<2
    phCol=2;
    ampCol=1;
    if size(RF,2)>2
        timeCol=3;
        if size(RF,2)>3
            gradCol=4;
        end
    end
else
    if length(colOrder) < size(RF,2)
        error('ERROR:  length of colOrder must be at least as large as the number of columns of the RF input!!  ABORTING!');
    else
        ampCol=colOrder(1);
        phCol=colOrder(2);
            if length(colOrder)>2
                timeCol=colOrder(3);
                    if length(colOrder)>3
                        gradCol=colOrder(4);
                    end
            end
    end
end

%store the phase as the 1st column:
rf(:,1)=RF(:,phCol); 

%store the amplitude as the 2nd column:
rf(:,2)=RF(:,ampCol); 

%Store the time waveform as the 3rd column:
if exist('timeCol')
    if timeCol
        rf(:,3)=RF(:,timeCol);
    else
        rf(:,3)=ones(length(RF(:,1)),1);
    end
else
    rf(:,3)=ones(length(RF(:,1)),1);
end

%Store the gradient waveform as the 4th column (if it exists):
if exist('gradCol')
    if gradCol
        rf(:,4)=RF(:,gradCol);
    end
end


info=[];


   