%io_readRFBruk.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% rf=io_readRFBruk(filename);
% 
% DESCRIPTION:
% Read a Bruker RF pulse file into matlab.  The resulting RF matrix will have 
% 2 columns specifying magnitude and phase.
% 
% INPUTS:
% filename   = filename of the .pta file to read in. 
%
% OUTPUTS:
% rf         = Input rf pulse waveform saved as a matlab array with 2
%               columns (magnitude and phase).

function rf=io_readRFBruk(filename)

%try to incorporate the header information

%FIND THE PULSE TYPE ('Inversion','Refocusing' or 'Excitation')
fid=fopen(filename);
line=fgets(fid);
exmode_index=findstr(line,'##$SHAPE_EXMODE');
while isempty(exmode_index)
    line=fgets(fid);
    exmode_index=findstr(line,'##$SHAPE_EXMODE');
end
equals_index=findstr(line,'=');
pulseType=line(equals_index+2:end)
fclose(fid);

%Find the flip angle in degrees
fid=fopen(filename);
line=fgets(fid);
rotation_index=findstr(line,'##$SHAPE_TOTROT');
while isempty(rotation_index)
    line=fgets(fid);
    rotation_index=findstr(line,'##$SHAPE_TOTROT');
end
equals_index=findstr(line,'=');
flipangle=str2num(line(equals_index+2:end))
fclose(fid);

%Find the R value
fid=fopen(filename);
line=fgets(fid);
R_index=findstr(line,'##$SHAPE_BWFAC');
while isempty(R_index)
    line=fgets(fid);
    R_index=findstr(line,'##$SHAPE_BWFAC');
end
equals_index=findstr(line,'=');
R=str2num(line(equals_index+2:end))
fclose(fid);

%Find the integral factor:
fid=fopen(filename);
line=fgets(fid);
integ_index=findstr(line,'##$SHAPE_INTEGFAC');
while isempty(integ_index)
    line=fgets(fid);
    integ_index=findstr(line,'##$SHAPE_INTEGFAC');
end
equals_index=findstr(line,'=');
integ=str2num(line(equals_index+2:end))
fclose(fid);

%Find the number of points:
fid=fopen(filename);
line=fgets(fid);
npts_index=findstr(line,'##NPOINTS');
while isempty(npts_index)
    line=fgets(fid);
    npts_index=findstr(line,'##NPOINTS');
end
equals_index=findstr(line,'=');
npts=str2num(line(equals_index+2:end))
fclose(fid);

%Find the last line of the header:
fid=fopen(filename);
line=fgets(fid);
xypoints_index=findstr(line,'##XYPOINTS');
while isempty(xypoints_index)
    line=fgets(fid);
    xypoints_index=findstr(line,'##XYPOINTS');
end
line=fgets(fid);


% Now begin to read the data.  Bruker files are comma delimited.  
%Search for the semicolon on each line and read only the 
%data that preceeds it.  
RF=zeros(2,0);

% Loop through npts lines
for linenum=1:npts
    comma_index=findstr(line,',');
    dataline1=line(1:comma_index);
    dataline2=line(comma_index+2:end);
    % Store the read values into rf array
    RF(1,linenum) = str2num(dataline1);
    RF(2,linenum) = str2num(dataline2);
    linenum = linenum + 1;
    line=fgets(fid);
end
fclose(fid);
   
RF=RF';
rf(:,1)=RF(:,2);
rf(:,2)=RF(:,1);
rf(:,3)=ones(length(RF(:,1)),1);
   