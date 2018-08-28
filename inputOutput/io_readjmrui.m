% io_readjmrui.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out=io_readjmrui(filename);
% 
% DESCRIPTION:
% Reads jMRUI .txt format into a Nx2 MATLAB array where the 2 rows are the
% time-domain and frequency domain data, respectively. 
% 
% INPUTS:
% filename   = filename of jMRUI .txt file.
%
% OUTPUTS:
% out        = Input dataset in FID-A structure format.

function [rf,info]=io_readjmrui(filename)

%try to incorporate the header information into a structure called 'info'
fid=fopen(filename);
line=fgets(fid);
line=fgets(fid);
colon_index=findstr(line,':');
if isempty(colon_index)
    line=fgets(fid);
    colon_index=findstr(line,':');
end

while ~isempty(colon_index)
    fieldname=line(1:colon_index-1);
    info.(genvarname(fieldname))=line(colon_index+1:end);
    line=fgets(fid);
    colon_index=findstr(line,':');
end

line=fgets(fid);
colon_index=findstr(line,':');

while ~isempty(colon_index)
    fieldname=line(1:colon_index-1);
    info.(genvarname(fieldname))=line(colon_index+1:end);
    line=fgets(fid);
    colon_index=findstr(line,':');
end
line=fgets(fid);
line=fgets(fid);

% Now begin to read the data.  PTA files have a semicolon marking each
% line of Data.  Search for the semicolon on each line and read only the 
%data that preceeds it.  
RF=zeros(1,4);
start_index=findstr(line,'Signal and FFT');
if isempty(start_index)
    line=fgets(fid);
    start_index=findstr(line,'Signal and FFT');
end
line=fgets(fid);
line=fgets(fid);
line=fgets(fid);

% If the line is empty skip it
for n=1:str2num(info.PointsInDataset)
    dataline=line(1:end);
    [A,count, errmsg, nextindex] = sscanf(dataline, '%f', inf);
    % If read failed, output the error     
    if ~isempty(errmsg)
       fclose(fid)
       error('READJMRUI failed with read error: %s', errmsg);
    end
    % Store the read values into rf array
    RF(n,:) = A;
    line=fgets(fid);
end
fclose(fid);
   




 rf(:,1)=RF(:,1)+i*RF(:,2);
 rf(:,2)=fftshift(ifft(rf(:,1)));
 
