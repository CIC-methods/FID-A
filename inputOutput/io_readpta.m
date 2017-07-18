%io_readpta.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% [rf,info]=io_readpta(filename)
% 
% DESCRIPTION:
% Read a Siemens .pta file into matlab.  The resulting RF matrix will have 
% 2 columns specifying magnitude and phase.
% 
% INPUTS:
% filename   = filename of the .pta file to read in.  
%
% OUTPUTS:
% rf        = Input rf pulse waveform saved as a matlab array with 2
%               columns (magnitude and phase).
% info      = Not used.  

function [rf,info]=io_readpta(filename)

%try to incorporate the header information into a structure called 'info'
fid=fopen(filename);
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


% Now begin to read the data.  PTA files have a semicolon marking each
% line of Data.  Search for the semicolon on each line and read only the 
%data that preceeds it.  
RF=zeros(2,0);
linenum=1;
semicol_index=findstr(line,';');
if isempty(semicol_index)
    line=fgets(fid);
    semicol_index=findstr(line,';');
end

% If the line is empty skip it
while ~isempty(semicol_index)
    dataline=line(1:semicol_index-2);
    [A,count, errmsg, nextindex] = sscanf(dataline, '%f', inf);
    % If read failed, output the error     
    if ~isempty(errmsg)
       fclose(fid)
       error('READPTA failed with read error: %s', errmsg);
    end
    % Store the read values into rf array
    RF(:,linenum) = A;
    linenum = linenum + 1;
    line=fgets(fid);
    semicol_index=findstr(line,';');
end
fclose(fid);
   
RF=RF';
rf(:,1)=RF(:,2)*180/pi;
rf(:,2)=RF(:,1);
rf(:,3)=ones(length(RF(:,1)),1);
   