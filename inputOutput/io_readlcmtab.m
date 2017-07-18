%io_readlcmtab.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% out = io_readlcmtab(filename) 
% 
% DESCRIPTION:
% Reads a LCModel .table output file and stores the metabolite 
% concentrations into a matlab structure array.
% 
% INPUTS:
% filename   = filename of the LCModel .table file.
%
% OUTPUTS:
% out        = A structure containing the LCmodel concentration estimates 
%               and CRLB values for each metabolite.

function out=io_readlcmtab(filename)

%try to incorporate the header information into a structure called 'info'
fid=fopen(filename);
dollar_index=[];

while isempty(dollar_index)
    line=fgets(fid);
    dollar_index=findstr(line,'$$');
end

line=fgets(fid);
line=fgets(fid);


% Now begin to read the data.  LCModel table files have a % sign marking each
% line of Data.  Search for the semicolon on each line and read only the 
%data that preceeds it.  


while length(line)>1;
    out.(genvarname(strtrim(line(24:end))))=str2num(line(1:9));
    out.(genvarname(['d_' strtrim(line(24:end))]))=str2num(line(11:13));
    line=fgets(fid);
end

FWHM_index=findstr(line,'FWHM');
while isempty(FWHM_index)
    line=fgets(fid);
    FWHM_index=findstr(line,'FWHM');
end

equals_indices=findstr(line,'=');
ppm_index=findstr(line,'ppm');
SN_index=findstr(line,'S/N');

out.FWHM=str2num(line(equals_indices(1)+1:ppm_index-1));
out.SNR=str2num(line(equals_indices(2)+1:end));

   