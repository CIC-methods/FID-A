%io_readlcmcoord.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% out = io_readlcmcoord(filename,metab) 
% 
% DESCRIPTION:
% Reads a LCModel .coord file and extracts the desired part.
% 
% INPUTS:
% filename   = filename of the LCModel .coord file.
% part       = Which metabolite fit to extract from the .coord file - 
%              The abbreviated metabolite name should be given (ie.
%              'Cr','PCr','Glu','GABA',etc.)
%
% OUTPUTS:
% out        = Desired metabolite component in simplified FID-A structure format.
% 
%Esin Ozturk Isik, Bogazici University, 01/06/2017.  Corrected the handling
%of non-existing metabolites.

function [out]=io_readlcmcoord(filename,metab)

fid=fopen(filename);
linenum=1;
line=fgets(fid);

ppmstr='points on ppm-axis =';
ppm_index=findstr(line,ppmstr);
while isempty(ppm_index)
    line=fgets(fid);
    ppm_index=findstr(line,ppmstr);
end

Npts=str2num(line(2:5));

line=fgets(fid);
while linenum<=(Npts/10)
    [A,count, errmsg, nextindex] = sscanf(line, '%f', inf);
    A;
    ppm((linenum-1)*10+1:linenum*10,1)=A;
    linenum=linenum+1;
    line=fgets(fid);
end


metab_str='';
metab_str(2:1+length(metab))=metab;
metab_str(11:17)='Conc. =';
metab_str
metab_index=findstr(line,metab);
while isempty(metab_index) & ~feof(fid)
    line=fgets(fid);
    
    metab_index=findstr(line,metab);
end
linenum=1;
line=fgets(fid);

if ~feof(fid)
while linenum<=(Npts/10) 
    [A,count, errmsg, nextindex] = sscanf(line, '%f', inf);
    A;
    specs((linenum-1)*10+1:linenum*10,1)=A;
    linenum=linenum+1;
    line=fgets(fid);
end

out.ppm=ppm;
out.specs=specs;
else
out.ppm=[];
out.specs=[];
end
