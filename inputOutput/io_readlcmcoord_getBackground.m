%io_readlcmcoord_getBackground.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% out = io_readlcmcoord_getBackground(filename,part) 
% 
% DESCRIPTION:
% Reads a LCModel .coord file and extracts the desired part.
% 
% INPUTS:
% filename   = filename of the LCModel .coord file.
% part       = Which part of the .coord file to extract - 'bg' extracts the
%             LCModel baseline signal, 'sp' extracts the spectrum, and 
%             'fit' extracts the fit.
%
% OUTPUTS:
% out        = Desired background component in simplified FID-A structure format.

function [out]=io_readlcmcoord_getBackground(filename,part)

fid=fopen(filename);
linenum=1;
line=fgets(fid)


%READ PPM AXIS
ppmstr='points on ppm-axis ='
ppm_index=findstr(line,ppmstr);
while isempty(ppm_index)
    line=fgets(fid);
    ppm_index=findstr(line,ppmstr);
end

Npts=str2num(line(2:5))

line=fgets(fid);
while linenum<=(Npts/10)
    [A,count, errmsg, nextindex] = sscanf(line, '%f', inf);
    A
    ppm((linenum-1)*10+1:linenum*10,1)=A;
    linenum=linenum+1;
    line=fgets(fid);
end



%GET BACKGROUND

linenum=1;
switch part
    case 'bg'
        bgstr='NY background values follow'
    case 'sp'
        bgstr='NY phased data points follow'
    case 'fit'
        bgstr='NY points of the fit to the data follow'
    otherwise
        error('ERROR:  part not found');
end
bg_index=findstr(line,bgstr);
while isempty(bg_index)
    line=fgets(fid);
    bg_index=findstr(line,bgstr);
end

line=fgets(fid);
while linenum<=(Npts/10)
    [A,count, errmsg, nextindex] = sscanf(line, '%f', inf);
    A;
    bg((linenum-1)*10+1:linenum*10,1)=A;
    linenum=linenum+1;
    line=fgets(fid);
end


out.ppm=ppm;
out.specs=bg;


