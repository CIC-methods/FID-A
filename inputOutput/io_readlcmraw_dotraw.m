% io_readlcmraw_dotraw.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out=io_readlcmraw_dotraw(filename);
% 
% DESCRIPTION:
% Reads LCModel .RAW model spectrum file into the FID-A data structure format in MATLAB.
% 
% INPUTS:
% filename   = filename of LCModel raw file.
%
% OUTPUTS:
% out        = Input metabolite basis spectrum in FID-A structure format.

function [out]=io_readlcmraw_dotraw(filename)

% Begin to read the data.  
fid=fopen(filename);
linenum=1;
line=fgets(fid);
        
sw_index=findstr(line,'Sweep Width');
while isempty(sw_index);
    line=fgets(fid);
    sw_index=findstr(line,'Sweep Width');
end
equals_index=findstr(line,'=');
Hz_index=findstr(line,'Hz');
spectralwidth=str2num(line(equals_index+2:Hz_index-1));

n_index=findstr(line,'Vector Size');
while isempty(n_index);
    line=fgets(fid);
    n_index=findstr(line,'Vector Size');
end
equals_index=findstr(line,'=');
pts_index=findstr(line,'points');
vectorsize=str2num(line(equals_index+2:pts_index-1));

Bo_index=findstr(line,'B0 Field');
while isempty(Bo_index)
    line=fgets(fid);
    Bo_index=findstr(line,'B0 Field');
end
equals_index=findstr(line,'=');
T_index=findstr(line,'T');
Bo=str2num(line(equals_index+2:T_index-1));
hzpppm=42.577*Bo;
               

hdrEnd_index=findstr(line,'$END');
while isempty(hdrEnd_index);
    line=fgets(fid);
    hdrEnd_index=findstr(line,'$END');
end

line=fgets(fid);

% If the line is empty skip it
while line~=-1
    %dataline=line(1:semicol_index-2);
    [A,count, errmsg, nextindex] = sscanf(line, '%f', inf);
    % If read failed, output the error     
    if ~isempty(errmsg);
       fclose(fid);
       error('READLCMRAW_BASIS failed with read error: %s', errmsg);
    end
    % Store the read values into rf array
    RF(linenum) = A(1)+i*A(2);
    linenum = linenum + 1;
    line=fgets(fid);
end
sz=[vectorsize 1];
out.fids=RF';
out.specs=fftshift(ifft(out.fids));
out.sz=[vectorsize 1 1 1];
out.spectralwidth=spectralwidth;
out.sz=sz;
out.Bo=Bo;
out.dwelltime=1/out.spectralwidth;

f=[(-spectralwidth/2)+(spectralwidth/(2*sz(1))):spectralwidth/(sz(1)):(spectralwidth/2)-(spectralwidth/(2*sz(1)))];
ppm=-f/(Bo*42.577);
ppm=ppm+4.65;
out.ppm=ppm;

t=[out.dwelltime:out.dwelltime:vectorsize*out.dwelltime];
out.t=t;

txfrq=hzpppm*1e6;
out.txfrq=txfrq;

out.date=date;
out.seq='';

out.dims.t=1;
out.dims.coils=0;
out.dims.averages=0;
out.dims.subSpecs=0;

out.averages=1;

out.flags.writtentostruct=1;
out.flags.gotparams=1;
out.flags.leftshifted=1;
out.flags.filtered=0;
out.flags.zeropadded=0;
out.flags.freqcorrected=0;
out.flags.phasecorrected=0;
out.flags.averaged=1;
out.flags.addedrcvrs=1;
out.flags.subtracted=1;
out.flags.writtentotext=1;
out.flags.downsampled=0;
out.flags.isISIS=0;

fclose(fid);
   
% RF=RF';
% rf(:,1)=RF(:,2)*180/pi;
% rf(:,2)=RF(:,1);
% rf(:,3)=ones(length(RF(:,1)),1);