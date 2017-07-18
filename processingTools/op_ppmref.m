% op_ppmref.m
% Jamie Near, McGill University 2015.
% 
% USAGE:
% [out,frqshift]=op_ppmref(in,ppmmin,ppmmax,ppmrefval,dimNum);
% 
% DESCRIPTION:
% Search for the peak located between ppmmin and ppmmax, and then give that
% peak a new ppm reference value.
% 
% INPUTS:
% in        = input data in matlab structure format.
% ppmmin    = minimum of ppm search range.
% ppmmax    = maximum of ppm search range.
% ppmrefval = new reference ppm value.
% dimNum    = which subspectrum to used for referencing (optional).
%
% OUTPUTS:
% out       = Output dataset following frequency shift.
% frqshift  = Frequency shift applied (in Hz).

function [out,frqshift]=op_ppmref(in,ppmmin,ppmmax,ppmrefval,dimNum);


if in.dims.coils>0
    error('ERROR:  Can not operate on data with multilple coils!  ABORTING!!')
end
if in.dims.averages>0
    error('ERROR:  Can not operate on data with multiple averages!  ABORTING!!');
end
if in.dims.extras>0
    error('ERROR:  Can not operate on data with extras dimension!  ABORTING!!');
end
if in.dims.subSpecs>0
    if nargin<5
        plot(in.ppm,in.specs);
        legend('subspec 1','subspec 2');
        dimNum=input('Input which subspectrum to use for referencing: ');
    end
else
    dimNum=1;
end

%Zeropad the data if it hasn't already been done
if ~in.flags.zeropadded
    in_zp=op_zeropad(in,10);
else
    in_zp=in;
end

%Find the ppm of the maximum peak magnitude within the given range:
ppmindex=find(abs(in_zp.specs(in_zp.ppm>ppmmin & in_zp.ppm<ppmmax,dimNum))==max(abs(in_zp.specs(in_zp.ppm>ppmmin & in_zp.ppm<ppmmax,dimNum))));
ppmrange=in_zp.ppm(in_zp.ppm>ppmmin & in_zp.ppm<ppmmax);
ppmmax=ppmrange(ppmindex);

%Now frequency shift the dataset so that the max peak appears at ppmrefval:
frqshift=(ppmmax-ppmrefval)*in.txfrq/1e6;
out=op_freqshift(in,(ppmmax-ppmrefval)*in.txfrq/1e6);
