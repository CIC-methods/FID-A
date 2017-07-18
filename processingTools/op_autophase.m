% op_autophase.m
% Jamie Near, McGill University 2015.
% 
% USAGE:
% [out,phaseShift]=op_autophase(in,ppmmin,ppmmax,ph,dimNum);
% 
% DESCRIPTION:
% Search for the peak located between ppmmin and ppmmax, and then phase the
% spectrum so that that peak reaches the desired phase.
% 
% INPUTS:
% in         = input data in matlab structure format.
% ppmmin     = minimum of ppm search range.
% ppmmax     = maximum of ppm search range.
% ph         = desired phase value in degrees [optional.  Default=0].
% dimNum     = which subSpec dimension to use for phasing? [Only for use in
%              data with multiple subSpectra].  
%
% OUTPUTS:
% out        = Output following automatic phasing.
% phaseShift = The phase shift (in degrees) that was applied.

function [out,phShft]=op_autophase(in,ppmmin,ppmmax,ph,dimNum);


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
        dimNum=input('Input which subspectrum to use for phasing: ');
        if nargin<4
            ph=0;
        end
    end
else
    if nargin<4
        ph=0;
    end
    dimNum=1;
end

%Zeropad the data if it hasn't already been done
if ~in.flags.zeropadded
    in_zp=op_zeropad(in,10);
else
    in_zp=in;
end

%Narrow the frequency range:
in_zp=op_freqrange(in_zp,ppmmin,ppmmax);

%Find the ppm of the maximum peak magnitude within the given range:
ppmindex=find(abs(in_zp.specs(:,dimNum))==max(abs(in_zp.specs(:,dimNum))));

%now do automatic zero-order phase correction (Use Creatine Peak):
ph0=-phase(in_zp.specs(ppmindex,dimNum))*180/pi;

%Now phase shift the dataset so that the desired peak has the correct phase:
phShft = ph0 + ph;
out=op_addphase(in,phShft);
