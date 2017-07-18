% op_phaseAlignAverages_fd.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% [out,phs]=op_phaseAlignAverages_fd(in,minppm,maxppm,Npts,avg,weighting)
% 
% DESCRIPTION:
% Perform time-domain spectral registration using a limited range of
% frequencies and using only phase adjustment (no frequency adjustment).
% This is rarely used.  
% 
% INPUTS:
% in         = Input data structure.
% minppm     = Minimum of frequency range (ppm).
% maxppm     = Maximum of frequnecy range (ppm).
% Npts       = Number of points in time domain to use for alignment.
% avg        = Align averages to the average of the averages ('y'), or the 
%              first average in the series ('n'); 
% weighting	 = (Optional) Apply less weight to the later points of the fid?
%
% OUTPUTS:
% out        = Output following alignment of averages.  
% phs        = Vector of phases (in degrees) used for alignment.

function [out,phs]=op_phaseAlignAverages_fd(in,minppm,maxppm,Npts,avg,weighting)

if ~in.flags.addedrcvrs
    error('ERROR:  I think it only makes sense to do this after you have combined the channels using op_addrcvrs.  ABORTING!!');
end

if nargin<6
    weighting='n';
end

if in.dims.subSpecs==0
    B=1;
else
    B=in.sz(in.dims.subSpecs);
end

phs=zeros(in.sz(in.dims.averages),B);
in_avg=op_averaging(in);
fids=zeros(size(in.fids));
if weighting=='y' || weighting=='Y'
    temp=op_freqrange(in_avg,minppm,maxppm);
    wgt_func=abs(temp.fids).^2;
else
    temp=op_freqrange(in_avg,minppm,maxppm);
    wgt_func=ones(size(temp.fids));
end

for m=1:B
    if avg=='y'||avg=='Y'
        base=in_avg;
        base=op_freqrange(base,minppm,maxppm);
        base=base.fids(1:Npts,m);
        begin=1;
    else
        base=op_freqrange(in,minppm,maxppm);
        base=base.fids(1:Npts,1,m);
        begin=2;
        fids(:,1,m)=in.fids(:,1,m);
    end
    for n=begin:in.sz(in.dims.averages)
        datarange=op_freqrange(in,minppm,maxppm);
        phsdiffs=(phase(base)-phase(datarange.fids(1:Npts,n,m)))*180/pi;
        phs(n,m)=mean(phsdiffs.*wgt_func(1:Npts,m))/mean(wgt_func(1:Npts,m));
        fids(:,n,m)=addphase(in.fids(:,n,m),phs(n,m));
    end
end


%re-calculate Specs using fft
specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);


%FILLING IN DATA STRUCTURE
out=in;
out.fids=fids;
out.specs=specs;

%FILLING IN THE FLAGS
out.flags=in.flags;
out.flags.writtentostruct=1;
out.flags.freqcorrected=1;

end
