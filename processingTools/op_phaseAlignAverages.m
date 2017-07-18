% op_phaseAlignAverages.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% [out,phs]=op_phaseAlignAverages(in,Npts,avg,weighting)
% 
% DESCRIPTION:
% Perform time-domain spectral registration using only phase adjustment 
% (no frequency adjustment).  This is rarely used.  
% 
% INPUTS:
% in         = Input data structure.
% Npts       = Number of points in time domain to use for alignment.
% avg        = Align averages to the average of the averages ('y'), or the 
%              first average in the series ('n'); 
% weighting	 = (Optional) Apply less weight to the later points of the fid?
%
% OUTPUTS:
% out        = Output following alignment of averages.  
% phs        = Vector of phases (in Degrees) used for alignment.

function [out,phs]=op_phaseAlignAverages(in,Npts,avg,weighting)

if ~in.flags.addedrcvrs
    error('ERROR:  I think it only makes sense to do this after you have combined the channels using op_addrcvrs.  ABORTING!!');
end

if nargin<4
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
    wgt_func=abs(in_avg.fids).^2;
    disp('weighting is on');
%     plot([1:length(in_avg.fids)],wgt_func,'.');
%     pause;
else
    wgt_func=ones(size(in_avg.fids));
    disp('weighting is off');
%     plot([1:length(in_avg.fids)],wgt_func,'.');
%     pause;
end

for m=1:B
    if avg=='y' || avg=='Y'
        disp('aligning all averages to the Average of the averages');
        base=in_avg;
        base=base.fids(1:Npts ,m);
        begin=1;
    else
        disp('aligning all averages to the first average');
        base=in.fids( 1:Npts ,1,m);
        begin=1;
        fids(:,1,m)=in.fids(:,1,m);
    end
    for n=begin:in.sz(in.dims.averages)
        phsdiffs=(phase(base)-phase(in.fids(1:Npts,n,m)))*180/pi;
        phs(n,m)=mean(phsdiffs.*wgt_func(1:Npts,m))/mean(wgt_func(1:Npts,m));
        if phs(n,m)>=0
           phs(n,m)=phs(n,m)-(360*fix((phs(n,m)+180)/360));
        else
           phs(n,m)=phs(n,m)-(360*fix((phs(n,m)-180)/360));
        end
        if n==1
           firstPhase(m)=phs(1,m);
        end
        phs(n,m)=phs(n,m)-firstPhase;
        fids(:,n,m)=addphase(in.fids(:,n,m),phs(n,m));
        %pause;
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