%op_rmNworstaverages.m
%Jamie Near, McGill University 2014.
%
%USAGE:
%[out,metric,badAverages]=op_rmNworstaverages(in,n);
%
%DESCRIPTION:
%Removes motion corrupted averages from a dataset containing multiple
%averages.  The N most badly motion corrupted averages are discarded.
%
%INPUTS:
%in         = input data in matlab structure format
%n          = number of bad averages to remove

function [out,metric,badAverages]=op_rmNworstaverages(in,n);

if in.flags.averaged
    error('ERROR:  Averaging has already been performed!  Aborting!');
end

if ~in.flags.addedrcvrs
    error('ERROR:  Receivers should be combined first!  Aborting!');
end

%first, make a metric by subtracting all averages from the first average, 
%and then taking the sum of all all the spectral points.  
if in.dims.subSpecs>0
    SS=in.sz(in.dims.subSpecs);
else
    SS=1;
end
infilt=op_filter(in,10);
inavg=op_averaging(infilt);
for n=1:in.sz(in.dims.averages)
    for m=1:SS
            metric(n,m)=sum((real(infilt.specs(:,n,m))-(real(inavg.specs(:,m))/in.sz(in.dims.averages))).^2);
    end
end
for m=1:SS
    P(m,:)=polyfit([1:in.sz(in.dims.averages)]',metric(:,m),2);
    figure
    plot([1:in.sz(in.dims.averages)],metric(:,m),'.',[1:in.sz(in.dims.averages)],polyval(P(m,:),[1:in.sz(in.dims.averages)]));
end

%find the average and standard deviation of the metric
avg=mean(metric);
stdev=std(metric);

%Now make a mask that represents the locations of the averages 
%whose metric values are more than nsd standard deviations away from the 
%mean metric value.

for n=1:SS
    %mask(:,n)=metric(:,n)>(avg(n)+(nsd*stdev(n))) | metric(:,n)<(avg(n)-(nsd*stdev(n)));
    %mask(:,n)=metric(:,n)>(polyval(P(n,:),[1:in.sz(in.dims.averages)])'+(nsd*stdev(n))) | metric(:,n)<(polyval(P(n,:),[1:in.sz(in.dims.averages)])'-(nsd*stdev(n)));
    %mask(:,n)=metric(:,n)>(polyval(P(n,:),[1:in.sz(in.dims.averages)])'+(nsd*stdev(n)));
    mask(:,n)=(metric(:,n)-polyval(P(n,:),[1:in.sz(in.dims.averages)])')==max((metric(:,n)-polyval(P(n,:),[1:in.sz(in.dims.averages)])'));
end

%Unfortunately, if one average is corrupted, then all of the subspecs
%corresponding to that average have to be thrown away.  Therefore, take the
%minimum intensity projection along the subspecs dimension to find out
%which averages contain at least one corrupted subspec:
if size(mask,2)>1
    mask=sum(mask')'>0;
end

%now the corrupted and uncorrupted average numbers are given by:
badAverages=find(mask);
goodAverages=find(~mask);

%make a new fids array containing only good averages
fids=in.fids(:,goodAverages,:,:);

%%re-calculate Specs using fft
specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);

%re-calculate the sz variable
sz=size(fids);

%FILLING IN DATA STRUCTURE
out=in;
out.fids=fids;
out.specs=specs;
out.sz=sz;
out.averages=length(goodAverages) * in.rawSubspecs;

%FILLING IN THE FLAGS
out.flags=in.flags;
out.flags.writtentostruct=1;
