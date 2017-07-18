% op_unfilter.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out=op_unfilter(in,lb);
% 
% DESCRIPTION:
% Multiply the fid by an inverted exponential decay function to undo the effects of
% filtering.
% 
% INPUTS:
% in     = input data in matlab structure format.
% lb     = line narrowing factor (Hz). 
%
% OUTPUTS:
% out    = Output dataset following application of inverted exponential filter.

function out=op_unfilter(in,lb);

if in.flags.filtered
    cont=input('WARNING:  Line Broadening has already been performed!  Continue anyway?  (y or n)','s');
    if cont=='y'
        %continue;
    else
        error('STOPPING');
    end
end

fids=in.fids;

lbt=1/lb;

%Create a lorentzian filter
lor=(1/pi)*(((lbt/2)./((in.t.^2)+((lbt/2)^2))));
lor=lor/max(lor);
lor=1./lor;
%plot(in.t,lor);

%first make a bunch of vectors of ones that are the same lengths as each of
%the dimensions of the data.  Store them in a cell array for ease of use.
for n=1:length(in.sz)
    p{n}=ones(in.sz(n),1);
end

%Now, now take the lorentzian filter vector that we made earlier (lor) and use it
%to populate a filter array that has the same dimensions as the data.  To 
%do this, we have to use the ndgrid function, which is essentially the same
%as the meshgrid function, except in multiple dimensions.  b, c, and d are 
%dummy variables and are not used.  
if length(in.sz)==1
    fil=lor;
end
if length(in.sz)==2
    [fil,b]=ndgrid(lor,p{2});
end
if length(in.sz)==3
    [fil,b,c]=ndgrid(lor,p{2},p{3});
end
if length(in.sz)==4
    [fil,b,c,d]=ndgrid(lor,p{2},p{3},p{4});
end
if length(in.sz)==5
    [fil,b,c,d,e]=ndgrid(lor,p{2},p{3},p{4},p{5});
end

%Now multiply the data by the filter array.  
fids=fids.*fil;

%re-calculate Specs using fft
specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);

%FILLING IN DATA STRUCTURE
out=in;
out.fids=fids;
out.specs=specs;

%FILLING IN THE FLAGS
out.flags=in.flags;
out.flags.writtentostruct=1;
out.flags.filtered=1;
