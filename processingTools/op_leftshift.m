% op_leftshift.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out=op_leftshift(in,ls);
% 
% DESCRIPTION:
% Remove leading datapoints from the fid to get rid of 1st order phase
% errors.
% 
% INPUTS:
% in     = input data in matlab strucure format.
% ls     = number of points to remove from the beginning of the fid.
%
% OUTPUTS:
% out    = Output following left shifting.  

function out=op_leftshift(in,ls);

if in.flags.leftshifted
    cont=input('WARNING:  Left shifting has already been performed!  Continue anyway?  (y or n)','s');
    if cont=='y'
        %continue;
    else
        error('STOPPING');
    end
end

fids=in.fids;

sz=in.sz;
sz(1)=sz(1)-ls;
fids=reshape(fids(ls+1:end,:),sz);

%re-calculate Specs using fft
specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);


%Calculate t and ppm arrays using the calculated parameters:
f=[(-in.spectralwidth/2)+(in.spectralwidth/(2*sz(1))):in.spectralwidth/(sz(1)):(in.spectralwidth/2)-(in.spectralwidth/(2*sz(1)))];
ppm=-f/(in.Bo*42.577);
ppm=ppm+4.65;

%t=[0:in.dwelltime:(sz(1)-1)*in.dwelltime];
t=[in.dwelltime:in.dwelltime:sz(1)*in.dwelltime];

    
%FILLING IN DATA STRUCTURE
out=in;
out.fids=fids;
out.specs=specs;
out.sz=sz;
out.ppm=ppm;  
out.t=t;    

%FILLING IN THE FLAGS
out.flags=in.flags;
out.flags.writtentostruct=1;
out.flags.leftshifted=1;
