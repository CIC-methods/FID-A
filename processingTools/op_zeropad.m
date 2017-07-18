% op_zeropad.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out=op_zeropad(in,zpFactor);
% 
% DESCRIPTION:
% Apply zeropadding (a.k.a. zero-filling) to MRS data.
% 
% INPUTS:
% in         = input data in matlab structure format.
% zpFactor   = the factor by which the number of points in the fid will be
%             increased.  ie.  if zpFactor =2, then the number of zeros 
%             added to the end of the fid will be equal to the number of 
%             points in the original spectrum.
%
% OUTPUTS:
% out        = Output dataset following zeropadding.

function out=op_zeropad(in,zpFactor);

if in.flags.zeropadded
    cont=input('WARNING:  Zero padding has already been performed!  Continue anyway?  (y or n)','s');
    if cont=='y'
        %continue;
    else
        error('STOPPING');
    end
end


%calculate how many zeros to add
zp=ceil((in.sz(1)*zpFactor)-in.sz(1));

%Add zeros using MATLAB array zeropadding function;
fids=padarray(in.fids,zp,'post');

%Calculate Specs using fft
specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);

%recalculate the sz vector
sz=size(fids);


%Now re-calculate t and ppm arrays using the calculated parameters:
f=[(-in.spectralwidth/2)+(in.spectralwidth/(2*sz(1))):...
    in.spectralwidth/(sz(1)):...
    (in.spectralwidth/2)-(in.spectralwidth/(2*sz(1)))];

ppm=-f/(in.Bo*42.577);
ppm=ppm+4.65;

t=[0:in.dwelltime:(sz(1)-1)*in.dwelltime];


%FILLING IN DATA STRUCTURE
out=in;
out.fids=fids;
out.specs=specs;
out.sz=sz;
out.ppm=ppm;  
out.t=t;   
out.n=sz(1);

%FILLING IN THE FLAGS
out.flags=in.flags;
out.flags.writtentostruct=1;
out.flags.zeropadded=1;
