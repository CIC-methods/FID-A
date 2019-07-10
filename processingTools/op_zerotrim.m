% op_zerotrim.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out=op_zerotrim(in,numPointsToTrim);
% 
% DESCRIPTION:
% Remove zeros (or even non-zero data points) from the end of the fid.
% 
% INPUTS:
% in                = input data in matlab structure format.
% numPointsToTrim   = The number of points to trim from the end of the fid.
%
% OUTPUTS:
% out               = Output dataset following truncation in time domain.  

function out=op_zerotrim(in,numPointsToTrim);

if ~in.flags.zeropadded
    disp('WARNING:  you are trimming points from the end of the FID even though zero padding has not been performed!');
end

%calculate how many zeros to leave
zp=(in.sz(1)-numPointsToTrim);

sz=in.sz;
sz(1)=sz(1)-numPointsToTrim;
%Trim zeros
fids=reshape(in.fids(1:zp,:),sz);

%Calculate Specs using fft
specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);

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

%FILLING IN THE FLAGS
out.flags=in.flags;
out.flags.writtentostruct=1;
out.flags.zeropadded=0;