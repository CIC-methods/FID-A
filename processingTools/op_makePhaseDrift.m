% op_makePhaseDrift.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% [out,phDrift]=op_makePhaseDrift(in,totalDrift,noise);
% 
% DESCRIPTION:
% Add phase drift to a dataset containing multiple averages.  This is generally 
% used to generate simulated datasets with phase drift. 
% 
% INPUTS:
% in         = input data in matlab structure format.
% totalDrift = total amount of phase drift (in degrees) to add over the whole scan.
% noise      = the standard deviation of noise to add to the phase drift
%             function.

function [out,phDrift]=op_makePhaseDrift(in,totalDrift,noise);
%out=op_makedrift(in,totalDrift);

%First make the matrices needed for multiplication
if totalDrift
    ph=[0:totalDrift/(in.sz(in.dims.averages)-1):totalDrift];
else
    ph=zeros(1,in.sz(in.dims.averages));
end
phnoise=[0 noise*randn(1,length(ph)-1)];
phDrift=ph+phnoise;
PH=repmat(phDrift,in.sz(1),1);

%Now apply the drift to the fids;
fids=in.fids.*exp(-1i*PH*pi/180);

%Now re-calculate specs using ifft
specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);

%FILLING IN DATA STRUCTURES
out=in;
out.fids=fids;
out.specs=specs;

%FILLING IN THE FLAGS
out.flags=in.flags;
