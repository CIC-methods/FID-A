%op_makeFreqDrift.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% [out,fDrift]=op_makeFreqDrift(in,totalDrift,noise);
% 
% DESCRIPTION:
% Add frequency drift to a dataset containing multiple averages.  This is 
% generally used to generate simulated datasets with phase drift. 
% 
% INPUTS:
% in         = input data in matlab structure format.
% totalDrift = total amount of frequency drift (in Hz) to add over the whole scan.
% noise      = the standard deviation of noise to add to the frequency drift
%             function.

function [out,fDrift]=op_makeFreqDrift(in,totalDrift,noise);
%

%First make the matrices needed for multiplication
T=repmat(in.t',[1 in.sz(2:end)]);
if totalDrift
    f=[0:totalDrift/(in.sz(in.dims.averages)-1):totalDrift];
else
    f=zeros(1,in.sz(in.dims.averages));
end
fnoise=[0 noise*randn(1,length(f)-1)];
fDrift=f+fnoise;
F=repmat(fDrift,in.sz(1),1);

%Now apply the drift to the fids;
fids=in.fids.*exp(1i*T.*F*2*pi);

%Now re-calculate specs using ifft
specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);

%FILLING IN DATA STRUCTURES
out=in;
out.fids=fids;
out.specs=specs;

%FILLING IN THE FLAGS
out.flags=in.flags;




