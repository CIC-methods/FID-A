% op_concatSubspecs.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out=op_concatSubspecs(in1,in2);
% 
% DESCRIPTION:
% Concatenate two scans along the subspecs dimension.  Two scans with 50
% averages each will now look like a single scan with 100 averages.
% 
% INPUTS:
% in1    = first input in matlab structure format.
% in2    = second input in matlab structure format.
%
% OUTPUTS:
% out    = Output following concatenation along the subspecs dimension.  


function out=op_concatSubspecs(in1,in2);

if in1.dims.subSpecs ~= in2.dims.subSpecs || in1.dims.t ~= in2.dims.t || in1.dims.coils ~= in2.dims.coils || in1.dims.averages ~=in2.dims.averages
    error('subSpecs dimensions must be the same for both inputs');
end

%if subspecs dimension is zero, make a new dimension for them
if in1.dims.subSpecs==0
    newSubspecsDim=max([in1.dims.t in1.dims.coils in1.dims.averages in1.dims.subSpecs in1.dims.extras])+1;
else
    newSubspecsDim=in1.dims.subSpecs;
end

fids=cat(newSubspecsDim,in1.fids,in2.fids);
specs=cat(newSubspecsDim,in1.specs,in2.specs);
sz=size(fids);

%FILLING IN DATA STRUCTURE
out=in1;
out.fids=fids;
out.specs=specs;
out.sz=sz;   
out.dims=in1.dims;
out.dims.subSpecs=newSubspecsDim;
    
out.rawSubspecs=in1.rawSubspecs+in2.rawSubspecs;
out.subspecs=in1.subspecs+in2.subspecs;

%FILLING IN THE FLAGS
out.flags=in1.flags;
out.flags.writtentostruct=1;
out.flags.subtracted=0;



