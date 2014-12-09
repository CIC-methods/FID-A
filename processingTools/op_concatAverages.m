%op_concatAverages.m
%Jamie Near, McGill University 2014.
%
%USAGE:
%out=op_concatAverages(in1,in2);
%
%DESCRIPTION:
%Concatenate two scans along the averages dimension.  Two scans with 50
%averages each will now look like a single scan with 100 averages.
%
%INPUTS:
%in1    = first input in matlab structure format.
%in2    = second input in matlab structure format.

function out=op_concatAverages(in1,in2);

if in1.dims.averages ~= in2.dims.averages
    error('averages dimensions must be the same for both inputs');
end

fids=cat(in1.dims.averages,in1.fids,in2.fids);
specs=cat(in1.dims.averages,in1.fids,in2.fids);
sz=size(fids);

%FILLING IN DATA STRUCTURE
out=in1;
out.fids=fids;
out.specs=specs;
out.sz=sz;   
out.averages=in1.averages+in2.averages;
out.rawAverages=in1.rawAverages+in2.rawAverages;

%FILLING IN THE FLAGS
out.flags=in1.flags;
out.flags.writtentostruct=1;



