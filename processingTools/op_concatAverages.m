% op_concatAverages.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out=op_concatAverages(in1,in2);
% 
% DESCRIPTION:
% Concatenate two scans along the averages dimension.  Two scans with 50
% averages each will now look like a single scan with 100 averages.  For
% looping, it is possible to enter a blank structure as the first input.
% In this case, the function will simply return the second input.
% 
% INPUTS:
% in1    = first input in matlab structure format.
% in2    = second input in matlab structure format.
%
% OUTPUTS:
% out    = Output following concatenation of inputs along the 
%           averages dimension.  

function out=op_concatAverages(in1,in2);


if isempty(fieldnames(in1)) %Check if the first argument to be a null string
    out=in2; %If so, return in2;
else
    if in1.dims.averages~=0 && in2.dims.averages~=0 && in1.dims.averages ~= in2.dims.averages
        %RULE:  Averages dimension must be the same for both inputs UNLESS one of
        %the inputs has an averages dimension of zero:
        error('averages dimensions must be the same for both inputs');
    end
    
    if in1.dims.averages==0 && in2.dims.averages==0
        dim=2;
    elseif in1.dims.averages==0
        dim=in2.dims.averages;
    elseif in2.dims.averages==0;
        dim=in1.dims.averages;
    else
        dim=in1.dims.averages;
    end
    
    fids=cat(dim,in1.fids,in2.fids);
    specs=cat(dim,in1.specs,in2.specs);
    sz=size(fids);
    
    %FILLING IN DATA STRUCTURE
    out=in1;
    out.fids=fids;
    out.specs=specs;
    out.sz=sz;
    out.averages=in1.averages+in2.averages;
    out.rawAverages=in1.rawAverages+in2.rawAverages;
    if in1.dims.averages==0
        out.dims.averages=dim;
    end
    
    %FILLING IN THE FLAGS
    out.flags=in1.flags;
    out.flags.writtentostruct=1;
    out.flags.averaged=0;
end



