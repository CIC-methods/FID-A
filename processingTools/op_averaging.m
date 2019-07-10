% op_averaging.m
% Jamie Near, McGill University 2014.
% 
% USAGE:
% out=op_averaging(in);
% 
% DESCRIPTION:
% Combine the averages in a scan by adding the averages together and then 
% dividing by the number of averages.
% 
% INPUTS:
% in	= input data in matlab structure format.
%
% OUTPUTS:
% out   = Output following averaging.  

function out=op_averaging(in);

if in.flags.averaged || in.averages<2
    %DO NOTHING
    disp('WARNING: No averages found. Returning input without modification!');
    return;
end

if in.dims.averages==0
    %DO NOTHING
    disp('WARNING: No averages found. Returning input without modification!');
    out=in;
    return;
else
    
    %add the spectrum along the averages dimension;
    fids=sum(in.fids,in.dims.averages);
    fids=squeeze(fids);
    fids=fids/in.sz(in.dims.averages); %divide by number of averages;
    
    %re-calculate Specs using fft
    specs=fftshift(ifft(fids,[],in.dims.t),in.dims.t);
    
    %change the dims variables.
    if in.dims.t>in.dims.averages
        dims.t=in.dims.t-1;
    else
        dims.t=in.dims.t;
    end
    if in.dims.coils>in.dims.averages
        dims.coils=in.dims.coils-1;
    else
        dims.coils=in.dims.coils;
    end
    dims.averages=0;
    if in.dims.subSpecs>in.dims.averages
        dims.subSpecs=in.dims.subSpecs-1;
    else
        dims.subSpecs=in.dims.subSpecs;
    end
    if in.dims.extras>in.dims.averages
        dims.extras=in.dims.extras-1;
    else
        dims.extras=in.dims.extras;
    end
    
    
    %re-calculate the sz variable
    sz=size(fids);
    
    
    %FILLING IN DATA STRUCTURE
    out=in;
    out.fids=fids;
    out.specs=specs;
    out.sz=sz;
    out.dims=dims;
    out.averages=1;
    
    %FILLING IN THE FLAGS
    out.flags=in.flags;
    out.flags.writtentostruct=1;
    out.flags.averaged=1;

end



