%op_CSIAverage.m
%Brenden Kadota, Sunnybrook, 2021
%
% USAGE:
% [out]=op_CSIAverage(in);
% 
% DESCRIPTION:
% Averages CSI FID-A structure along average dimension.
% 
% INPUTS:
% in   = Input CSI FID-A structure
%
% OUTPUTS:
% out = Averaged CSI FID-A structure


function out = op_CSIAverage(in)
    if(in.dims.averages == 0)
        error('no dims to average');
    end
    if(in.flags.averaged)
        error('already averaged')
    end
    out = in;
    out.specs = squeeze(mean(in.specs, in.dims.averages));
    out.flags.averaged = 1;
    
    fn = fieldnames(in.dims);
    %updating MRSI parameters
    for names = 1:length(fn)
        if(out.dims.(fn{names}) > out.dims.averages)
            out.dims.(fn{names}) = out.dims.(fn{names}) - 1;
        end
    end
    out.dims.averages = 0;
    out.sz = size(out.specs);
end