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


function MRSIStruct = op_CSIAverage(MRSIStruct)
    checkArguments(MRSIStruct);
    data = getData(MRSIStruct);

    averageDimension = getDimension(MRSIStruct, 'averages');
    data = squeeze(mean(data, averageDimension));

    MRSIStruct = setData(MRSIStruct, data);
    MRSIStruct = setFlags(MRSIStruct, 'averaged',  1);
    MRSIStruct = removeDimension(MRSIStruct, 'averages');
end

function checkArguments(in)
    if(in.dims.averages == 0)
        error('no dims to average');
    end
    if(in.flags.averaged)
        error('already averaged')
    end
end