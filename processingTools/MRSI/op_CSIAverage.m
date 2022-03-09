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
    flag = checkArguments(MRSIStruct);
    %if check arguments returns true, abort
    if(flag); return; end

    data = getData(MRSIStruct);

    averageDimension = getDimension(MRSIStruct, 'averages');
    data = squeeze(mean(data, averageDimension));

    MRSIStruct = setData(MRSIStruct, data);
    MRSIStruct = setFlags(MRSIStruct, 'averaged',  1);
    MRSIStruct = removeDimension(MRSIStruct, 'averages');
end

% check if dims have been averaged already or dim dimension exists.
function flag = checkArguments(in)
    flag = false;
    if(in.dims.averages == 0)
        disp('No average Dim! Aborting')
        flag = true;
    end
    if(in.flags.averaged)
        disp('Already averaged! Aborting')
        flag = true;
    end
end