%op_CSIRemoveWater.m
%Brenden Kadota, Sunnybrook 2021.
%
% USAGE:
% in=op_CSIRemoveWater(in);
%
% DESCRIPTION:
% This function uses op_removeWater to remove water signal from spectrum.
% op_removeWater ruses the HSVD method to remove the water which is
% described by H. BARKHUIJSEN et al. 1987.
%
% INPUTS:
% in        = CSI FID-A data structure
%
% OUTPUTS:
% out       = output
function MRSIStruct = op_CSIRemoveWater(MRSIStruct, wlim, kInitial, M)
    arguments
        MRSIStruct (1, 1) struct
        wlim (1, 2) double = [4.4 5]
        kInitial (1, 1) double = 30
        M (1, 1) double = floor(MRSIStruct.sz(MRSIStruct.dims.t)*.75);
    end
    residual_error = 0;
    [MRSIStruct, prevPermute, prevShape]= reshapeDimensions(MRSIStruct, {'t', 'y', 'x'});
    data = zeros(getSizeFromDimensions(MRSIStruct, {'t', 'y', 'x', 'extras'}));
    for e = 1:getSizeFromDimensions(MRSIStruct, {'extras'})
        for x = 1:getSizeFromDimensions(MRSIStruct, {'x'})
            for y = 1:getSizeFromDimensions(MRSIStruct, {'y'})
                voxel = op_CSItoMRS(MRSIStruct, x, y, 'Extra', e);
                voxel_supressed = op_removeWater(voxel, wlim, kInitial, M);
                data(:, y, x, e) = voxel_supressed.specs;
                residual_error = residual_error + voxel_supressed.watersupp.residual_error;
            end
        end
    end
    MRSIStruct = setData(MRSIStruct, data);
    MRSIStruct = reshapeBack(MRSIStruct, prevPermute, prevShape);


    MRSIStruct.watersupp.residual_error = residual_error;
end
