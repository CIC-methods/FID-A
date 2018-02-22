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
function out = op_CSIRemoveWater(in)
if(in.flags.spectralFT == 0)
    error('FID-A Error: Input type invalid. Please fourier transform along the spectral dimension');
end
residual_error = 0;
    for x = 1:in.sz(in.dims.x)
        for y = 1:in.sz(in.dims.y)
            voxel = CSItoMRS(in, x, y);
            voxel_supressed = op_removeWater(voxel, [4.4, 5], 20, 1500, 0);
            in.specs(:, x, y) = voxel_supressed.specs;
            residual_error = residual_error + voxel_supressed.watersupp.residual_error;
        end
    end
    out = in;
    out.watersupp.residual_error = residual_error;

end