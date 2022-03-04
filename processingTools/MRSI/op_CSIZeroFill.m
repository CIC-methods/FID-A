% op_CSIZeroFill.m
% Brenden Kadota, SunnyBrook Hosptial 2021.
%
% USAGE:
% [out]=io_loadspec_rda(rda_filename);
% 
% DESCRIPTION:
% Zero fills in the spatial domain
% 
% op_CSIZeroFill pads zeros to the spectral domain. If the padded zeros
% result in an odd number of points, it will always pad onto the left side.
% 
% INPUTS:
% in   = FID-A CSI structure
% x    = Number of points to pad in the x direction
% y    = Number of points to pad in the y direction
%
% OUTPUTS:
% out        = Output CSI structure with zero padding
function out = op_CSIZeroFill(in, x, y)
    arguments
        in (1,1) struct 
        x (1,1) double {mustBeInteger, mustBeGreaterThan(x, -1)} = 0
        y (1,1) double {mustBeInteger, mustBeGreaterThan(y, -1)} = 0
    end
    if(getFlags(in, 'spatialft') == true)
        error('Please keep the data in k space before zero filling')
    end
    dim_size = getSize(in);
    x_split_left = ceil(x/2);
    x_split_right = floor(x/2);
    y_split_top = ceil(y/2);
    y_split_bottom = floor(y/2);
    data = getData(in);
    
    xDimension = getDimension(in , 'kx');
    data = pad(dim_size, xDimension, x_split_left, x_split_right, data);
    
    yDimension = getDimension(in, 'ky');
    dim_size = size(data);
    data = pad(dim_size, yDimension, y_split_top, y_split_bottom, data);
    

    out = in;
    out = setData(out, data);
    fovX = getFov(out, 'x');
    fovY = getFov(out, 'y');
    xCoordinates = createCoordinates(fovX/2, fovX/getSizeFromDimensions(out, {'kX'}));
    yCoordinates = createCoordinates(fovY/2, fovY/getSizeFromDimensions(out, {'ky'}));
    out = setVoxelSize(out, 'x', xCoordinates(2) - xCoordinates(1));
    out = setVoxelSize(out, 'y', yCoordinates(2) - yCoordinates(1));
    out = setCoordinates(out, 'x', xCoordinates);
    out = setCoordinates(out, 'y', yCoordinates);
    
    prev_affine_scale_matrix = eye(4);
    prev_affine_scale_matrix(1, 1) = getVoxSize(in, 'x');
    prev_affine_scale_matrix(2, 2) = getVoxSize(in, 'y');
    prev_affine_scale_matrix(3, 3) = getVoxSize(in, 'z');
    prev_affine_scale_matrix(4, 4) = 1;
    
    affine_scale_matrix = eye(4);
    affine_scale_matrix(1, 1) = getVoxSize(in, 'x');
    affine_scale_matrix(2, 2) = getVoxSize(in, 'y');
    affine_scale_matrix(3, 3) = getVoxSize(in, 'z');
    affine_scale_matrix(4, 4) = 1;

    affineMatrix = getAffineMatrix(out);
    newAffineMatrix = affineMatrix*inv(prev_affine_scale_matrix)*affine_scale_matrix;
    out = setAffineMatrix(out, newAffineMatrix);
    out = setFlags(out, 'zeropadded', 1);
end

function padding = getPadding(dataSize, dimension, splitAmount)
    paddingSize = dataSize;
    paddingSize(dimension) = splitAmount;
    padding = zeros(paddingSize);
end

function data = pad(dim_size, dimension, leftPad, rightPad, data)
    paddingLeft = getPadding(dim_size, dimension, leftPad);
    paddingRight = getPadding(dim_size, dimension, rightPad);
    data = cat(dimension, paddingLeft, data, paddingRight);
end
