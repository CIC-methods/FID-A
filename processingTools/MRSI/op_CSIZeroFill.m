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
        x (1,1) double {mustBeNumeric,mustBeInteger, mustBeGreaterThan(x, -1)} = 0
        y (1,1) double {mustBeNumeric,mustBeInteger,mustHaveDims(in,x,y), mustBeGreaterThan(y,-1)} = 0
    end
    dim_size = in.sz;
    x_split_left = ceil(x/2);
    x_split_right = floor(x/2);
    y_split_top = ceil(y/2);
    y_split_bottom = floor(y/2);
    fids = in.fids;
    
    padding_dim = dim_size;
    padding_dim(in.dims.x) = x_split_left;
    padding = zeros(padding_dim);
    fids = cat(in.dims.x, padding, fids);
    
    padding_dim = size(fids);
    padding_dim(in.dims.x) = x_split_right;
    padding = zeros(padding_dim);
    fids = cat(in.dims.x, fids, padding);
    
    padding_dim = size(fids);
    padding_dim(in.dims.y) = y_split_top;
    padding = zeros(padding_dim);
    fids = cat(in.dims.y, padding, fids);
    
    padding_dim = size(fids);
    padding_dim(in.dims.y) = y_split_bottom;
    padding = zeros(padding_dim);
    fids = cat(in.dims.y, fids, padding);
    
    

    
    out = in;
    out.fids = fids;
    out.sz = size(fids);
    out.deltaX = out.fovX/size(fids,in.dims.x);
    out.deltaY = out.fovY/size(fids,in.dims.y);
    
    prev_affine_scale_matrix = eye(4);
    prev_affine_scale_matrix(1,1) = in.deltaX;
    prev_affine_scale_matrix(2,2) = in.deltaY;
    prev_affine_scale_matrix(3,3) = in.deltaZ;
    
    affine_scale_matrix = eye(4);
    affine_scale_matrix(1,1) = out.deltaX;
    affine_scale_matrix(2,2) = out.deltaY;
    affine_scale_matrix(3,3) = out.deltaZ;
    
    out.affine_matrix = out.affine_matrix*inv(prev_affine_scale_matrix)*affine_scale_matrix;
    
    
    
end

function mustHaveDims(in, x, y)
    if (in.dims.x == 0 && x > 0)
        eidType = 'mustExist:XDimensionDoesNotExist';
        msgType = 'Input of x must be zero or x dimension must exist in CSI strucutre';
        throwAsCaller(MException(eidType,msgType))
    elseif(in.dims.x == 0 && y > 0)
        eidType = 'mustExist:YDimensionDoesNotExist';
        msgType = 'Input of y must be zero or y dimension must exist in CSI strucutre';
        throwAsCaller(MException(eidType,msgType))
    end
end
