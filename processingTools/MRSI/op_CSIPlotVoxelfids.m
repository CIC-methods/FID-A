% op_CSIPlotVoxelfids.m
% Brenden Kadota, SunnyBrook Hosptial 2021.
%
% USAGE:
% [out]=op_CSIPlotVoxelfids(in, x, y, ppmmin, ppmmax, xlab, ylab, title);
% 
% DESCRIPTION:
% Plots the fids of one voxel in the CSI FID-A structure. Additional
% arguments augments the output plot. Uses the op_plotfid function.
% 
% INPUTS:
% in   = FID-A CSI structure
% x    = x index
% y    = y index
% tmax = maximum time to plot
% xlab = x label for the plot
% ylab = y label for the plot
% title = title
%
% OUTPUTS:
% fig   = figure handle
function fig = op_CSIPlotVoxelfids(in, x, y, tmax, xlab, ylab, title)
    if getFlags(in,'spatialFT') == 0
        error('please fourier transfrom along the spatial and spectral dimensions')
    end
    if ~exist('x', 'var') || ~exist('y', 'var')
        error('please provide an x and y coordinate');
    end
    
    single_vox = op_CSItoMRS(in, x, y);
    if(~exist('tmax', 'var'))
        fig = op_plotfid(single_vox);
    elseif(~exist('xlab', 'var'))
        fig = op_plotfid(single_vox, tmax);
    elseif(~exist('ylab', 'var'))
        fig = op_plotfid(single_vox, tmax, xlab);
    elseif(~exist('title', 'var'))
        fig = op_plotfid(single_vox, tmax, xlab, ylab);
    else
        fig = op_plotfid(single_vox, tmax, xlab, ylab, title);
    end
    
end
