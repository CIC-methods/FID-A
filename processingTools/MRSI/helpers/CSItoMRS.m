% CSItoTwix.m
% converts a CSI twix to a single voxel FID-A object wich allows for the single 
% voxel FID-A porcessing tools to be used
% 
% USAGE:
% out = CSItoMRS(in, xCoordinate, yCoordinate)
%
% INPUT: 
% in = FID-A CSI object
% xCoordinate = X coordinate wanted to be converted
% yCoordinate = Y coorindate wanted to be converted
%
% OUTPUT:
% out = single voxel FID-A object

function out = CSItoMRS(in, xCoordinate, yCoordinate, averageNum)
    if(in.flags.spatialFT == 0)
        error('please fourier transform along the spatial and spectral dimension before using this function');
    end
    out = in;
    idx = cell(length(in.sz),1);
    for i = 1:length(idx)
        idx{i} = ':';
    end
    idx{in.dims.x} = xCoordinate;
    idx{in.dims.y} = yCoordinate;
    
    if(in.dims.averages && exist('averageNum', 'var'))
        idx{in.dims.averages} = averageNum;
    end
    out.specs = in.specs(idx{:});
    out.specs = squeeze(out.specs);
    if(in.flags.spectralFT == 1)
    %inverse fourrier transform to get fids at x and y coordiantes
        out.fids = squeeze(ifft(ifftshift(out.specs, in.dims.t), [], in.dims.t));
    else
        out.fids = out.specs;
        out.specs = fftshift(fft(out.specs, [], in.dims.t), in.dims.t);
    end
    %adjust structure to match FID-A object
    out.sz = size(out.fids);
    fn = fieldnames(in.dims);
    for name = 1:length(fn)
        if(out.dims.(fn{name}) > in.dims.x)
            out.dims.(fn{name}) = out.dims.(fn{name}) - 1;
        end
        if(out.dims.(fn{name}) > in.dims.y)
            out.dims.(fn{name}) = out.dims.(fn{name}) - 1;
        end
    end
    out.dims.x = 0;
    out.dims.y = 0;
    
end
