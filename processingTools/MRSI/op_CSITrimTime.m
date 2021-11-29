% op_CSIZeroFill.m
% Brenden Kadota, SunnyBrook Hosptial 2021.
%
% USAGE:
% [out]=io_loadspec_rda(rda_filename);
% 
% DESCRIPTION:
% Removes the first the first time points in all fids until new_start_idx is reached
% 
% INPUTS:
% in   = FID-A CSI structure
% new_start_idx    = Number of points to pad in the x direction
%
% OUTPUTS:
% out        = Output CSI structure with zero padding

function out = op_CSITrimTime(in, new_start_idx)
    if(in.flags.spectralFT == 1)
        fids = ifft(fftshift(in.specs, in.dims.t), [], in.dims.t);
    else
        fids = in.specs;
    end
    permute_dims = nonzeros([in.dims.t, in.dims.x, in.dims.y, in.dims.coils, in.dims.averages]);
    permute_back(permute_dims) = 1:length(permute_dims);
    fids = permute(fids, permute_dims);
    sz = size(fids);
    fids = reshape(fids, size(fids, 1), []);
    fids = fids(new_start_idx:end, :);
    sz(1) = length(new_start_idx:sz(1));
    fids = reshape(fids, sz);
    fids = permute(fids, permute_back);
    
    if(in.flags.spectralFT == 1)
        fids = fftshift(fft(fids, [], in.dims.t), in.dims.t);
    end
    in.specs = fids;
    out = in;
    out.t = out.t(new_start_idx:sz(in.dims.t));
    out.sz = size(fids);
    
    lb = (-out.spectralwidth/2)+(out.spectralwidth/(2*out.sz(out.dims.t)));
    ub = (out.spectralwidth/2)-(out.spectralwidth/(2*out.sz(out.dims.t)));
    step = out.spectralwidth/(size(out.specs,1));
    
    f=lb:step:ub;
    
    %calculating the ppm
    ppm=-f/(in.Bo*42.577);
    ppm=ppm+4.65;
    out.ppm = flip(ppm);
end