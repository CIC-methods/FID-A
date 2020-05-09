%CSItoTwix.m
%converts a CSI twix to a single voxel twix file.
%allows for the single voxel FID-A porcessing tools to be used

function out = CSItoMRS(in, xCoordinate, yCoordinate)
    out = in;
  
    if(in.dims.coils ~= 0)
        out.specs = in.specs(:,:,xCoordinate, yCoordinate);
    else
        out.specs = in.specs(:,xCoordinate, yCoordinate);
    end
    
    out.fids = ifft(ifftshift(out.specs, in.dims.t), [], in.dims.t);
    out.sz = size(out.fids);
    out.dims.x = 0;
    out.dims.y = 0;
    
end
