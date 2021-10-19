% op_CSIFOurierTransform.m
% Brenden Kadota, SunnyBrook Hosptial 2021.
%
% Does spectral and spatial fourier transform on MRSI structure. The fast
% fourier transform is done on the spectral dimension. If the spatial k
% space is cartesian, the fast fourier transform is applied.
% USAGE:
% out = op_CSIFourierTransform(out, spatialFT, spectralFT, isCartesian)
%
%
% input:    in         = Twix object of CSI data
%           spatialFT  = Boolean flag (1 or 0) to compute fourier
%                        transfrom along spatial dimension
%           spectralFT = Same as above but for the spectral dimension
%           isCartesian= boolean flag for fft or slow fourier tranform alon
%                        the spatial dimension
%
% output:   in         = Twix object with new out.specs field of the
%                         fourier transformed data
function [in] = op_CSIFourierTransform(in, k_file, fourier_dimension)
arguments
    in (1,1) struct
    k_file (1,:) char = ""
    fourier_dimension.spatial (1,1) logical {mustHaveSpatial(fourier_dimension.spatial, in)}
    fourier_dimension.spectral (1,1)  logical {mustHaveSpectral(fourier_dimension.spectral, in)} 
end

%if no specs field is made create a new one
if ~isfield(in, 'specs')
    in.specs = in.fids;
end
if(~isfield(fourier_dimension, 'spatial'))
    if(in.flags.spatialFT)
        fourier_dimension.spatial = 0;
    else
        fourier_dimension.spatial = 1;
    end
end
if(~isfield(fourier_dimension, 'spectral'))
    if(in.flags.spectralFT)
        fourier_dimension.spectral = 0;
    else
        fourier_dimension.spectral = 1;
    end
end


if(in.flags.spatialFT == 1 && fourier_dimension.spatial == 1)
    error('already done spatial fourier transform!!!')
end
if(in.flags.spectralFT == 1 && fourier_dimension.spectral == 1)
    error('already done spectral fourier transform!!!')
end
%   fourier transform along the spectral dimension


% spatial dimension fourier transform
if (fourier_dimension.spatial == 1)
    disp('Calculating spatial dimension');
    
    %calculating the x and y coordinates of the spatial dimensions
    xCoordinates = -in.fovX/2+in.deltaX/2:in.deltaX:in.fovX/2-in.deltaX/2;
    xCoordinates = xCoordinates - in.imageOrigin(1);
    yCoordinates = -in.fovY/2+in.deltaY/2:in.deltaY:in.fovY/2-in.deltaY/2;
    yCoordinates = yCoordinates - in.imageOrigin(2);
    
    
    if(k_file == "")
        %applying the fast fourier transform if k space is cartesian
        disp('Applying fast fourier transform');
        
        if(mod(size(in.specs, in.dims.x),2) == 1)
            in.specs = circshift(in.specs,1, in.dims.x);
        end
        in.specs = fftshift(fft(fftshift(in.specs, in.dims.x), [], in.dims.x), in.dims.x);
        if(mod(size(in.specs, in.dims.y),2) == 1)
            in.specs = circshift(in.specs,1, in.dims.y);
        end
        in.specs = fftshift(fft(fftshift(in.specs, in.dims.y), [], in.dims.y), in.dims.y);
    else
        %applying the slow fourier transform if the k space is non
        %cartesian
        %calculating kSpace trajectory for fourier transform (needs to
        %be updated for non cartesian coordinates)
        
        
        [x,y] = meshgrid(xCoordinates, yCoordinates);
        %trajectory of the cartesian trajectory. (needs to be found if
        %slow fourier transform is done on non cartesian data)
        out_traj = [x(:), y(:)];
        k_table = readtable(k_file);
        k_traj = [k_table.K_x, k_table.K_y];
        
        num_TR = max(k_table.TR, [], 'all');
        spatial_points = height(k_table)/num_TR;
        temporal_points = floor(size(in.specs,1)/spatial_points);
        
        %creating fourier transform matrix for the spatial domain
        %fourier transform applied by out = sftOperator*(fids at values at
        %inTraj)
        sftOperator = sft2_Operator(k_traj, out_traj, 0);
        
        %permute so first 3 dimensions are x, y and t
        
        in.specs = permute(in.specs, nonzeros([in.dims.t, in.dims.x, in.dims.y, in.dims.coils, in.dims.averages]));
        shape = size(in.specs);
        %reshape to vector or x,y, and t along the first diimension
        %TODO: figure out a way to do this more elegantly
        if(~ismatrix(in.fids))
            reshaped_specs = reshape(in.specs, in.sz(in.dims.t), in.sz(in.dims.x), in.sz(in.dims.y), []);
            
        else
            reshaped_specs = reshape(in.specs, in.sz(in.dims.t), in.sz(in.dims.x), []);
            
        end
        specs = zeros(temporal_points, size(x, 1), size(x, 2), size(in.specs, 4));
        
        for i = 1:size(in.specs, 4)
            counter = 1;
            for j = 1:spatial_points:in.sz(in.dims.t)
                if(j + spatial_points - 1 < in.sz(in.dims.t))
                    spatial_slice = reshaped_specs(j:j+spatial_points-1, :, :, i);
                    flat_slice = spatial_slice(:);
                    flat_ft = sftOperator*flat_slice;
                    ft = reshape(flat_ft, size(x));
                    specs(counter, :, :, i) = ft;
                    counter = counter + 1;
                end
            end
        end
        in.specs = reshape(specs, [size(specs) shape(4:end)]);
    end
    
    in.flags.spatialFT = 1;
    in.xCoordinates = xCoordinates;
    in.yCoordinates = yCoordinates;
end

if (fourier_dimension.spectral)
    disp('Calculating spectral dimension');
    
    %fourier transform in the spectral domain
    in.specs = fftshift(fft(in.specs,[],in.dims.t),in.dims.t);
    
    %lower bounds of frequency
    lb = (-in.spectralwidth/2)+(in.spectralwidth/(2*in.sz(in.dims.t)));
    %upper bounds of frequency
    ub = (in.spectralwidth/2)-(in.spectralwidth/(2*in.sz(in.dims.t)));
    %frequency step
    step = in.spectralwidth/(size(in.specs,1));
    
    %calculating the frequency
    f=lb:step:ub;
    
    %calculating the ppm
    ppm=-f/(in.Bo*42.577);
    ppm=ppm+4.65;
    in.ppm = flip(ppm);
    in.flags.spectralFT = 1;
end

end

function sft2_Oper = sft2_Operator(InTraj,OutTraj,Ift_flag)
%
% sft2 Perform 2-dimensional slow Fourier transform
%
% This function was written by Bernhard Strasser, April 2018.
%
%
% The function performs a slow Fourier transform in two spatial dimensions by using the definition of the Fourier Transform
% FT(f)(a_j) = sum(k=-N/2;N/2-1) f(x_k)exp(-2pi i x_k a_j)
% (For even N. For odd, change sum slightly)
%
%
% sft2Operator = sft2(A,Ifft_flag)
%
% Input:
% -         Ift_flag                    ...     Flag of 0 or 1 determining if an inverse Fourier transform or normal should be performed.
% -         InputSize                   ...     Size of Input data.
%
% Output:
% -         sft2Operator                ...     Fourier transformation matrix which can be applied to the data by Out = sft2Operator*In
%
%
% Feel free to change/reuse/copy the function.
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: ?
% 0. Preparations

% Define Exponent of exp
if(~Ift_flag)
    Expy = -2*pi*1i;
else
    Expy = 2*pi*1i;
end

% Define Output Size
NOut = size(OutTraj,1);

% 1. FT

sft2_Oper = zeros([size(OutTraj,1) size(InTraj,1)]);
for j=1:NOut
    sft2_Oper(j,:) = exp(Expy*(OutTraj(j,1)*InTraj(:,1) ...
        + OutTraj(j,2)*InTraj(:,2)));
end

if(Ift_flag)
    sft2_Oper = sft2_Oper / (size(InTraj,1));
end

end

function mustHaveSpatial(a, in)
if isfield(in, 'spatialFT') && (a == true && in.spatialFT == 1)
    eidType = 'spatialTransform:noSpatialTransform';
    msgType = 'Spaitial Fourier Transform already done!';
    throwAsCaller(MException(eidType,msgType))
end
end

function mustHaveSpectral(a, in)
if isfield(in, 'spectral') && (a == true && in.spectral == 1)
    eidType = 'spectralTransform:noSpectralTransform';
    msgType = 'Spectral Fourier Transform already done!';
    throwAsCaller(MException(eidType,msgType))
end
end


