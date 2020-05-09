function [OutArray,mask] = EllipticalFilter(InArray,ApplyAlongDims,EllipsoidCoefficients,InputIskSpace_flag)
%
% EllipticalFilter_x_y Apply an elliptical filter to k-space data
%
% This function was written by Bernhard Strasser, July 2012.
%
%
% The function masks the data in k-space, so that k-space values outside of an ellipsoid-shaped mask are set to zero. The mask can be a
% 3d-ellipsoid, or an 2d-ellipse. The equation for the mask is
% mask = {(x,y,z) E R³ | (x/a)² + (y/b)² + (z/c)² <= R²}
% a, b, c, and R can be chosen by the user.
%
%
% [OutArray,mask] = EllipticalFilter_x_y(InArray,ApplyAlongDims,EllipsoidCoefficients,InputIskSpace_flag)
%
% Input: 
% -     InArray                     ...    Input array to which the filter should be applied
% -     ApplyAlongDims              ...    Along these dimensions the filter is applied. If this vector has two elements, a two dimensional 
%                                          Filter is applied. Otherwise, a 3d filter is used.
% -     EllipsoidCoefficients       ...    The values for [a b c R], which determine the shape and size of the ellipsoid. For two dimensional
%                                          Filter, set c = 1;
% -     InputIskSpace_flag          ...    If it is 0, the image gets Fourier transformed to k-space before applying the filter, 
%                                          and transformed back to image domain afterwards
%
% Output:
% -     OutArray                    ...     The filtered/masked output array
% -     mask                        ...     The mask of the filter
%
%
% Feel free to change/reuse/copy the function. 
% If you want to create new versions, don't degrade the options of the function, unless you think the kicked out option is totally useless.
% Easier ways to achieve the same result & improvement of the program or the programming style are always welcome!
% File dependancy: myrepmat

% Further remarks: 



%% 0. Declarations, Preparations, Definitions


% 0.1 Preparations

if(~exist('ApplyAlongDims','var'))
   ApplyAlongDims = [1 2]; 
end
if(~exist('EllipsoidCoefficients','var'))
   EllipsoidCoefficients = [1 1 1 32];                  % This is a ball in x-y-z-space with radius 32
end
EllipsoidCoefficients(EllipsoidCoefficients(1:3) == 0) = 1;  % There are no zeros allowed, because through these values gets divided. 
if(~exist('InputIskSpace_flag','var'))
   InputIskSpace_flag = true; 
end 


% 0.2 Declarations


% 0.3 Definitions
    
OutArray = InArray;
kSpaceCenter = [0 0 0];
for Dim = 1:numel(ApplyAlongDims)
    kSpaceCenter(Dim) = floor(size(OutArray,ApplyAlongDims(Dim))/2) + 1;
end

% Due to small numerical errors, inputting EllipsoidCoefficients = [1 1 1 R] gives different results than EllipsoidCoefficients = [R R R 1], which should be equivalent. This Epsilon makes them same again.
Epsilon = 0.000001;         

 





%% 1. FFT to k-space

if(~InputIskSpace_flag)

    for filter_dim = ApplyAlongDims
        OutArray = ifftshift(OutArray,filter_dim);
        OutArray = ifft(OutArray,[],filter_dim);
        OutArray = fftshift(OutArray,filter_dim);
    end  

end




%% 2. Compute Mask


if(numel(ApplyAlongDims) == 2)
    [x_grid,y_grid,z_grid] = ndgrid(1:size(OutArray,ApplyAlongDims(1)), 1:size(OutArray,ApplyAlongDims(2)),0);
else
    [x_grid,y_grid,z_grid] = ndgrid(1:size(OutArray,ApplyAlongDims(1)), 1:size(OutArray,ApplyAlongDims(2)), 1:size(OutArray,ApplyAlongDims(3)));   
end

mask = ( (x_grid - kSpaceCenter(1)) / EllipsoidCoefficients(1) ).^2 + ( (y_grid - kSpaceCenter(2)) / EllipsoidCoefficients(2) ).^2 + ( (z_grid - kSpaceCenter(3)) / EllipsoidCoefficients(3) ).^2 <= EllipsoidCoefficients(4)^2+Epsilon;
mask = myrepmat(mask,size(OutArray));




%% 3. Apply Mask

OutArray = mask .* OutArray;

if(nargout < 2)
    clear mask
end




%% 4. FFT to ImageSpace


if(~InputIskSpace_flag)
    
    for filter_dim = ApplyAlongDims
        OutArray = ifftshift(OutArray,filter_dim);
        OutArray = fft(OutArray,[],filter_dim);
        OutArray = fftshift(OutArray,filter_dim);
    end  
    
end




%% 5. Postparations







