
% rf2d_shapes.m
% Jason Rock, Jamie Near, Elena Osipyan, 2025.
%
% USAGE:
%   img = rf2d_shapes(instance)
%
% DESCRIPTION:
%   Generate a synthetic 2-D target pattern (binary or grayscale) for
%   2DRF excitation testing/visualization. The pattern is selected by
%   'instance' (enumerated below). Except for the Gaussian (case 3,
%   1024×1024), all outputs are 256×256. Outputs are of class double;
%   most cases are in [0,1].
%
% INSTANCES:
%   1  Square (64×64 block centered in 256×256; value = 1 inside).
%   2  External image 'triangle.png' → grayscale → rescaled to [0,1] →
%      resized to 256×256. triangle.png can be changed to any desired
%      excitation image.

%   3  2-D isotropic Gaussian (σ = 120 px) on a 1024×1024 grid with peak
%      ≈ 1/(2πσ^2). 
%   4  Isosceles triangle (binary) blurred with imgaussfilt(·, 4).
%   5  Filled circle (radius = 256/6) centered at (128,128).
%   6  Annular “lipid layer”: r ∈ (r0, 2r0), r0 = 256/6, blurred σ = 4.
%   7  External mask 'mask_only.png' → grayscale in [0,1].
%   9  Shifted square: case 1 translated by (+30, +30) px.
%   10 Square cross: orthogonal bars with arm_length = 64, arm_width = 32.
%
% INPUTS:
%   instance   - Integer selector ∈ {1,2,3,4,5,6,7,9,10}.
%
% OUTPUTS:
%   image      - Double matrix:
%                  • 256×256 for all cases except instance = 3
%                  • 1024×1024 for instance = 3 (Gaussian)
%                Values are binary or grayscale as described above.
%
% DEPENDENCIES:
%   Requires Image Processing Toolbox functions: imgaussfilt, imresize.
%   Files required on path for certain cases:
%     • 'circle.png' (instance 2)
%     • 'mask_only.png' (instance 7)
%
% NOTES:
%   • Unsupported instance values trigger: error('Invalid case_id').

function[image] = rf2d_shapes(instance)
    
        if instance == 3 % Gaussian
            image = zeros(1024, 1024);
            sig = 120;
            for i = 1:1024
                for j = 1:1024
                    image(i, j) = 1/(2*pi*sig^2) * exp(-((i - 512)^2 + (j - 512)^2) / (2 * sig^2));
                end
            end
        end

         if instance == 1 % Square
            image = zeros(256, 256);
            for i = 96:160
                for j = 96:160
                    image(i, j) = 1;
                end
            end
         end

        if instance == 2 % Load and resize image
            image_filename = 'triangle.png';  % Change this to your desired image
            img = imread(image_filename);
            if size(img, 3) == 3
                img = mean(img, 3); % Convert RGB to grayscale
                img = (img - min(img(:)))/(max(img(:)) - min(img(:)));
                %img = img/max(abs(img));
            end
            image = imresize(img, [256, 256]);
        end
        if instance== 4 % Triangle
            image = zeros(256, 256);
            counter = 1;
            for i = 101:154
                for j = 1:(2*counter)
                    col = 128 - counter + j;
                    if col >= 1 && col <= 256
                        image(i, col) = 1;
                    end
                end
                counter = counter + 1;
            end
            image = imgaussfilt(image, 4);
        end
        if instance == 5 % Circle
            image = zeros(256, 256);
            radius = 256 / 6;
            for i = 1:256
                for j = 1:256
                    if (i - 128)^2 + (j - 128)^2 < radius^2
                        image(i, j) = 1;
                    end
                end
            end
        end
        if instance == 6 % Lipid layer
            image = zeros(256, 256);
            radius = 256 / 6;
            radius2 = 2 * radius;
            for i = 1:256
                for j = 1:256
                    r2 = (i - 128)^2 + (j - 128)^2;
                    if r2 > radius^2 && r2 < radius2^2
                        image(i, j) = 1;
                    end
                end
            end
            image = imgaussfilt(image, 4);
        end
        if instance == 7 % Mask image
            image_filename = 'mask_only.png';
            img = imread(image_filename);
            if size(img, 3) == 3
                img = mean(img, 3);
            end
            image = double(img) / max(img(:));
        end

        if instance == 9 % Shifted square
            image = zeros(256, 256);
            shift_sq_x = 30;
            shift_sq_y = 30;
            for i = (96 - shift_sq_x):(160 - shift_sq_x)
                for j = (96 - shift_sq_y):(160 - shift_sq_y)
                    new_i = i;
                    new_j = j;
                    if new_i >= 1 && new_j >= 1 && new_i <= 256 && new_j <= 256
                        image(new_i, new_j) = 1;
                    end
                end
            end
        end 

        if instance == 10 % Square cross
            image = zeros(256, 256);
            arm_length = 64;
            arm_width = 32;
            center = 128;
            image(center - arm_width/2 + 1:center + arm_width/2, center - arm_length + 1:center + arm_length) = 1;
            image(center - arm_length + 1:center + arm_length, center - arm_width/2 + 1:center + arm_width/2) = 1;
        end
        if ~ismember(instance, [1,2,3,4,5,6,7,9,10])
            error('Invalid case_id');
        end

end