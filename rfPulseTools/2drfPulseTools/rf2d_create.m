% rf2d_create.m
% Jason Rock, Elena Osipyan and Jamie Near, 2025.
%
% USAGE:
% rf2d = rf2d_create.makepulse(params)
%
% DESCRIPTION:
% This class provides functions to design, generate, and export 2D RF (2DRF) 
% pulses with corresponding gradient trajectories. It implements spiral k-space 
% trajectories with gradient and slew rate constraints, Gaussian k-space filtering, 
% B1 peak calibration, and export of Siemens-compatible RF and gradient text files. 
% The class also supports visualization of excitation targets, k-space trajectories, 
% and the resulting RF pulse amplitude/phase profiles.
%
% INPUTS (params struct):
% params.MAX_GRAD     = Maximum gradient amplitude [T/m] (default = 0.08 T/m)
% params.MAX_SLEW     = Maximum gradient slew rate [T/m/ms] (default = 0.2 T/m/ms)
% params.RF_RASTER    = RF raster time [s] (default = 5e-6 s)
% params.GRAD_RASTER  = Gradient raster time [s] (default = 1e-5 s)
% params.type         = Type of RF pulse (default = '2drf', anything else yields exception)
% params.shape        = Ex citation shape (defined in rf2d_shapes.m). Example: shapes(1) is
% a square.

% params.FOV          = Field of view [cm] (default = 12 cm)
% params.tp           = Total pulse duration [s] (default = 0.015 s)
% params.k_us         = Transmit k-space trajectory undersampling rate (default = 4)
% params.flip_angle   = Flip angle [rad] (default = 10 deg / 0.174533 rad)
% params.ppm          = Off-resonance/chemical-shift for simulation [ppm] (on-resonance default  = 0)
% params.voxels       = Number of spatial samples per axis for design/Bloch sim grid (default = 128)

% params.calibration  = Calibration for calculating the B1 peak. Options:
% 'center' (gets the center pixel of the image), 'brightest' (gets the
% brightest pixel of the image), 'custom' (lets user choose custom pixel of
% image)
%
% OUTPUTS:
% rf2d struct containing:
%   .waveform   = RF pulse [amplitude, phase, 0]
%   .type       = '2drf'
%   .gradients  = Gradient waveforms [gx, gy, gz]
%   .b1_peak    = Calibrated B1 peak [Tesla]
%   .flip_angle = Flip angle [radians]
%   .tp         = Total pulse duration [s]
%
% FILE OUTPUTS:
% 2DRFpulse.txt   = RF amplitude and phase
% 2DGpulse.txt    = Gradient waveforms
% rf2d.mat        = MATLAB struct containing pulse and gradient information


classdef rf2d_create
    methods(Static)
        function[res] =  shift(b1, x0, k)
    res = b1 * exp (-1i * dot(x0, k));
end

function val = first_zero_crossing(arr)
    val = NaN;  % Default return value
    if length(arr) < 2
        return;
    end
        arr = arr(isfinite(arr));
    if length(arr) < 2
        return;
    end
        for i = 1:(length(arr) - 1)
        if (arr(i) >= 0 && arr(i+1) < 0) || (arr(i) < 0 && arr(i+1) >= 0)
            val = i;
            return;
        end
    end
    
    % If no zero crossing found, return a reasonable default
    val = min(20, floor(length(arr)/4));
end

% ===== Helpers (boolean-returning) =====

function [ok, b1_peak] = attempt_center(image_in, flip_angle, amp, phase, tp)
    % center is (0,0); reject if center is black (no data)
    b1_peak = [];
    % Get dimensions 
    [rows, cols] = size(image_in); % Compute center pixel coordinates
    cy = round(rows/2);
    cx = round(cols/2);

    isBlack = rf2d_create.check_coord(image_in, cx, cy);  % true => black/invalid
    if isBlack
        ok = false; return;
    end
    offset_x = 0; offset_y = 0;
    [ok, b1_peak] = rf2d_create.compute_b1_peak(offset_x, offset_y, flip_angle, amp, phase, tp);
end

function [ok, b1_peak] = attempt_brightest(image_in, flip_angle, amp, phase, tp, x_axis, y_axis)
    b1_peak = [];
    try
        [cx, cy] = rf2d_create.brightest(image_in);
         isBlack = rf2d_create.check_coord(image_in, cx, cy);  % true => black/invalid
         if isBlack
             ok = false; return;
         end
    catch
        ok = false; return;
    end
    [ok, b1_peak] = rf2d_create.compute_b1_peak(cx/1e6, cy/1e6, flip_angle, amp, phase, tp);
end

function [ok, b1_peak] = attempt_custom(image_in, flip_angle, amp, phase, tp, x_axis, y_axis)
    b1_peak = [];
    try
        fprintf('Click your preferred pixel for B1 calibration.\n');
        [cx, cy] = rf2d_create.pick_point(image_in);
         isBlack = rf2d_create.check_coord(image_in, cx, cy);  % true => black/invalid
         if isBlack
             ok = false; return;
         end
    catch
        ok = false; return;
    end
    [ok, b1_peak] = rf2d_create.compute_b1_peak(cx/1e6, cy/1e6, flip_angle, amp, phase, tp);
end

function [ok, b1_peak] = compute_b1_peak(offset_x, offset_y, flip_angle, amp, phase, tp)
    % Core boolean computation: no throws, returns ok=false on any invalidity
    b1_peak = [];
    ok = false;

    % Guard: amp must be nonzero to normalize
    if isempty(amp) || all(~isfinite(amp)) || max(abs(amp)) == 0
        return;
    end

    % Run calibration sweep via your existing API
    try
        [b1_peaks, M_z_cal] = rf2d_create.b1_peak_calc( ...
            flip_angle, amp, phase, tp, offset_x, offset_y);
    catch
        return;
    end

    % Find first zero crossing
    try
        idx = rf2d_create.find_zero_crossings(M_z_cal, flip_angle);
    catch
        idx = [];
    end
    if isempty(idx) || any(~isfinite(idx))
        return;
    end

    cand = b1_peaks(idx(1));
    if isempty(cand) || ~isfinite(cand) || isnan(cand)
        return;
    end

    b1_peak = cand;
    ok = true;
end


function[B1_shifted] = modify_rf_pulse_for_shift(B1, kx, ky, shift_x, shift_y)
    % RF pulse can be modified to account for object shift, handling phase wrapping.
    phase_shift = exp(-1i * (kx * shift_x + ky * shift_y));
    phase_shift_wrapped = angle(phase_shift);  % Wraps the phase to -pi to pi
    phase_shift_corrected = exp(1i * phase_shift_wrapped);

    B1_shifted = B1 * phase_shift_corrected;
end 

function idx = find_zero_crossings(arr, flip_angle)
        idx = find(diff(sign(arr - cos(flip_angle))));
end

function isBlack = check_coord(image_in, x, y)
% CHECK_COORD  Return true if a given pixel in a double image is zero (black).
%   isBlack = check_coord(image_in, x, y)
%
%   Inputs:
%     image_in : 2D or 3D numeric array (double image)
%     x        : column index (cx)
%     y        : row index (cy)
%   Output:
%     isBlack  : logical true/false

    % If RGB double image, convert to grayscale
    if ndims(image_in) == 3 && size(image_in,3) == 3
        img = mean(image_in, 3);  % simple grayscale conversion
    elseif ndims(image_in) == 2
        img = image_in;
    else
        error('Unsupported input: must be 2D grayscale or 3D RGB double image.');
    end

    % Ensure indices are valid integers inside the image
    [rows, cols] = size(img);
    if x < 1 || x > cols || y < 1 || y > rows
        error('Coordinates (x=%d, y=%d) out of image bounds.', x, y);
    end

    % Logical check (allow small tolerance for floating-point)
    tol = 1e-12;
    isBlack = (abs(img(round(y), round(x))) < tol);

end
function [cx, cy, isBlack] = pick_point(imgIn)
% PICK_POINT  Let the user click a location on the image. Returns (cx, cy).
%   [cx, cy] = pick_point(imgIn)
%   - imgIn: image array or filename (grayscale or RGB)
%   - cx, cy: clicked coordinates (x=column, y=row), double precision

    % Load image if filename given
    if ischar(imgIn) || isstring(imgIn), img = imread(imgIn); else, img = imgIn; end

    % Convert per your rule
    if ndims(img) == 2
        I = im2gray(img);
    elseif ndims(img) == 3 && size(img,3) == 3
        I = rgb2gray(img);
    else
        error('Unsupported image format: use 2D grayscale or 3D RGB.');
    end

    % Show image and get one click
    hFig = figure('Name','Click a point'); 
    imshow(I, []); axis image; hold on;
    title('Click your preferred point (press Enter/Esc to cancel)');
    try
        [x, y, btn] = ginput(1);   % one click
    catch
        close(hFig); error('Click canceled.');
    end
    if isempty(x) || isempty(y) || ~isscalar(btn)
        close(hFig); error('Click canceled.');
    end

    % Clip to image bounds and return
    [H, W] = size(I);
    cx = min(max(x, 1), W);
    cy = min(max(y, 1), H);
    % Visual feedback
    plot(cx, cy, 'r+', 'MarkerSize', 15, 'LineWidth', 2);
    drawnow;
end

function result = mod_custom(p)
    result = p - 2 * pi * floor(p / (2 * pi));
end

function[ft] =  calculate_2dft(image_in)
   image_in = im2gray(image_in); % Convert to grayscale if RGB
   ft = ifftshift(image_in);
   ft = fft2(ft);
   ft = fftshift(ft);
end

function[b1_peaks, M_z_cal] = b1_peak_calc(flip_angle, amp, phase, tp, offset_x, offset_y)
    b1_min = 1e-8;
    b1_max = 25e-5;
    disp(['*** Flip angle: ', num2str(rad2deg(flip_angle)), ' deg'])

    dt = tp / length(amp); % tp should be defined earlier
    M_z_cal = zeros(1, 500);
    b1_peaks = linspace(b1_min, b1_max, 500);
    gamma = 267.522*10^6;

    for i = 1:length(b1_peaks)
        b1 = (amp / max(amp)) * b1_peaks(i);  % Normalized B1 amplitude
        M_temp = [0;0;1];
        for m = 1:length(b1)
            phi1 = phase(m);  % phase in radians
            B1_eff = sqrt(b1(m)^2 + (offset_x + offset_y)^2);
            alpha1 = atan2(b1(m), (offset_x + offset_y));
            theta1 = gamma * B1_eff * dt;

            M_temp = rf2d_bloch_sim.rot_z(M_temp, phi1);
            M_temp = rf2d_bloch_sim.rot_y(M_temp, alpha1);
            M_temp = rf2d_bloch_sim.rot_z(M_temp, theta1);
            M_temp = rf2d_bloch_sim.rot_y(M_temp, -alpha1);
            M_temp = rf2d_bloch_sim.rot_z(M_temp, -phi1);
        end

        M_z_cal(i) = M_temp(3);  % Store z-component
    end
end

function[var] = G(T, t, N, k_max)  % grad trajectory -  Rf pulse =  Rf(interpolation) * |G|
    var = sqrt((2*N*pi*(1 - t/T)).^2 + 1);
end

function kx = k_x(t, A, n, T, variable_density)
% k_x - Computes x-component of a spiral k-space trajectory
%
% Inputs:
%   t - time array
%   A - amplitude (k_max)
%   n - number of spiral turns
%   T - total RF pulse durationmath
%   variable_density - 1x2 logical+numeric array: [true/false, density_param]
%
% Output:
%   kx - k-space x-component trajectory

    if nargin < 5
        variable_density = [false, 0];  % default: uniform density
    end

    is_vd = variable_density(1);        % true or false
    density_param = variable_density(2);

    if is_vd
        kx = A * (1 - t/T) .* cos(2 * pi * t * n / T) .* exp(density_param * t / T);
    else
        kx = A * (1 - t/T) .* cos(2 * pi * t * n / T);
    end

end

function ky = k_y(t, A, n, T, variable_density)
% k_y - Computes y-component of a spiral k-space trajectory
%
% Inputs:
%   t - time array
%   A - amplitude (k_max)
%   n - number of spiral turns
%   T - total RF pulse duration
%   variable_density - 1x2 logical+numeric array: [true/false, density_param]
%
% Output:
%   ky - k-space y-component trajectory

    if nargin < 5
        variable_density = [false, 0];  % default: uniform density
    end

    is_vd = variable_density(1);        % true or false
    density_param = variable_density(2);

    if is_vd
        ky = A * (1 - t/T) .* sin(2 * pi * t * n / T) .* exp(density_param * t / T);
    else
        ky = A * (1 - t/T) .* sin(2 * pi * t * n / T);
    end

end

function[conv] = G_x(kx, dt)
    % gamma in rad/(s*T) = 267.522e6 rad/(s*T)
    gamma = 267.522e6;
    
    % Convert k-space trajectory to gradient
    % kx is in cm^-1, need to convert to m^-1 (multiply by 100)
    % Then divide by gamma to get gradient in T/m
    %conv = (gradient(kx * 1000)/dt) / gamma / pi;  % T/m

    conv =  1/gamma / (pi) * gradient(1000*kx, dt);
    
    % Convert to mT/m for display (multiply by 1000)
    %conv = conv * 1000;  % mT/m
end

function[conv] = G_y(ky, dt)
    % Fixed gradient calculation for Y  
    % gamma in rad/(s*T) = 267.522e6 rad/(s*T)
    gamma = 267.522e6;
    
    % Convert k-space trajectory to gradient
    % ky is in cm^-1, need to convert to m^-1 (multiply by 100)
    % Then divide by gamma to get gradient in T/m
    %conv = (gradient(ky * 1000)/dt) / gamma / pi;  % T/m
    conv =  1/gamma / (pi) * gradient(1000*ky, dt);
    % Convert to mT/m for display (multiply by 1000)
    %conv = conv * 1000;  % mT/m
end

function [cx, cy] = brightest(imgIn)
% BRIGHTEST: Mark the brightest location in an image with a box.
%   [cx, cy, peakVal] = brightest (imgIn, boxSizePx)
%   - imgIn: image array or filename
%   - boxSizePx: optional, side length of the box in pixels (auto if omitted)
%   Returns centroid (cx,cy) in image coordinates (x=col, y=row) and peak value.

    % --- Load
    if ischar(imgIn) || isstring(imgIn), img = imread(imgIn); else, img = imgIn; end

    % --- Convert per your rule
    if ndims(img) == 2
        I = im2gray(img);
    elseif ndims(img) == 3 && size(img,3) == 3
        I = rgb2gray(img);
    else
        error('Unsupported image format. Use 2D grayscale or 3D RGB.');
    end

    % --- Safe double image in [0,1], handle NaNs
    I = im2double(I);
    I(isnan(I)) = -inf;  % exclude NaNs from the max

    % --- Find brightest location(s)
    peakVal = max(I(:));
    if ~isfinite(peakVal)
        error('Image has no finite pixels.');
    end
    [rows, cols] = find(I == peakVal);   % all ties
    cy = mean(rows);                      % average if multiple maxima
    cx = mean(cols);
end


function[absolute_values, phase_values] = write_b1(amp, pha)
   absolute_values = amp(:);
   phase_values = pha(:);
   lenRF = length(absolute_values);
   file_path = fullfile('2DRFpulse.txt');
   fid = fopen(file_path, 'w');
   fprintf(fid, '%d\n', lenRF);  % write the header line
   for i = 1:lenRF
     fprintf(fid, '%.6f\t%.6f\n', absolute_values(i), phase_values(i));  % fixed 6 decimal places
   end
   fclose(fid);
  
   fprintf("Complex B1 has been written to %s\n", file_path);
end

function[gx,gy] = write_grad(gx, gy)
   % Siemens requires gradients to be normalized
   gx = gx(:);  % Ensure column vector
   gy = gy(:);  % Ensure column vector
  
   gx_max_raw = max(abs(gx));
   gy_max_raw = max(abs(gy));
   gx = gx / gx_max_raw;
   gy = gy / gy_max_raw;
   gz = zeros(length(gx), 1);
   % For Siemens requirements - ensure last elements are zero
   if length(gx) > 0
       gx(end) = 0.00;
       gy(end) = 0.00;
       gz(end) = 0.00;
   end
   Gx_max = round(gx_max_raw, 3);  % in mT/m
   Gy_max = round(gy_max_raw, 3);
   % Preparing the variables for header
   lenGrad = int32(length(gx));
   header = [Gx_max, Gy_max, 0];
   header_row = reshape(header, 1, numel(header));
   % Main data table
   data = [gx(:), gy(:), gz(:)];
   %Define the file path
   file_path = fullfile('2DGpulse.txt');
   % % Write lenGrad as the first row
    fid = fopen(file_path, 'w');  % Overwrite to start clean
    fprintf(fid, '%d\n', lenGrad);
   %
   % % Write header row
   fprintf(fid, '%.6f\t%.6f\t%.6f\n', header_row);
  
   % Write data
   fprintf(fid, '%.6f\t%.6f\t%.6f\n', data');
   fclose(fid);
  
   fprintf("Gradient data written to %s\n", file_path);
      end


  function[gx,gy] = write_files(amp, pha, b1, kx, ky, dt, withRamp, gx, gy, T, RF_RASTER, GRAD_RASTER)

   RF_samples = round(T / RF_RASTER);
   Grad_samples = round(T / GRAD_RASTER);
   % Ensure inputs are column vectors
   amp = amp(:);
   pha = pha(:); 
   gx = gx(:)*1000;
   gy = gy(:)*1000;
   len_padding =  rf2d_create.first_zero_crossing(gy);
   if withRamp
       len_padding_global = len_padding;
       b1_new_amp = zeros(1, len_padding + length(amp));
       b1_new_amp(len_padding+1:end) = amp;

       b1_new_pha = zeros(1, len_padding + length(pha));
       b1_new_pha(len_padding+1:end) = pha;

       % Resampling
       RF_samples = round(T / RF_RASTER);
       b1_new_amp = b1_new_amp / max(b1_new_amp);
       rf2d_create.write_b1(b1_new_amp, b1_new_pha);
   
       gx_new = zeros(len_padding + length(gx), 1);
       gx_new(len_padding + 1:end) = gx;
       gy_new = zeros(len_padding + length(gy), 1);
      
       % Create ramp part - ensure we don't exceed array bounds
       ramp_length = min(len_padding, length(gy));
       if ramp_length > 0
           ramp_part = flipud(gy(1:ramp_length));
           gy_new(1:ramp_length) = ramp_part;
       end
      
       % Fill the rest of gy_new with original gy values
       remaining_start = len_padding + 1;
       remaining_end = length(gy_new);
       source_start = 1;
       source_end = min(length(gy), remaining_end - remaining_start + 1);
      
       if source_end > 0 && remaining_start <= remaining_end
           gy_new(remaining_start:remaining_start + source_end - 1) = gy(source_start:source_end);
       end
     
       gx = resample(gx_new, (Grad_samples + floor(len_padding * Grad_samples / (T/RF_RASTER))), length(gx_new));
       gy = resample(gy_new, (Grad_samples + floor(len_padding * Grad_samples / (T/RF_RASTER))), length(gy_new));
      
       rf2d_create.write_grad(gx, gy);


   else
       % Non-ramping version
       Grad_samples = int32(round(T / GRAD_RASTER));
       % === Resample b1 amplitude ===
       
       gx = resample(gx, Grad_samples, length(gx));
       gy = resample(gy, Grad_samples, length(gy));

       pha(pha > 2*pi) = 2*pi;
       pha(pha < 0) = 0;

       %plot(pha);
       %title('Pha After Write Files');
       rf2d_create.write_b1(amp, pha);
       rf2d_create.write_grad(gx, gy);  % Non-ramping
   end
end
function[rf2d] = makepulse(params)
    clear;
    tic;
    load params
    image_in = params.shape;
    image_in = fliplr(image_in);
    FOV = params.FOV; % cm. FOV = dr * Nx
    dr = (FOV / (2 * size(image_in, 2)));
    k_us = params.k_us;
    MAX_GRAD = params.MAX_GRAD; %T/m, could be much higher. Up for the user to change! 
    MAX_SLEW = params.MAX_SLEW; % T/m/ms figure out, will affect spiral density N
    % raster times - minimum hardware discretization
    RF_RASTER = params.RF_RASTER;
    GRAD_RASTER = params.GRAD_RASTER;
    T = params.tp; % seconds - total pulse duration
% high spatial frequencies are hard to capture - multiply k space by gaussian.
% if we don't filter signal -  lots of ripples in image. hence we gauss blur the k space
    save_it = true;
    withRamp = true;
    b1_peak_algo = true;
    if ~strcmp(params.type, '2drf')
    ME = MException('ForceStop:ExecutionTerminated', ...
        'You are not creating a 2DRF Pulse! Program stopped.');
    throw(ME);
    end
t = linspace(0, T, T/RF_RASTER); 

dt = t(2) - t(1);
image_in = fliplr(image_in);

ft = rf2d_create.calculate_2dft(image_in);
[rows, cols] = size(ft);

%Create a meshgrid of the same size as the Fourier transform
x1 = linspace(-floor(cols/2), floor(cols/2), cols);
y1 = linspace(-floor(rows/2), floor(rows/2), rows);
[X1, Y1] = meshgrid(x1, y1);


% the k space trajectory is calculated through the parametric equation. the gradient of the spiral
no_amp_restriction = false;
% gets you a set of gradient. there is a max amplitude, the rate of change of the gradient is also limited.
% always check that the slope of gradient is the max
% ex. n=26, kmax is prop to res
variable_density = [false, 1];

if (no_amp_restriction)
    N = 26; % Arbitrary number
    k_max = 1 / (2 * dr); % cm^-1
    % k_max = N / FOV
    % num of spirals too low while kmax is high = aliasing

    kx_axis = linspace(-k_max, k_max, size(image_in, 2));
    ky_axis = linspace(-k_max, k_max, size(image_in, 1));
    [KX, KY] = meshgrid(kx_axis, ky_axis);
    kx = rf2d_create.k_x(t, k_max / k_us, N, T, variable_density);
    ky = rf2d_create.k_y(t, k_max / k_us, N, T, variable_density);

else
    % restriction with grad amplitude and slew
    % can either restrict number of spirals and k_max for fixed FOV

    restrict_N = true;

    if (restrict_N)

        % for brain/circle, effective FOV = 1
       max_of_max_slew = 0.95;  % 80% of MAX SLEW
        % N_list = np.linspace(30,2, 30-2+1)
       N_list = 30:-0.2:1.8; % list of spirals, increments of .2

       k_max = 1 / (2 * dr);  % cm^-1
       kx_axis = linspace(-k_max, k_max, size(image_in, 2));
       ky_axis = linspace(-k_max, k_max, size(image_in, 1));
       [KX, KY] = meshgrid(kx_axis, ky_axis);


       for i = 1:(length(N_list))
           kx = rf2d_create.k_x(t, k_max/k_us, N_list(i), T);
           ky = rf2d_create.k_y(t, k_max/k_us, N_list(i), T);
            if (max(rf2d_create.G_x(kx, dt)) < MAX_GRAD && max(rf2d_create.G_y(ky, dt)) < MAX_GRAD)
                %grad_gx = gradient(G_x(kx, dt), dt);
                if (max(abs(gradient(rf2d_create.G_x(kx, dt), dt)))/1000) < MAX_SLEW*max_of_max_slew && max(abs(gradient(rf2d_create.G_y(ky, dt), dt))/1000) < MAX_SLEW*max_of_max_slew
                    N = N_list(i);
                    %plot(kx); hold on; plot(ky);
                    fprintf("restricted k_x, k_y found with: N = %d\n", N);
                    break
                end
            end
       end 
    end 
    
end


% Define the standard deviation for the Gaussian
% 0.03 more blurry than 0.05
% Adjust this value as needed, here scaled with the size
sigma = 0.08 * min(rows, cols);

% Create the Gaussian filter
gaussian_filter = exp(-(X1.^2 + Y1.^2) / (2 * sigma^2));

% Apply the Gaussian filter to the Fourier transform
ft = ft .* gaussian_filter;
% --- Setup ---
figure;
% Manually move title upward by increasing the y-position
sgtitle('RF Pulse Design');

% --- Parameters ---
kx = kx(:);
ky = ky(:);

% --- Axis vectors (assumes square image) ---
x_axis = linspace(-dr * size(image_in, 2), dr * size(image_in, 2), size(image_in, 2));
y_axis = x_axis;
[X, Y] = meshgrid(x_axis, y_axis);

% --- Left Plot: Original Image ---
nexttile;
imagesc(x_axis, y_axis, rot90(image_in, 2));
axis equal tight;
colormap gray;
xlabel('[cm]');
ylabel('[cm]');
title('Target Excitation Region');

% --- Right Plot: k-space & Trajectory Overlay ---
nexttile;

% Convert k-space coordinates to pixel units
kx_pixels = ((kx / k_us) * (dr * size(image_in, 1) / 2)); 
ky_pixels = ((ky / k_us) * (dr * size(image_in, 2)/ 2)); 

% Determine bounding box of the k-space trajectory
ft_logabs = log(abs(ft));
% Display Fourier transform of image within trajectory bounds
imagesc([min(x_axis(1,:)), max(x_axis(1,:))], [min(y_axis(1,:)), max(y_axis(1,:))], ft_logabs);

% Display Fourier transform of image within trajectory bounds
hold on;
plot(kx_pixels/3, ky_pixels/3, 'r-', 'LineWidth', 2);
xlim([min(x_axis(1,:)) max(x_axis(1,:))]);
axis equal tight;
colormap gray;
alpha(0.8);  % Make FT image semi-transparent
xlabel('$k_x$ [cm$^{-1}$]', 'Interpreter', 'latex');
ylabel('$k_y$ [cm$^{-1}$]', 'Interpreter', 'latex');
title('k-space trajectory over $\mathcal{FT}$', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');
set(gcf, 'Color', 'w');


grad_gx = gradient(rf2d_create.G_x(kx, dt), dt);
max_slew = max(abs(grad_gx)) / 1000;  % T/m/ms
fprintf("Slew rate MAX for G_x: %.3f T/m/ms\n", max_slew);
% Print info
%fprintf('Slew rate MAX for G_x: %.3f T/m/ms\n', max(abs(grad_gx)) / 1000);
fprintf('Slew rate MAX: 0.200 T/m/ms\n');
fprintf('Gradient frequency: N / T = %.2f Hz\n', N/T);
fprintf('Forbidden gradient frequency ranges:\n');
fprintf('    * 1140 Hz +/- 110 Hz\n');
fprintf('    * 590 Hz +/- 50 Hz\n');

if (N/T <= 1250 && N/T >= 1030) || ...
   (N/T <= 640  && N/T >= 540)
    fprintf('You are within a forbidden gradient frequency range!\n');
end

% Reconstruct matching k-space grid for interpolation
kx_axis = linspace(-k_max, k_max, size(ft, 2));  % 64 points
ky_axis = linspace(-k_max, k_max, size(ft, 1));  % 64 points
[KX, KY] = meshgrid(kx_axis, ky_axis);           % KX, KY are 64x64

% Flatten data for griddata
points = [KX(:), KY(:)];   % (4096, 2)
values = ft(:);            % (4096, 1)

% Confirm shape match
assert(size(points,1) == length(values), 'griddata input mismatch');

%b1 = griddata(points(:,1), points(:,2), values, kx(:), ky(:), 'linear');

% Apply G modulation
G_vector =  rf2d_create.G(T, t, N, k_max);
% Use 1D query points
b1 = griddata(points(:,1), points(:,2), values, kx(:), ky(:), 'cubic');

%b1_unfixed = griddata(points(:,1), points(:,2), values, kx(:), ky(:), 'nearest');
%b1_unfixed = b1_unfixed .* G_vector;

b1 = b1(:, 1);
time_axis = linspace(0, T, length(b1));

% Plot b1 against time_axis (converted to ms)


G_vector = reshape(rf2d_create.G(T, t, N, k_max), [], 1);
if isvector(G_vector) && isequal(numel(G_vector), numel(b1))
    b1 = b1 .* G_vector;
elseif isscalar(G_vector)
    b1 = b1 * G_vector;
else
    error('Size mismatch between b1 and G(t, T, N, k_max)');
end


% ----
len_padding_global = 0;

gx = rf2d_create.G_x(kx, dt); % good here
gy = rf2d_create.G_y(ky, dt);
amp = abs(b1) / max(abs(b1));
pha =  rf2d_create.mod_custom(unwrap(angle(b1)));

if (save_it)
    
    [g_x, g_y] =  rf2d_create.write_files(amp, pha, b1, kx, ky, dt, withRamp, gx, gy, T, RF_RASTER, GRAD_RASTER);
end 

% Open and read the file
pulsedata = readmatrix('2DRFpulse.txt');
amp_data = pulsedata(:,1);
phase_data = pulsedata(:,2); 

hold off;

tp = T;

flip_angle = params.flip_angle;
phase = phase_data;
amp = amp_data;
b1_peak = 8e-6;
b1 = ((amp) / max(amp))* b1_peak;
b1 = reshape(b1.', 1, []);  % transpose first, then reshape into a row vector
offset_x = 0;
offset_y = 0;
% ===== Entry point =====
if b1_peak_algo
    method = lower(strtrim(params.calibration));  % 'center' | 'custom' | 'brightest'
    valid = any(strcmp(method, {'center','custom','brightest'}));
    assert(valid, 'Invalid params.calibration: %s', params.calibration);

    done = false;
    max_loops = 3;   % upper bound to avoid accidental infinite loops
    loops = 0;

    while ~done
        loops = loops + 1;
        if loops > max_loops
            error('B1Peak:AllFailed','B1 Peak cannot be found based on provided excitation shape. Switch excitation shape.');
        end

        switch method
            case 'center'
                [ok, b1_peak] = rf2d_create.attempt_center(image_in, flip_angle, amp, phase, tp);
                if ok
                    done = true;
                else
                    resp = input('B1 Calibration based on center voxel failed. Select the brightest pixel (b) or a custom pixel (c) to perform the b1 peak calculation: ', 's');
                    switch resp
                        case {'b','B'}
                            fprintf('You chose brightest pixel option.\n');
                            method = 'brightest';
                % --- call your brightest pixel routine here ---
                        case {'c','C'}
                            fprintf('You chose custom pixel option.\n');
                            method = 'custom';
                % --- call your custom pixel picker here --
                        
                        otherwise
                            fprintf('Invalid input. Must be b or c.\n')
                    end
                end

            case 'brightest'
                [ok, b1_peak] = rf2d_create.attempt_brightest(image_in, flip_angle, amp, phase, tp, x_axis, y_axis);
                if ok
                    done = true;
                else
                    % fallback: brightest → center → custom
                    resp = input('B1 Calibration based on brightest voxel failed. Select the center pixel (b) or a custom pixel (c) to perform the b1 peak calculation: ', 's');
                    switch resp
                        case {'b','B'}
                            fprintf('You chose brightest pixel option.\n');
                            method = 'brightest';
                % --- call your brightest pixel routine here ---
                        case {'c','C'}
                            fprintf('You chose center pixel option.\n');
                            method = 'center';                        
                        otherwise
                            fprintf('Invalid input. Must be b or c.\n')
                    end
                end

            case 'custom'
                [ok, b1_peak] = rf2d_create.attempt_custom(image_in, flip_angle, amp, phase, tp, x_axis, y_axis);
                if ok
                    done = true;
                else
                    % rule: custom failure ⇒ custom error (no more fallbacks)
                    error('B1Peak:CustomFailed', ...
                          'Custom pixel B1 peak calculation failed as the selected point has no excitation. Retry again with a different excitation shape or select a different point.');
                end
        end

        % If we just tried center or brightest and it failed twice (A↔B), drop to custom:
        if ~done && loops == 2 && ~strcmp(method,'custom')
            method = 'custom';
        end
    end
else
    b1_peak = 2.45e-6;  % fallback
end
waveform = [];
waveform(:,1) = amp_data;
waveform(:,2) = phase_data;
waveform(:,3) = zeros(length(amp_data(:,1)), 1);
all_gradients = [];
all_gradients(:,1) = g_x;
all_gradients(:,2) = g_y;
all_gradients(:,3) = zeros(length(g_x(:,1)), 1);

rf2d = struct();
rf2d.waveform = waveform;
rf2d.type = '2drf';
rf2d.gradients = all_gradients;
rf2d.b1_peak = b1_peak;
rf2d.flip_angle = flip_angle;
rf2d.tp = tp;
save('rf2d.mat', 'rf2d');   % write to disk
elapsed = toc;
fprintf('Elapsed time: %.3f seconds\n', elapsed);
%---------------- ANIMATION --------------------------------------

    end


    end
end