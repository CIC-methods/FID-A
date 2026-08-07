% rf2d_bloch_sim.m
% Jason Rock, Elena Osipyan and Jamie Near, 2025.
%
% USAGE:
% rf2d_bloch_sim.runblochsim(params)
%
% DESCRIPTION:
% This class performs Bloch simulations of two-dimensional RF (2DRF) pulses. 
% It supports both simulated pulses generated with rf2d_create.makepulse 
% and Siemens-provided RF/gradient waveform data. The simulator models 
% magnetization dynamics under off-resonance conditions and can visualize 
% excitation maps, magnetization components of 
% the Bloch simulation.
%
% INPUTS:
% params.RF_RASTER   = RF raster time [s]
% params.GRAD_RASTER = Gradient raster time [s]
% params.tp          = Total pulse duration [s]
% params.flip_angle  = Flip angle [degrees] (default = 10°)
% cond               = Integer flag to select pulse type (e.g., cond=7 for 
%                      rf2d_create.m-generated pulses)
% ppm                = Spin system frequency offset definition. Desired off-resonance offset [Hz]. Based off the ppm
% value in params.mat
% file_path          = Path to RF pulse text file (default = '2DRFpulse.txt')
% b1_peak            = Peak B1 amplitude [Tesla] (default ≈ 8e-6 for 2DRF)
%
% OUTPUTS:
% result      = Complex transverse magnetization (Mxy) across spatial grid
% M_z_range   = Final z-component magnetization values
% M_x_range   = Final x-component magnetization values
% M_y_range   = Final y-component magnetization values
% (All outputs can be visualized as images or trajectories of magnetization)
%
% INTERNAL METHODS:
% rot_y()     - Rotation matrix about the y-axis
% rot_z()     - Rotation matrix about the z-axis
% Bloch_Sim() - Core Bloch simulation routine, computes Mx, My, Mz
% runblochsim()- Main execution function, loads pulse/gradients and runs sim
%
% VISUALIZATIONS:
% - 3D plot of magnetization components (Mx, My, Mz)
% - Transverse excitation map (|Mxy|)
% - RF amplitude/phase vs time
% - Gradient waveforms vs time


classdef rf2d_bloch_sim
    methods(Static)
      
% Rotation Functions
function[res] = rot_y(M_temp,angle)
    angle = angle(1);
    R = [cos(angle),0,-sin(angle);
         0, 1, 0;
         sin(angle), 0, cos(angle)];
    res = R * M_temp(:);
    res = res';
end 

function[res] = rot_z(M_temp, angle)
    angle = angle(1);
    R = [cos(angle), -sin(angle), 0;
         sin(angle), cos(angle), 0;
         0, 0, 1];
    res = R * M_temp(:);  % Ensure proper matrix multiplication
    res = res';
end


function[result, M_z_range, M_y_range, M_x_range, Mxy] = Bloch_Sim(b0_map_flag, gamma, x_pos, y_pos, grad_x, grad_y, offset, b1, phase, amp, b1_peak)
    RF_RASTER = 5*10^-6;
    tp = round(length(amp)*RF_RASTER, 5); 
    b0_map = fliplr(readNPY('b0map.npy'));
    b0_map = imresize(b0_map, [length(y_pos), length(x_pos)]);
    figure;
    imagesc(b0_map);           % Scales and displays matrix as image
    colormap jet;              % Apply the 'jet' colormap
    colorbar;                  % Add colorbar
    title('loaded b0 map');    % Set title
    axis image;                % Optional: makes pixel aspect ratio 1:1
    hold off;
    

    fprintf("**Bloch Sim Starting with this offset:" + offset); 
    for n = 1:length(x_pos)
       offset_x = grad_x .* x_pos(n) / 1000;
        for ny = 1:length(y_pos)
            offset_y =  grad_y .* y_pos(ny) / 1000; % + 300/gamma #for off res
            offset_b0 = 2*pi*b0_map(ny,n)/gamma;
            M_new = [0 0 1];
           for q = 1:length(b1)-1
             phi = phase(q);

             if(b0_map_flag)
                B1_eff = sqrt(b1(q).^2 + (offset_x(q) + offset_y(q) + offset/gamma + offset_b0).^2);
                alpha = atan2(b1(q), (offset_x(q) + offset_y(q) + offset_b0 + offset/gamma));
          
             else
                B1_eff = sqrt(b1(q).^2 + (offset_x(q) + offset_y(q) + offset/gamma).^2);
                alpha = atan2(b1(q), (offset_x(q) + offset_y(q) + offset/gamma));
   
             end
            dt = tp / length(amp);
            theta = gamma * B1_eff * dt;
            M_new = rf2d_bloch_sim.rot_z(M_new, phi);
            M_new = rf2d_bloch_sim.rot_y(M_new, alpha);
            M_new = rf2d_bloch_sim.rot_z(M_new, theta);
            M_new = rf2d_bloch_sim.rot_y(M_new, -alpha);
            M_new = rf2d_bloch_sim.rot_z(M_new, -phi);
    
           end
           M_z_range(ny,n) = M_new(3);
           M_x_range(ny,n) = M_new(1);
           M_y_range(ny,n) = M_new(2);
        end  
    end
   % Compute transverse magnetization magnitude

N = length(M_x_range);
idx = 1:N;

figure;

% X component along x-axis
pMx = plot3(idx, M_x_range, zeros(1, N), 'r', 'DisplayName', 'M_x'); 
hold on;

% Y component along y-axis
pMy = plot3(idx, zeros(1, N), M_y_range, 'g', 'DisplayName', 'M_y');

% Z component rotated to lie along z-axis
pMz = plot3(idx, M_z_range, zeros(1, N), 'b', 'DisplayName', 'M_z');

xlabel('M_z');
ylabel('M_y');
zlabel('M_x');
title('Magnetization Components vs Sample Index');

legend([pMx(1) pMy(1) pMz(1)], 'Location', 'northeastoutside');  % compact legend
grid on;
view(45, 30);  % adjust viewing angle

Mxy = M_x_range + 1i * M_y_range;
% --- figure: Transverse excitation map ---
figure; 

imagesc(100 * [min(x_pos), max(x_pos)], ...
         100 * [min(y_pos), max(y_pos)], ...
         (abs(Mxy)));
colormap('gray');
axis image;
set(gca, 'FontSize', 20);
colorbar;
title('Transverse Excitation', 'Interpreter', 'latex');
xlabel('X position (cm)');
ylabel('Y position (cm)');
result = Mxy; 
end

function[] = runblochsim(params)
   tic;
   load params
   load rf2d
   if ~strcmp(rf2d.type, '2drf')
    ME = MException('ForceStop:ExecutionTerminated', ...
        'You are not creating a 2DRF Pulse! Program stopped.');
    throw(ME);
   end
   RF_RASTER = params.RF_RASTER;
   GRAD_RASTER = params.GRAD_RASTER;
   b1_peak = rf2d.b1_peak; %flip_angle/gamma/tp/integral_factor #Tesla #8 for regular pulses
   chris_pulse_test = false;
   cond = 7;
   [amp, phase, g_x, g_y, g_z] = rf2d_pulse_readin.read_pulse_for_bloch(cond);
   ppm = params.ppm; 
   offset = 123.25*ppm*2*pi;
   gamma = 267.522*10^6; 

T = params.tp;  % seconds - total pulse duration.
t = linspace(0, T, 6000);  % rf pulse legnth: 10ms, # of samples
% 6000 = lots of samples

off_resonance_testing = false;
b0_map_flag = false;

% Gradients resampled to be same length as the phase and amplitude waveforms.
g_x = resample(g_x, length(amp), length(g_x));
g_y = resample(g_y, length(amp), length(g_y));

b1_peak_approx = 8e-6;
b1 = ((amp) / max(amp))* b1_peak_approx;
b1 = reshape(b1.', 1, []);  % transpose first, then reshape into a row vector


tp = round(length(amp)*RF_RASTER, 5);
t_axis = linspace(0, tp * 1000, length(amp));  % ms
% % Plot RF amplitude
% figure;
% plot(t_axis, amp, 'LineWidth', 1);
% title("RF Pulse Amplitude");
% xlabel("Time (ms)");
% ylabel("Amplitude");
% grid on;

% % Plot RF phase
% figure;
% plot(t_axis, phase, 'LineWidth', 1);
% title("RF Pulse Phase");
% xlabel("Time (ms)");
% ylabel("Phase (radians or a.u)");
% grid on;


tp = round(length(amp)*RF_RASTER, 5); 
dt = tp / length(amp);


if(chris_pulse_test)
    data = load('sim_noB0.npz');
    data = load('jason_square_B0_23.npz');
    lst = data.files;
    new_amp = [];
    for i = 1:length(lst)
        item = lst{i};
        new_amp{end+1} = data.(item);
    end
    amp = new_amp{1};
    phase = new_amp{2};
    data = [amp, phase];  % Concatenate as two columns
    
    lenRF = length(amp);
    file_path = '2DRFpulse.txt';
    
    %Write the data to a text file with tab separation
    
     fid = fopen(file_path, 'w');
     fprintf(fid, '%d\n', lenRF);
     fclose(fid);
     % Append the amplitude and phase data
     writematrix(file_path, data, '-append', 'delimiter', '\t', 'precision', '%.6f');
    
    fprintf('Complex B1 has been written to {file_path}');

end


RF_lobe_correction = zeros(size(amp));
signed_RF = zeros(size(amp));

realAmp = 0;
imagAmp = 0;

if cond == 1 || cond == 2
    for i = 1:length(phase)
        if(phase(i) == 0.0)
            RF_lobe_correction(i) = 1;
        end
        if(phase(i) == 180.0)
            RF_lobe_correction(i) = -1;
        end
        signed_RF(i) = amp(i)*RF_lobe_correction(i);
        %signed_RF has the arbitrary unit scaling
    integral_factor = sum(signed_RF)/(length(signed_RF)*max(signed_RF));
    end
end
if cond == 3 || cond == 4
    %integral_factor = 82
    for i = 1:length(amp)
        realAmp = realAmp + amp(i) * cos(phase(i));
        imagAmp = imagAmp + amp(i) * sin(phase(i));
        
    end
end

integral_factor = sqrt(realAmp^2 + imagAmp^2);


if(cond == 1 || cond == 2)
  grad_x = 0.5595 * ones(size(b1));  % mT/m
  grad_y = 0.5595 * ones(size(b1));
end 
if(cond == 3 || cond == 7) %Siemens sequence requirements
    grad_x = resample(g_x, length(b1), length(g_x));
    grad_y = resample(g_y, length(b1), length(g_y));
end
    
    
   
if(cond == 4 || cond == 5 || cond == 6)
   t_grad = linspace(0, 1, length(g_x));
   t_b1 = linspace(0, 1, length(b1));
   % Interpolate gradients to match b1 length
   grad_x = g_x;  % Or 'spline'
   grad_y = interp1(t_grad, g_y, t_b1, 'spline');
end 

flip_angle = params.flip_angle;

x_fov = params.FOV/100; %meters
y_fov = params.FOV/100; 
linspace(x_fov/2, -x_fov/2, params.voxels); 
M0 = [0,0,1];

b1_unsgn = b1; %(np.array(amp) / max(amp)) * b1_peak
a = 0.40;
beta = pi / 2;

M_temp = [0 0 1];
amp = amp(:);
phase = phase(:);

tp = RF_RASTER * length(amp);  % update time duration


x_pos = linspace(x_fov/2, -x_fov/2, params.voxels); 
y_pos = linspace(y_fov/2, -y_fov/2, params.voxels); 

images = zeros(length(amp), length(x_pos), length(y_pos), 'single');
images_pha = zeros(length(amp), length(x_pos), length(y_pos), 'single');

blochs_simside = {};
if(off_resonance_testing)
    ORs = [0, 66.5, 194.4, 217, 263] .* (2 * pi);
    for i = 1:length(ORs)
        bloch_simside = rf2d_bloch_sim.Bloch_Sim(b0_map_flag, gamma, x_pos, y_pos, grad_x, grad_y, offset, b1, phase, amp, b1_peak);
    end

    instance = 7;
    % instance 2 = circle 
    image = shape(instance);
    image = fliplr(image); 
    mse_onres, roi_integrals, frac_outside = off_res_metrics(image, blochs_simside, '.', FA = 1);
    
    for i = 1:length(ORs)
        figure;
        offset_Hz = 0;
        titleStr = sprintf('Off-resonance Transverse Excitation M_{xy}, nOffset = %g', offset);
        Mxy = M_x_range + 1i*M_y_range;
        % Plot the masked image
        imagesc(100 * [min(x_pos), max(x_pos)], ...
        100 * [min(y_pos), max(y_pos)], ...
        masked_image);
        axis image;
        colormap gray;
        
        c = colorbar;
        c.Label.String = 'Fraction of \bfM in M_{xy}';
        title(titleStr, 'Interpreter','tex')
        xlabel("x position (cm)")
        ylabel("y position (cm)")
    end 
    
else
    blochs_simside = rf2d_bloch_sim.Bloch_Sim(b0_map_flag, gamma, x_pos, y_pos, grad_x, grad_y, offset, b1, phase, amp, b1_peak);
end

elapsed = toc;
fprintf('Elapsed time of Bloch Sim: %.3f seconds\n', elapsed);
end
end
end