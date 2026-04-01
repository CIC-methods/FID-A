% rf2d_samplescript.m
% Jason Rock, Elena Osipyan and Jamie Near, 2025.
%
%
% DESCRIPTION:
% This script demonstrates how to generate a 2D radiofrequency (2DRF) pulse
% and simulate its excitation profile using a Bloch simulator. The script 
% initializes key scan parameters, constructs the RF waveform using 
% rf2d_create.makepulse, and then runs a Bloch simulation via 
% rf2d_bloch_sim.runblochsim.
%
% INPUTS (as fields of params structure):
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
% 'center' (gets the center pixel of the image), 'brightest' (gets the brightest pixel of the image), 'custom' (lets user choose custom pixel of image).
%
% OUTPUTS:
% rf2d        = Struct containing generated 2DRF pulse and gradient waveforms
% Bloch sim   = Visualization of excitation pattern from bloch_sim_2drf.runblochsim

params = struct();

% Hardware limits (must not exceed scanner specs; e.g., 80 mT/m and 200 T/m/s on many 3T systems)
params.MAX_GRAD     = 0.08;    % [T/m]  Maximum gradient amplitude (sets k-space reach per unit time)

% Slew is in T/m/ms per header; equals 200 T/m/s. Controls ramp times and achievable spiral density.
params.MAX_SLEW     = 0.2;     % [T/m/ms] Maximum gradient slew rate (limits gradient ramps; impacts N/variable-density)

% Rasterization (DAC) constraints — minimum hardware time-steps for RF and gradients
params.RF_RASTER    = 5e-6;    % [s]   RF raster time (discrete RF dwell; impacts waveform length quantization)
params.GRAD_RASTER  = 1e-5;    % [s]   Gradient raster time (discrete grad dwell; enforces timing quantization)

% Pulse type and spatial target
params.type         = '2drf';  % Type of RF pulse ('2drf' expected; other values are rejected)
params.shape        = rf2d_shapes(2); % 2D target excitation pattern (matrix from rf2d_shapes; e.g., 1=square). 

% Spatial/temporal design knobs
params.FOV          = 12;      % [cm]  Field of view for the excitation pattern (design/grid FOV)
params.tp           = 0.015;   % [s]   Total RF pulse duration (sets time–bandwidth, SAR, and B1 efficiency)
params.k_us         = 4;       % [-]   Transmit k-space undersampling factor (>=1; 1 = fully sampled; higher = faster/aliasing trade-offs)

% Flip angle and frequency offset
params.flip_angle   = 0.174533;  % [rad] Nominal small-tip flip angle (10°); used in B1 peak calibration
params.ppm          = 0;       % [ppm] Off-resonance/chemical-shift for simulation (0 = on-resonance)

% Discretization for design/simulation
params.voxels       = 128;     % [-]   Number of spatial samples per axis for design/Bloch sim grid (power-of-two convenient)

% B1 calibration policy for setting rf2d.b1_peak (affects small-tip scaling and validation)
params.calibration  = 'center';% {'center','brightest','custom'} method for B1 peak estimation (with black-pixel guard checks)


save('params.mat', 'params');
load params

rf2d_create.makepulse(params); % Creates the 2DRF pulse

load rf2d
rf2d_plotWaveform(rf2d,'all'); % Sample use of function to plot amplitude, phase and gradients


rf2d_bloch_sim.runblochsim(rf2d); % Runs the Bloch simulator with the 2DRF pulse


