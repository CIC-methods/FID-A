% rf2d_pulse_readin.m
% Jason Rock, Elena Osipyan and Jamie Near, 2025.
%
% USAGE:
% [amp, phase, g_x, g_y, g_z] = rf2d_pulse_readin.read_pulse_for_bloch(cond)
%
% DESCRIPTION:
% This class provides utility functions for reading in 2D RF pulse and gradient
% waveforms from different file formats. It supports multiple input conditions
% corresponding to pulse formats from various sources (e.g., Vespa, McGill,
% Peyshan, Pauly, custom B1_maker). The output waveforms can then be used in 
% Bloch simulations or other MR sequence modeling workflows.
%
% INPUTS:
% cond        = integer specifying which input case to use:
%               1 - Vespa format text file ("MAO_vespaFormat.txt")
%               2 - Jamie’s RF pulse (loaded from "RFpulse.mat")
%               3 - Peyshan’s pulse (RF: "MAO_vespaFormat.txt", Gradients: "2DGpulse.txt")
%               4 - McGill EPI pulse (RF: "RF_EPI_McGill.txt", Gradients: "grad_EPI_McGill.txt")
%               5 - Custom RF + gradient (RF: "2DRFpulse.txt", Gradients: "2DGpulse.txt")
%               6 - Pauly analytical pulse (generated internally, no files required)
%               7 - Peyshan’s sequence-compatible RF + gradient files
%
% FILE FORMATS:
% - RF text files: contain amplitude and phase (in degrees or radians).
% - Gradient text files: contain x, y, (optionally z) gradient values per line.
% - .mat file: should contain struct/array "RFpulse" with phase and amplitude columns.
%
% OUTPUTS:
% amp         = RF amplitude waveform
% phase       = RF phase waveform
% g_x         = Gradient waveform in X direction
% g_y         = Gradient waveform in Y direction
% g_z         = Gradient waveform in Z direction (may be zero if not provided)


classdef rf2d_pulse_readin
    methods(Static)
    function lines = readfiles(f)
      lines = {};
      tline = fgetl(f);
      while ischar(tline)
          lines{end+1} = tline;
          tline = fgetl(f);  % <- use 'f' here instead of 'fid'
      end
      fclose(f);
    end
    function [amp, phase, g_x, g_y, g_z] = read_pulse_for_bloch(cond)
    amp = [];
    phase = [];

    g_x = [];
    g_y = [];
    g_z = [];
    if(cond == 1)
        f = fopen('MAO_vespaFormat.txt', 'r');
        lines = rf2d_pulse_readin.readfiles(f);
         for i = 1:length(lines)
                tokens = str2double(strsplit(lines{i}));
                amp(end+1) = tokens(1);
                phase(end+1) = tokens(2);
         end
    end
    if(cond == 2)  % Jamie's pulse
        data = load('RFpulse.mat');  % Loads struct
        RFpulse = data.RFpulse;      % Extract the variable

        % Assuming RFpulse is an Nx2 array: column 1 = phase, column 2 = amp
        phase = RFpulse(:, 1);
        amp = RFpulse(:, 2);
    end

            
    if(cond == 3) %Peyshan's pulse
        f = fopen('MAO_vespaFormat.txt', 'r');
        lines = rf2d_pulse_readin.readfiles(f);
            for i = 1:length(lines)
                if i > 0
                    tokens = str2double(strsplit(lines{i}));
                    amp(end+1) = tokens(1);
                    phase(end+1) = tokens(2);
                end
            end
       
        f2 = fopen('2DGpulse.txt', 'r');
        lines = rf2d_pulse_readin.readfiles(f2);

        for i = 1:length(lines)
                if i > 0
                    tokens = str2double(strsplit(lines{i}));
                    g_x(end+1) = 20 * tokens(1);
                    g_y(end+1) = 20 * tokens(2);
                    g_z(end+1) = 10 * tokens(3);
                end
        end
    end
    if(cond == 4) %for Jamie McGill pulse
        f = fopen('RF_EPI_McGill.txt', 'r');
        lines = rf2d_pulse_readin.readfiles(f);
            for i = 1:length(lines)
               tokens = str2double(strsplit(strtrim(lines{i}), ','));
               if ~isempty(tokens) && ~any(isnan(tokens))
                   amp(end+1) = tokens(1);
                   phase(end+1) = tokens(2);
               end
            end

         f2 = fopen('grad_EPI_McGill.txt', 'r');
         lines = rf2d_pulse_readin.readfiles(f2);
             for i = 1:length(lines)
               tokens = str2double(strsplit(lines{i}, ','));
               if ~isempty(tokens) && ~all(isnan(tokens))
                   g_x(end+1) = 9.63 * tokens(1);
                   g_y(end+1) = 2.10 * tokens(2);
                   g_z(end+1) = 0 * tokens(3);  % or just g_z(end+1) = 0;
               end
             end
    end
    if cond == 5  % My custom RF from B1_maker
    Gx_max = 0;
    Gy_max = 0;
    
    % --- Read RF pulse file ---
    data = readmatrix('2DRFpulse.txt');
    amp = data(:, 1);  % Extract column 1
    phase = data(:, 2);  % Extract column 2
    % fid = fopen(filename, 'r');
    % if fid == -1
    %     error('Could not open file: %s', filename);
    % end
    % try
    %     fgetl(fid);  % Skip header
    %     line_count = 0;
    %     while ~feof(fid)
    %         line = fgetl(fid);
    %         if ischar(line) && ~isempty(strtrim(line))
    %             line_data = str2num(line); %#ok<ST2NM>
    %             if length(line_data) >= 2
    %                 amp(end+1) = line_data(1);
    %                 phase(end+1) = line_data(2);
    %             end
    %         end
    %     end
    % catch ME
    %     fclose(fid);
    %     rethrow(ME);
    % end
    % fclose(fid);

    % --- Read gradient pulse file ---
    gradfile = '2DGpulse.txt';
    f2 = fopen(gradfile, 'r');
    if f2 == -1
        error('Could not open file: %s', gradfile);
    end

    lines1 = {};
    tline = fgetl(f2);
    while ischar(tline)
        lines1{end+1} = tline;
        tline = fgetl(f2);
    end
    fclose(f2);

    % Process gradient file
    for i = 1:length(lines1)
        tokens = str2double(split(strtrim(lines1{i})));
        if i == 2
            Gx_max = tokens(1);
            Gy_max = tokens(2);
        elseif i > 2 && length(tokens) >= 2
            g_x(end+1) = Gx_max * tokens(1);
            g_y(end+1) = Gy_max * tokens(2);
        end
    end
    end

    if cond == 7 % Peyshan's sequence compatible pulse
                Gx_max = 0;
                Gy_max = 0;
                
                % --- Read RF pulse file ---
                filename = '2DRFpulse.txt';
                fid = fopen(filename, 'r');
                if fid == -1
                    error('Could not open file: %s', filename);
                end
                
                try
                    % Read all lines first
                    lines = rf2d_pulse_readin.readfiles(fid);
                    
                    % Process RF data (skip header if present)
                    start_idx = 1;
                    if ~isempty(lines) && isempty(str2num(lines{1}))
                        start_idx = 2; % Skip header
                    end
                    
                    for i = start_idx:length(lines)
                        line = strtrim(lines{i});
                        if ~isempty(line)
                            line_data = str2num(line); %#ok<ST2NM>
                            if length(line_data) >= 2
                                amp(end+1) = line_data(1);
                                phase(end+1) = line_data(2);
                            end
                        end
                    end
                    
                catch ME
                    if fid ~= -1
                        fclose(fid);
                    end
                    rethrow(ME);
                end
                
                % --- Read gradient pulse file ---
                gradfile = '2DGpulse.txt';
                f2 = fopen(gradfile, 'r');
                if f2 == -1
                    error('Could not open file: %s', gradfile);
                end
                
                try
                    lines1 = rf2d_pulse_readin.readfiles(f2);
                    
                    % Process gradient file
                    for i = 1:length(lines1)
                        line = strtrim(lines1{i});
                        if ~isempty(line)
                            tokens = str2double(split(line));
                            tokens = tokens(~isnan(tokens)); % Remove NaN values
                            
                            if i == 2 && length(tokens) >= 2
                                Gx_max = tokens(1);
                                Gy_max = tokens(2);
                            elseif i > 2 && length(tokens) >= 2
                                g_x(end+1) = Gx_max * tokens(1);
                                g_y(end+1) = Gy_max * tokens(2);
                            end
                        end
                    end
                    
                catch ME
                    if f2 ~= -1
                        fclose(f2);
                    end
                    rethrow(ME);
                end
                
                % Data validation and debugging
                fprintf('RF Data: %d amplitude points, %d phase points\n', length(amp), length(phase));
                fprintf('Gradient Data: %d Gx points, %d Gy points\n', length(g_x), length(g_y));
                fprintf('Amplitude range: [%.4f, %.4f]\n', min(amp), max(amp));
                fprintf('Phase range: [%.4f, %.4f]\n', min(phase), max(phase));
                
                % Convert phase from degrees to radians if needed
                if max(abs(phase)) > 10 % Likely in degrees
                    fprintf('Converting phase from degrees to radians\n');
                    phase = phase * pi / 180;
                end
      end
        

    if(cond == 6) %pauly pulse exact
        gamma = 267.522*10^6;
        T = 8; %ms
        Gmax = 6; %mT/m 
        alpha = 1; %90 degree flip 
        beta = 2;
        n = 8;
        t_i = linspace(0, T, 1000);
        for i = 1:1000
            b1_temp = alpha * gamma / T * exp(-beta^2 * (1 - t_i(i)/T)^2) * sqrt((2 * pi * n * (1 - t_i(i)/T))^2);
            amp(i) = b1_temp;
            phase(i) = 0;
        end
            
        for i = 1:length(t_i)
             gx_val = -1 * (2*pi*n * (1 - t_i(i)/T) * sin(2*pi*n*t_i(i)/T) + cos(2*pi*n*t_i(i)/T));
             gy_val =      (2*pi*n * (1 - t_i(i)/T) * cos(2*pi*n*t_i(i)/T) - sin(2*pi*n*t_i(i)/T));
             g_x(end+1) = gx_val;
             g_y(end+1) = gy_val;
        end
        g_x = g_x / max(g_y) * 6;
        g_y = g_y / max(g_y) * 6;
    end
  end
 end
end

