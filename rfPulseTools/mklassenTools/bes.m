function [mvector,scan] = bes(RF, pw, scantype, var1, var2, var3, var4, phase, M0, T1, T2)

% Bloch Equation Simulator
%
% Required Arguments
%  RF        4xN definition of radio frequency pulse (phase, amplitude,
%               duration, gradient)
%  pw        duration of pulse in milliseconds
%  scantype  scan type: 'frequency', 'B1max', 'temporal'
%  B1max scan
%     var1     frequency / position (if gradient supplied in RF)
%     var2     upper bound B1
%     var3     lower bound B1
%     var4     steps
%  Frequency scan
%     var1     B1max
%     var2     upper bound frequency / position (if gradient supplied)
%     var3     lower bound frequency / position (if gradient supplied)
%     var4     steps
%  Temporal scan
%     var1     B1max
%     var2     frequency / position (if gradient supplied)
%     var3     ignored
%     var4     ignored
% Optional Arguments
%  phase    phase at which the RF pulse is applied (default 0)
%  M0       three component vector specifying the initial magnetization 
%              a 3xSTEPS matrix is used to define a different MO for each
%              steps. Default is [0; 0; 1]
%  T1       T1 decay time (default INF)
%  T2       T2 decay time (default INF)
%
% Output
%  mvector  3xSTEPS matrix specifying the magnetization at each point
%  scan     vector of frequency, B1max, or time steps modelled
%
% Units: B1 (kHz), frequencies (kHz), Phase (degrees), Time (ms),
%        Gradients (G/cm), Position (cm), gyro = 4.25763888 kHz/G

% Author: L. Martyn Klassen
% Copyright 2003 Robarts Research Institute
% This program is copyrighted worldwide by the Robarts Research
% Institute.  Any distribution, copying or redistribution is expressly
% forbidden.
% 03/04/24 Added support for gradient changes during RF

if nargin < 5
   error('BES requires five input arguments.');
end

% Initialize constants
scantype = lower(scantype(1));
grad_flag = 0;

if nargin < 11
   T2 = Inf;
   if nargin < 10
      T1 = Inf;
      if nargin < 9
         M0 = [0; 0; 1];
         if nargin < 8
            phase = 0;
            if (nargin < 7) && (scantype ~= 't')
               error('BES requires seven input arguments for specified scan type.');
            end
         end
      end
   end
end

% If only one column is defined in the RF pulse it is an amplitude
% modulated pulse, and the phase, duration columns have to be added
if ndims(RF) > 2
   error('BES: RF pulse cannot have more than 2 dimensions.');
end

sizRF = size(RF);
if all(sizRF > 4)
   error('BES: RF pulse has too many columns (phase, amplitude, duration, grad)')
elseif sizRF(2) > 4
   RF = RF';
end

switch size(RF, 2)
   case 1
      % Single column is RF magnitude so you have to add phase, duration,
      % and gradients
      RF(:,2) = RF;
      RF(:,1) = phase;
      RF(:,3) = 1;
      RF(:,4) = 0;
   case 2
      % Two columns are phase and magnitude, add duration.
      % Add the desired phase to the waveform:
      RF(:,1) = RF(:,1) + phase;
      RF(:,3) = 1;
   case 3
      % Three columns are phase, magnitude and duration.
      % Add the desired phase to the waveform:
      RF(:,1) = RF(:,1) + phase;
   case 4
      % Convert from G/cm to kHz/cm
      RF(:,4) = RF(:,4) .* 4.25763888;
      grad_flag = 1;
      % Add the desired phase to the waveform:
      RF(:,1) = RF(:,1) + phase;
end
      
% Convert phase to radians
RF(:,1) = RF(:,1) .* (pi / 180);

% Find out how many cyles are required to step through the RF pulse
cycles = size(RF, 1);

% Scale the RF amplitude to be between -1 and 1
RF(:,2) = RF(:,2)./max(abs(RF(:,2)));

% Kill any points with zero duration
RF = RF(RF(:,3) ~= 0, :);

% The smallest duration allowed is 1.0
minduration = min(RF(:,3));
maxduration = max(RF(:,3));
if (minduration < 1)
   error('BES: Pulse step duration should be zero or positive integer.');
end

% Adjust durations to be integers
RF(:,3) = round(RF(:,3));

% Pulse width has to be a scalar
pw = pw(1);

% Define the smallest time interval used in the pulse
dt = pw./sum(RF(:,3));

% A time scan uses the pulse squence to define the progression
% through time and not steps
if scantype == 't'
   steps = 1;
else
   steps = var4;
end

% Define M0
if isempty(M0)
   % Default is alignment along z axis
   M0 = [0; 0; 1];
   if scantype ~= 't'
      mvector = M0(:,ones(steps, 1));
   end
elseif prod(size(M0)) == 3
   % One vector must be repeated to size of the scan
   % Check to make sure that M0 is not too large
   if sum(M0.^2) > 1.0
      error('BES: Initial magnetization is greater than 1.0')
   else
      M0 = M0(:);
      if scantype ~= 't'
         mvector = M0(:,ones(steps, 1));
      end
   end 
elseif size(M0, 2) ~= steps
   if scantype ~= 't'
      % Otherwise M0 must be defined at each point in the scan
      error('BES: Initial M0 vector does not agree with number of steps (3xsteps).')
   else
      % Time requires it to 3x1 or 1x3 - previous if statement
      error('BES: Initial MO for time scan must be 3x1 vector')
   end
else
   if scantype == 't'
      error('BES: Initial MO for time scan must be 3x1 vector')
   end
   if any(sum(M0.^2, 1) > 1.000001)
      error('BES: Initial magnetization is greater than 1.0')
   end
   mvector = M0;
end

% Make sure lower limit is less than upper limit
% switch if they are the wrong order
if var2 < var3
   temp = var2;
   var2 = var3;
   var3 = temp;
   clear temp;
end

if (T2 > T1)
   error('BES: First Law Violation, T2 must be less than T1.');
end

% Calculate the matrices for the T1, T2 decay
if T1 == 0
   if T2 == 0
      A = [0 0 0; 0 0 0; 0 0 1];
   else
      A = [exp(-dt/T2) 0 0; 0 exp(-dt/T2) 0; 0 0 1];
   end
   B = [0; 0; 0];
else
   if T2 == 0
      A = [0 0 0; 0 0 0; 0 0 1-exp(-dt/T1)];
   else
      A = [exp(-dt/T2) 0 0; 0 exp(-dt/T2) 0; 0 0 exp(-dt/T1)];
   end
   B = [0; 0; 1-exp(-dt/T1)];
end

% Setup for the three different types of simulations
switch scantype
case 'f'   % Scan through frequency
    
   % Fixed value is B1max which has to be scaled for the RF pulse description
   % Also convert it to rad/millisecond
   B1max = var1.*(2*pi);
 
   % Limits are in frequency
   interval = (var2 - var3)/(steps-1);
   scan = [var3:interval:var2];

   % Rotate the matrix to align 1st RF with X
   cosa = cos(-RF(1,1));
   sina = sin(-RF(1,1));
   mvector = [cosa sina 0; -sina cosa 0; 0 0 1]*mvector;
   
   % Small tip angle vector
   theta = RF(:,2)*B1max.*dt;
   cosa = cos(theta);
   sina = sin(theta);
   
   % Free precession angle (rad)
   freq = scan .* (2*pi*dt);
   phi = freq; 
   cos2a = cos(phi);
   sin2a = sin(phi);
   
   % Alignment angles
   dRF = -diff(RF(:,1));
   dRF(end+1) = RF(end, 1);

   % Loop through the RF pulse
   for t = 1:cycles
      % Compute the rotation matrix
      Rx=[1 0 0; 0 cosa(t) sina(t); 0 -sina(t) cosa(t)];
      
      % Added in the decay matrix
      Rx = A * Rx;
      
      if grad_flag
         phi = freq .* RF(t,4);
         cos2a = cos(phi);
         sin2a = sin(phi);
      end
      
      % For durations longer than one, repeat the small angle 
      % approximation the prerquisite number of times
      for m = 1:RF(t,3)
         % Using small tip angle approximation to calculate the
         % new magnetization vector
    
         % Rotate the matrix by theta about X and let decay
         mvector = Rx*mvector;
         
         % Add the constant decay term
         mvector(3,:) = mvector(3,:) + B(3);
         
         % Free Precession angle with RF alignment on last increment
         if (m == RF(t,3)) & (dRF(t) ~= 0)
            phi2 = phi + dRF(t);
            cos3a = cos(phi2);
            sin3a = sin(phi2);
            
            temp = cos3a.*mvector(1,:) + sin3a.*mvector(2,:);
            mvector(2,:) = -sin3a.*mvector(1,:) + cos3a.*mvector(2,:);
            mvector(1,:) = temp;
         else
            % Rotate the vector by phi about Z
            temp = cos2a.*mvector(1,:) + sin2a.*mvector(2,:);
            mvector(2,:) = -sin2a.*mvector(1,:) + cos2a.*mvector(2,:);
            mvector(1,:) = temp;       
         end
      end
   end
   
case 'b'    % Scan through B1max
    
   % Fixed value is Frequency (rad/ms)
   freq = var1*2*pi;
   
   % Limits are in B1max which has to be scaled for the RF pulse description
   % and converted to rad/ms.
   interval = (var2 - var3)/(steps-1);
   B1max = [var3:interval:var2].*(2*pi);

   % Rotate the matrix to align 1st RF with X
   cosa = cos(-RF(1,1));
   sina = sin(-RF(1,1));
   mvector = [cosa sina 0; -sina cosa 0; 0 0 1]*mvector;

   % Free Precession rotation angle (rad)
   freq = freq*dt;
   phi = freq;
   % Compute the free precession matrix
   cosa = cos(phi);
   sina = sin(phi);         
   Rfree =[cosa sina 0; -sina cosa 0; 0 0 1];
   
   % Add the decay term as well
   Rfree = A * Rfree;
   
   % Calculate RF alignment angles
   dRF = -diff(RF(:,1));
   % Last one must realign to zero
   dRF(end+1) = RF(end,1);
   cos2a = cos(phi + dRF);
   sin2a = sin(phi + dRF);
    
   % Loop through the RF pulse
   for t = 1:cycles
      % Compute theta for given RF power
      theta = B1max*RF(t,2)*dt;
      cosa = cos(theta);
      sina = sin(theta);
   
      if grad_flag
         phi = freq * RF(t,4);
      end
      
      % For durations longer than one, repeat the small angle 
      % approximation the prerequisite number of times
      for m = 1:RF(t,3)
         % Using small tip angle approximation to calculate the
         % new magnetization vector
         
         % Rotate the matrix by theta about X
         temp = cosa.*mvector(2,:) + sina.*mvector(3,:);
         mvector(3,:) = -sina.*mvector(2,:) + cosa.*mvector(3,:);
         mvector(2,:) = temp;
         
         % On last time through apply free precession and rotate the matrix
         % to undo the previous alignment with the RF and also align the next
         % RF with X, otherwise just apply free precession
         if grad_flag
            % Compute the free precession matrix and apply
            cosa = cos(phi + dRF(t));
            sina = sin(phi + dRF(t));   
            mvector = (A * [cosa sina 0; -sina cosa 0; 0 0 1])*mvector;
         elseif (m == RF(t,3)) & (dRF(t) ~= 0)
            mvector=(A * [cos2a(t) sin2a(t) 0; -sin2a(t) cos2a(t) 0; 0 0 1])*mvector;
         else
            mvector = Rfree*mvector;
         end
         mvector(3,:) = mvector(3,:) + B(3);
      end
   end
  
   % Convert scan from rad/ms to kHz
   scan = B1max./(2*pi);
   
case 't'
   % Scan through time
   
   % Fixed value are B1max and freq
   B1max = var1.*(2*pi);
   freq = var2.*2*pi;
   
   mvector = zeros(3, sum(RF(:,3))+1);
   count = 1;
   
   % Store the starting magnetization
   mvector(:,count) = M0;
   
   % Free Precession rotation
   freq = freq*dt;
   phi = freq;
         
   % Loop through the RF pulse
   for t = 1:cycles
      if grad_flag
         phi = freq .* RF(t,4);
      end
      
      % Rotate the matrix to align RF with X
      cosa = cos(-RF(t,1));
      sina = sin(-RF(t,1));
      Rz1 = [cosa sina 0; -sina cosa 0; 0 0 1];
      
      % Compute the RF rotation matrix
      theta = B1max*RF(t,2)*dt;
      cosa = cos(theta);
      sina = sin(theta);
      Rx=[1 0 0; 0 cosa sina; 0 -sina cosa];

      % Also rotate the matrix to undo the previous alignment with the RF
      % and free precession
      cosa = cos(phi + RF(t,1));
      sina = sin(phi + RF(t,1));         
      Rz2=[cosa sina 0; -sina cosa 0; 0 0 1];
      
      % Add in the decay terms
      Rz2 = A * Rz2;
      
      % For durations longer than one, repeat the small angle 
      % approximation the prerquisite number of times
      for m = 1:RF(t,3)
         % Using small tip angle approximation to calculate the
         % new magnetization vector

         % RF alignment
         M0 = Rz1*M0;
         
         % Rotate the matrix by theta about X	
         M0 = Rx*M0;
         
         % Free precesion, realignment and decay
         M0 = Rz2*M0 + B;
        
         % Store the magnetization vector
         count = count + 1;
         mvector(:,count) = M0;
      end
   end
   
   scan = [0:sum(RF(:,3))].*dt;
   
otherwise
   error('BES: Scan type is frequency, B1max, or time');
end

return
