function [rf, info] = readrfvnmr(filename,verbosity)

% RF = READRFVNMR(FILENAME) reads a Varian VNMR format RF pulse
% description in FILENAME and generates an RF matrix with 3 columns
% specifying angle, magnitude, and duration.
%
% FILENAME can be in working directory, ~/vnmrsys/shapelib, or 
% $vnmrsystem/shapelib and can be specified with or with trailing .RF

% Author: L. Martyn Klassen
% Copyright 2003 Robarts Research Institute
% This program is copyrighted worldwide by the Robarts Research
% Institute.  Any distribution, copying or redistribution is expressly
% forbidden.
% Created: April 17, 2003 14:22
% Modified: N/A

if nargin < 2
   verbosity = 0;
   if nargin < 1
      error('READRFVNMR requires RF file name.');
   end
end

% info defaults to an empty structure
info = struct;

% Get the vnmrsystem path
[status, vnmrpath] = unix('echo $vnmrsystem');
if status == 0
   % Strip end-of-line character
   vnmrpath = vnmrpath(1:end-1);
else
   % Make a guess - good for RRI
   vnmrpath = '/vnmr';
end

fid = findfile(filename, vnmrpath);
if fid < 0
   % Check if user forgot .RF or .rf ending
   orgfilename = filename;
   filename = [orgfilename '.RF'];
   fid = findfile(filename, vnmrpath);
   if fid < 0
      filename = [orgfilename '.rf'];
      fid = findfile(filename, vnmrpath);
      if fid < 0
         error('READRFVNMR unable to find %s', orgfilename);
      end
   end
end

if verbosity > 0
   % Report what file was actually opened for reading
   disp(sprintf('READRFVNMR reading file %s', fopen(fid)));
end

% Get the first line of the file
fileline = fgets(fid);
linenumber = 1;
rf = zeros(3,0);
while ischar(fileline)
   % Strip out any comment indicated by #
   idxcomment = strfind(fileline,'#');
   if ~isempty(idxcomment)
      % Look for information normally stored in header
      % It must be commented so we don't need to do this every line
      idx = strfind(fileline, 'VERSION');
      if ~isempty(idx)
         info.version = sscanf(fileline(idx+7:end), '%s', 1);
      end
      idx = strfind(fileline, 'TYPE');
      if ~isempty(idx)
         info.type = sscanf(fileline(idx+4:end), '%s', 1);
      end
      idx = strfind(fileline, 'MODULATION');
      if ~isempty(idx)
         info.modulation = sscanf(fileline(idx+10:end), '%s', 1);
      end
      idx = strfind(fileline, 'EXCITEWIDTH');
      if ~isempty(idx)
         info.excitewidth = sscanf(fileline(idx+11:end), '%f', 1);
      end
      idx = strfind(fileline, 'INVERTWIDTH');
      if ~isempty(idx)
         info.invertwidth = sscanf(fileline(idx+11:end), '%f', 1);
      end
      idx = strfind(fileline, 'INTEGRAL');
      if ~isempty(idx)
         info.integral = sscanf(fileline(idx+8:end), '%f', 1);
      end
      idx = strfind(fileline, 'RF_FRACTION');
      if ~isempty(idx)
         info.aref = sscanf(fileline(idx+11:end), '%f', 1);
      end
      idx = strfind(fileline, 'AREF');
      if ~isempty(idx)
         info.aref = sscanf(fileline(idx+4:end), '%f', 1);
      end
      idx = strfind(fileline, 'b1excite');
      if ~isempty(idx)
         info.b1excite = sscanf(fileline(idx+8:end), '%f', 1);
      end
      idx = strfind(fileline, 'b1invert');
      if ~isempty(idx)
         info.b1invert = sscanf(fileline(idx+8:end), '%f', 1);
      end
      
      % Strip out everything after the comment indicator
      fileline = fileline(1:idxcomment(1));
   end
   
   % Kill the last character: end-of-line or #
   fileline(end) = [];
   
   % If the line is empty skip it
   if ~isempty(fileline)
      % if it blank then skip
      if any(fileline ~= ' ')
         [A,count, errmsg, nextindex] = sscanf(fileline, '%f', inf);
         
         % If read failed, output the error
         if ~isempty(errmsg)
            fclose(fid)
            error('READRFVNMR failed with read error: %s', errmsg);
         end
         
         % Each line is to contain 3 and only 3 numbers
         if (count > 0) & (count ~= 3)
            fclose(fid);
            error('READRFVNMR malformed line number %d', linenumber);
         end
         
         % Store the read values into rf array
         rf(:,linenumber) = A;
         linenumber = linenumber + 1;
      end
   end
   % Get the next line
   fileline = fgets(fid);
end
fclose(fid);

% Add a default aref value if not found
if ~isfield(info, 'aref')
   info.aref = 0.5;
end

% Switch from row to column orientation
rf = rf';
return


function fid = findfile(filename, vnmrpath)
% Search local, user vnmrsys, and global vnmrsys for RF file
fid = fopen(filename, 'r');
if fid < 0
   orgfilename = filename;
   filename = ['~/vnmrsys/shapelib/' orgfilename];
   fid = fopen(filename, 'r');
   if fid < 0
      filename = [vnmrpath '/shapelib/' orgfilename];
      fid = fopen(filename, 'r');
   end
end
return
