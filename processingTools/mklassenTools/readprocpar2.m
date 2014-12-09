function [par,fpath] = readprocpar2(fpath, ext)

% READPROCPAR reads in parameter used in FID acquisition
%
%   PAR = READPROCPAR(FPATH) reads the parameters in the the
%   procpar file in the FID directory FPATH
%
% see also WRITEPROCPAR, READFID, WRITEFID


% Written by L Martyn Klassen
% Copyright 2009 Robarts Research Institute

if nargin < 1
   error('CFMM:readprocpar','One input argument required.');
end

if nargin < 2
   ext = [];
end

% Check for validity of file
fp = fopen([fpath '/procpar' ext], 'r', 'ieee-be');
if fp == -1
   fp = fopen([fpath '.fid/procpar' ext], 'r', 'ieee-be');
   if fp == -1
      fp = fopen([fpath '.par/procpar' ext], 'r', 'ieee-be');
      if (fp == -1)
         fp = fopen([fpath ext], 'r', 'ieee-be');
         if fp == -1
            error('CFMM:readprocpar2','Unable to open %s for reading.', fpath);
         end            
      else
         fpath = [fpath '.par'];
      end
   else
      fpath = [fpath '.fid'];
   end
end

% VnmrJ support special characters in variable names but MATLAB does not
% support them in field names. A special field must be created
specialfield = 'zzzzz_not_a_valid_MATLAB_field_name';

% Read in all the parameters
while (true)
   [field, count] = fscanf(fp, '%s', 1);
   
   if (count == 0)
      break
   end
   
   param = getparameter(field, fp);
   
   if (~isempty(regexp(field, '\W', 'once')) || regexp(field, '[a-zA-Z]','once') ~= 1)
      % Variable name is not a valid field name in matlab
      if isfield(par, specialfield)
         par.(specialfield)(end+1) = param;
      else
         par.(specialfield) = param;
      end
   else
      par.(field) = param;
   end
end

% Close file
fclose(fp);

% Order the fields alphabetically
par = orderfields(par);

% Do some special parsing for array and rcvr
if isfield(par, 'rcvrs')
   % rcvr is not allowed to be arrayed
   par.rcvrs.number = sum(par.rcvrs.value{1} == 'y');
end

if isfield(par, 'array')
   % array is not allowed to be arrayed
   buffer = par.array.value{1};
   % Parse the data string
   index1 = 1;
   index2 = 1;
   incr = 0;
   par.array.loops{1}{1} = [];
   for o = 1:length(buffer)
      switch buffer(o)
      case '('
         incr = incr + 1;
         if incr > 1
            error('CFMM:readprocpar2','invalid coupling in array paramater.');
         end
      case ')'
         incr = incr - 1;
         if incr == 0
            index2 = 1;
         elseif incr < 0
            error('CFMM:readprocpar2','invalid coupling in array paramater.');
         end
      case ','
         if incr > 0
            index2 = index2 + 1;
         else
            index1 = index1 + 1;
         end
         par.array.loops{index1}{index2} = [];
      otherwise
         par.array.loops{index1}{index2} = [par.array.loops{index1}{index2} buffer(o)];
      end
   end
end
return

function [param] = getparameter(field, fp)

   param.name        = field;
   A = fscanf(fp, '%d %d %e %e %e %d %d %d %d %d', 10);
   param.subtype     = A(1);
   param.basictype   = A(2);
   param.maxvalue    = A(3);
   param.minvalue    = A(4);
   param.stepsize    = A(5);
   param.Ggroup      = A(6);
   param.Dgroup      = A(7);
   param.protection  = A(8);
   param.active      = A(9);
   param.intptr      = A(10);
  
   % Skip the rest of the line
   fgets(fp);
   loop = fscanf(fp,'%d',1);
   switch param.basictype
      case 1
         param.value = sscanf(fgets(fp), '%f').';
         % skip the last line
         fgets(fp);
      case {0, 2}
         fscanf(fp,'%[^"]');
         for i=1:loop
            string = fscanf(fp, '"%[^"]"');
            quotes = string;
            while (~isempty(quotes) && (string(end) == '\'))
               string(end) = '"';
               quotes = fscanf(fp, '%[^"]"');
               string = [string quotes]; %#ok<AGROW>
            end
            param.value{i} = string;
            fgets(fp);
         end
         % Read in the enumeration values
         loop = fscanf(fp, '%d', 1);
         for i=1:loop
            fscanf(fp,'%[^"]');
            string = fscanf(fp, '"%[^"]"');
            if isempty(string)
               fscanf(fp,'"');
            end
            quotes = fscanf(fp, '%[^ ]');
            param.enumeral{i} = [string quotes];
         end
         fgets(fp);
   end
 return  