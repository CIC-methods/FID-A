function writeprocpar2(fpath, par)

% WRITEPROCPAR writes in parameter used in FID acquisition
%
%   WRITEPROCPAR(FPATH,PAR) writes the parameter structure PAR, to the 
%   procpar in the FID directory FPATH.
%
% see also READPROCPAR, WRITEFID, READFID 


% Written by L Martyn Klassen
% Copyright 2009 Robarts Research Institute

if nargin < 2
   error('cfmm:fid:name','WRITEPROCPAR: Invalid arguments, requires path and parameter set.');
end

% Check for validity of file
fp = fopen([fpath '/procpar'], 'w+', 'ieee-be');
if fp == -1
   fp = fopen([fpath '.fid/procpar'], 'w+', 'ieee-be');
   if fp == -1
      error('cfmm:fid:name',['WRITEPROCPAR unable to open ' fpath '/procpar for writing.']);
   end
end

names = fieldnames(par);

for i = 1:length(names)
   field = par.(names{i});
   
   % Special variable that cannot be field names are stored in a structure
   % array so each element of the field structure must be processed
   % separately
   for j = 1:numel(field)
      writeparam(field(j),fp);
   end
   
end  % End of while loop

% Close file
fclose(fp);

return

function writeparam(field, fp)

if (~isfield(field, 'subtype')    && ...
    ~isfield(field, 'basictype')  && ...
    ~isfield(field, 'maxvalue')   && ...
    ~isfield(field, 'minvalue')   && ...
    ~isfield(field, 'stepsize')   && ...
    ~isfield(field, 'Ggroup')     && ...
    ~isfield(field, 'Dgroup')     && ...
    ~isfield(field, 'protection') && ...
    ~isfield(field, 'active')     && ...
    ~isfield(field, 'intptr')     && ...
    ~isfield(field, 'value'))
   return
end

fprintf(fp, '%s %d %d', field.name, field.subtype, field.basictype);

writeRounded(fp, field.maxvalue);
writeRounded(fp, field.minvalue);
writeRounded(fp, field.stepsize);

fprintf(fp, ' %d %d %d %d %d\n',field.Ggroup, field.Dgroup, field.protection, field.active, field.intptr);

switch field.subtype
   case {1, 3, 5, 6, 7}
      % Write numeric data
      buffer = makeoutput(field.value);
      fprintf(fp, '%s \n0 \n', buffer);
   case {2, 4}
      % Write string data
      loop = length(field.value);
      fprintf(fp, '%d ', loop);
      for j = 1:loop
         string = field.value{j};
         idx = strfind(string, '"');
         for k = idx(end:-1:1)
            string = [string(1:k-1) '\' string(k:end)];
         end
         fprintf(fp, '"%s"\n', string);
      end

      if isfield(field, 'enumeral')
         loop = length(field.enumeral);
         fprintf(fp, '%d', loop);
         for j = 1:loop
            fprintf(fp, ' "%s"', field.enumeral{j});
         end
         fprintf(fp,' \n');
      else
         fprintf(fp, '0 \n');
      end
end
return

%-------------------------------------------------------
function buffer = makeoutput(value)
% Make a string output for all the values
n = length(value);
buffer = sprintf('%d', n);
for m = 1:n
   if value(m) ~= 0
      p = floor(abs(log10(value(m))));
      if (p <= -5) || (p >= 12)
         buffer = [buffer sprintf([' %0.' num2str(numprecision(value(m)*10^(-p))) 'e'], value(m))]; %#ok<AGROW>
      else
         buffer = [buffer sprintf([' %0.' num2str(numprecision(value(m))) 'f'], value(m))]; %#ok<AGROW>
      end
   else
      buffer = [buffer ' 0']; %#ok<AGROW>
   end
end
return


function p = numprecision(v)
% Make the number positive, precision doesn't care about sign
v = abs(v);
poffset = 0;
if v >= 1
   % Strip out everything outside +/- 1
   remove = floor(v + eps);
   poffset = ceil(log10(remove));
   v = abs(v - remove);
end

% If the remainder is less than eps,
% then there is nothing left
if v < eps*10^poffset
   p = 0;
   return
end

% Strip out any zeros
pstart = -1*floor(log10(v));

% Maximum precision of 12 digits
for p = pstart:12
   check = v.*10^(p);
   check = check - floor(check + eps*10^(p+poffset));
   if abs(check) < 10^(p+poffset)*eps
      break
   end
end
return

function writeRounded(fp, value)
if value(1) ~= 0
   p = floor(abs(log10(value(1))));
   if (p <= -5) || (p >= 12)
      fprintf(fp, [' %0.' num2str(numprecision(value(1)*10^(-p))) 'e'], value(1));
   else
      fprintf(fp, [' %0.' num2str(numprecision(value(1))) 'f'], value(1));
   end
else
   fprintf(fp, ' 0');
end
return
%-------------------------------------------------------
