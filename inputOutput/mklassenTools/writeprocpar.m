function writeprocpar(fpath, par, cleanup)

% WRITEPROCPAR writes in parameter used in FID acquisition
%
%   WRITEPROCPAR(FPATH,PAR) writes the parameter structure PAR, to the 
%   procpar in the FID directory FPATH. CLEANUP is a flag to remove 
%   previous data from write location (default true)
%
% see also READPROCPAR, WRITEFID, READFID 


% Written by L Martyn Klassen
% Copyright 2003 Robarts Research Institute

if nargin < 2
   error('CFMM:writeprocpar','Invalid arguments, requires path and parameter set.');
end

if nargin < 3
   cleanup = true;
end

% Check for validity of file
procparfile = strcat(fpath,'/procpar');
fp = fopen(procparfile, 'r+', 'ieee-be');
if fp == -1
   procparfile = strcat(fpath,'.fid/procpar');
   fp = fopen(procparfile, 'r+', 'ieee-be');
   if (fp == -1)
      procparfile = strcat(fpath,'.par/procpar');
      fp = fopen(procparfile, 'r+', 'ieee-be');
      if (fp == -1)
         error('CFMM:writeprocpar','Unable to open %s for writing.',fpath);
      end
   end
end

procparSwap = strcat(procparfile,'.out');
fp2 = fopen(procparSwap, 'w+', 'ieee-be');
if fp2 == -1
   fclose(fp);
   error('CFMM:writeprocpar','Unable to create swap file %',procparSwap);
end

% This function extracts the first letter of each field name
% in order to provide a quickly discard any line which does
% not being with one of the parameters of interest.
names = fieldnames(par);
value = names{1}(1);
for n = 2:size(names,1)
   if ~sum(value == names{n}(1))
      value = [value names{n}(1)]; %#ok<AGROW>
   end
end
clear n names;

buffer = fgets(fp);
 
% Parse the ASCII procpar file
while (buffer ~= -1)
   % Check to see if the first letter in the buffer matches
   % the first character of any parameter of interest
   % This provides a 3 to 4 fold speed up.
   if (strfind(value, buffer(1)))
      
      % Get only the first word of the buffer, the parameter name
      ind = strfind(buffer, ' ');
      lenb = ind(1)-1;
      cmpbuffer = buffer(1:lenb);
   
      % Read in required parameters
      if (lenb == 2)
         if any(strcmp(cmpbuffer, {'z1','z2','z3','z4','z5','z6', ...
               'z7','z8','x1','y1','xz','yz','xy','x3','y3','x4', ...
               'y4','nD','nf','ni','np','ns','nv','ne','ss','nt'}))
            fprintf(fp2, '%s', buffer);
            fgets(fp);
            buffer = ['1' sprintf(' %d', par.(cmpbuffer)(1))];
            fprintf(fp2, '%s \n', buffer);
            buffer = fgets(fp);
         elseif any(strcmp(cmpbuffer, {'ti','te','tr','sw','at', ...
               'bt','d1','SR','r1','r2','B0','SR','Po'}))
            fprintf(fp2, '%s', buffer);
            fgets(fp);
            buffer = makeoutput(par.(cmpbuffer));
            fprintf(fp2, '%s \n', buffer);
            buffer = fgets(fp);
         elseif (strcmp(cmpbuffer, {'dp','tn'}))
            fprintf(fp2, '%s', buffer);
            fgets(fp);
            fprintf(fp2, '1 "%s" \n', par.(cmpbuffer));
            buffer = fgets(fp);
         end
      elseif (lenb == 3),
         if any(strcmp(cmpbuffer, {'z1c','z2c','z3c','z4c','xz2', ...
               'yz2','zxy','z3x','z3y','zx3','zy3','z4x','z4y', ...
               'z5x','z5y','ssc','nv2','nv3','ne1','ne2'}))
            fprintf(fp2, '%s', buffer);
            fgets(fp);
            buffer = ['1' sprintf(' %d', par.(cmpbuffer)(1))];
            fprintf(fp2, '%s \n', buffer);
            buffer = fgets(fp);
         elseif any(strcmp(cmpbuffer, {'tep','esp','lro','pro', ...
               'lpe','ppe','pss','phi','psi','thk','gro','gss', ...
               'tof','fov','res','gro','etl','Psl','dtg','tpe'}))
            fprintf(fp2, '%s', buffer);
            fgets(fp);
            buffer = makeoutput(par.(cmpbuffer));
            fprintf(fp2, '%s \n', buffer);
            buffer = fgets(fp);
         elseif (strcmp(cmpbuffer, 'nav')) || ...
                  (strcmp(cmpbuffer, 'seg'))
            fprintf(fp2, '%s', buffer);
            fgets(fp);
            fprintf(fp2, '1 "%s" \n', par.(cmpbuffer));
            buffer = fgets(fp);
         end
      elseif (lenb == 4)
         if any(strcmp(cmpbuffer, {'x2y2','z2xy','z3xy','z2x3', ...
               'z2y3','z3x3','z3y3','z4xy'}))
            fprintf(fp2, '%s', buffer);
            fgets(fp);
            buffer = ['1' sprintf(' %d', par.(cmpbuffer)(1))];
            fprintf(fp2, '%s \n', buffer);
            buffer = fgets(fp);
         elseif any(strcmp(cmpbuffer, {'lpe2','ppe2','thk2','pos1', ...
               'pos2','pos3','vox1','vox2','vox3','vpsi','vphi', ...
               'gmax','flip','gss2','sfrq','nnav','tnav','grox', ...
               'groy','nseg','npix','gimp','grof','freq','tpwr', ...
               'dpwr'}))
            fprintf(fp2, '%s', buffer);
            fgets(fp);
            buffer = makeoutput(par.(cmpbuffer));
            fprintf(fp2, '%s \n', buffer);
            buffer = fgets(fp);
         elseif any(strcmp(cmpbuffer, {'gain','cntr'}))
            fprintf(fp2, '%s', buffer);
            fgets(fp);
            buffer = [sprintf('%d', length(par.(cmpbuffer))) sprintf(' %d', par.(cmpbuffer))];
            fprintf(fp2, '%s \n', buffer);
            buffer = fgets(fp);
         elseif (strcmp(buffer, 'date')) || ...
                  (strcmp(cmpbuffer, 'fast'))
            fprintf(fp2, '%s', buffer);
            fgets(fp);
            fprintf(fp2, '1 "%s" \n', par.(cmpbuffer));
            buffer = fgets(fp);
         end
      elseif (lenb == 5)
         if any(strcmp(cmpbuffer, {'zx2y2','celem'}))
            fprintf(fp2, '%s', buffer);
            fgets(fp);
            buffer = ['1' sprintf(' %d', par.(cmpbuffer)(1))];
            fprintf(fp2, '%s \n', buffer);
            buffer = fgets(fp);
         elseif any(strcmp(cmpbuffer, {'tproj','trise','tpwr1', ...
               'tpwr2','tpwri','theta','resto','nhomo','nfreq', ...
               'dtmap','nzseg','state','flip1','tpwrf','dpwrf'}))
            fprintf(fp2, '%s', buffer);
            fgets(fp);
            buffer = makeoutput(par.(cmpbuffer));
            fprintf(fp2, '%s \n', buffer);
            buffer = fgets(fp);
         elseif (strcmp(cmpbuffer, 'nproj'))
            fprintf(fp2, '%s', buffer);
            fgets(fp);
            buffer = [sprintf('%d', length(par.(cmpbuffer))) sprintf(' %d', par.(cmpbuffer))];
            fprintf(fp2, '%s \n', buffer);
            buffer = fgets(fp);
         elseif any(strcmp(cmpbuffer, {'rcvrs','intlv'}))
            fprintf(fp2, '%s', buffer);
            fgets(fp);
            fprintf(fp2, '1 "%s" \n', par.(cmpbuffer));
            buffer = fgets(fp);
         elseif (strcmp(cmpbuffer, 'array'))
            fprintf(fp2, '%s', buffer);
            fgets(fp);
            buffer = '1 "';
            for o = 1:length(par.array)
               if length(par.array{o}) > 1
                  buffer = [buffer '(']; %#ok<AGROW>
                  for p = 1:length(par.array{o})
                     buffer = [buffer par.array{o}{p} ',']; %#ok<AGROW>
                  end
                  buffer = [buffer(1:end-1) '),'];
               else
                  buffer = [buffer par.array{o}{1} ',']; %#ok<AGROW>
               end
            end
            buffer = [buffer(1:end-1) '"'];
            fprintf(fp2, '%s\n', buffer);
            buffer = fgets(fp);
         end
      elseif (lenb == 6)
         if any(strcmp(cmpbuffer, {'z2x2y2','z3x2y2','z4x2y2'}))
            fprintf(fp2, '%s', buffer);
            fgets(fp);
            buffer = ['1' sprintf(' %d', par.(cmpbuffer)(1))];
            fprintf(fp2, '%s \n', buffer);
            buffer = fgets(fp);
         elseif any(strcmp(cmpbuffer, {'vtheta','fpmult','nturns', ...
               'tfirst'}))
            fprintf(fp2, '%s', buffer);
            fgets(fp);
            buffer = makeoutput(par.(cmpbuffer));
            fprintf(fp2, '%s \n', buffer);
            buffer = fgets(fp);
         elseif any(strcmp(cmpbuffer, {'seqfil','seqcon','orient'}))
            fprintf(fp2, '%s', buffer);
            fgets(fp);
            fprintf(fp2, '1 "%s" \n', par.(cmpbuffer));
            buffer = fgets(fp);
         end
      elseif (lenb == 7)
         if any(strcmp(cmpbuffer, {'espincr','rfphase','rfdelay','cluster'}))
            fprintf(fp2, '%s', buffer);
            fgets(fp);
            buffer = makeoutput(par.(cmpbuffer));
            fprintf(fp2, '%s \n', buffer);
            buffer = fgets(fp);
         elseif (strcmp(cmpbuffer, 'shimset'))
            fprintf(fp2, '%s', buffer);
            fgets(fp);
            buffer = ['1' sprintf(' %d', par.(cmpbuffer)(1))];
            fprintf(fp2, '%s \n', buffer);
            buffer = fgets(fp);
         elseif any(strcmp(cmpbuffer, {'rfspoil','varflip','petable', ...
               'console','profile','bipolar'}))
            fprintf(fp2, '%s', buffer);
            fgets(fp);
            fprintf(fp2, '1 "%s" \n', par.(cmpbuffer));
            buffer = fgets(fp);
         end
      elseif (lenb == 8)
         if any(strcmp(cmpbuffer, {'evenecho','polarity'}))
            fprintf(fp2, '%s', buffer);
            fgets(fp);
            buffer = [sprintf('%d', length(par.(cmpbuffer))) sprintf(' %d', par.(cmpbuffer))];
            fprintf(fp2, '%s \n', buffer);
            buffer = fgets(fp);
         elseif any(strcmp(cmpbuffer, {'fliplist','gradfrac'}))
            fprintf(fp2, '%s', buffer);
            fgets(fp);
            buffer = makeoutput(par.(cmpbuffer));
            fprintf(fp2, '%s \n', buffer);
            buffer = fgets(fp);
         elseif (strcmp(cmpbuffer, 'arraydim'))
            fprintf(fp2, '%s', buffer);
            fgets(fp);
            buffer = ['1' sprintf(' %d', par.(cmpbuffer)(1))];
            fprintf(fp2, '%s \n', buffer);
            buffer = fgets(fp);
         elseif any(strcmp(cmpbuffer, {'waveform','contrast','flipprep',...
                 'readaxis','ky_order'}))
            fprintf(fp2, '%s', buffer);
            fgets(fp);
            fprintf(fp2, '1 "%s" \n', par.(cmpbuffer));
            buffer = fgets(fp);
         end
      elseif (lenb == 9)
         if any(strcmp(cmpbuffer, {'direction','weightfit'}))
            fprintf(fp2, '%s', buffer);
            fgets(fp);
            buffer = [sprintf('%d', length(par.(cmpbuffer))) sprintf(' %d', par.(cmpbuffer))];
            fprintf(fp2, '%s \n', buffer);
            buffer = fgets(fp);
         elseif any(strcmp(cmpbuffer, {'threshold'}))
            fprintf(fp2, '%s', buffer);
            fgets(fp);
            buffer = makeoutput(par.(cmpbuffer));
            fprintf(fp2, '%s \n', buffer);
            buffer = fgets(fp);
         elseif any(strcmp(cmpbuffer, {'timescale','navigator','alternate',...
               'weightfit'}))
            fprintf(fp2, '%s', buffer);
            fgets(fp);
            fprintf(fp2, '1 "%s" \n', par.(cmpbuffer));
            buffer = fgets(fp);
         end
      elseif (lenb == 10)
         if any(strcmp(cmpbuffer, {'spiral_tep','randomseed'}))
            fprintf(fp2, '%s', buffer);
            fgets(fp);
            buffer = makeoutput(par.(cmpbuffer));
            fprintf(fp2, '%s \n', buffer);
            buffer = fgets(fp);
         end 
      elseif (lenb == 11)
         if any(strcmp(cmpbuffer, {'spiral_gmax'}))
            fprintf(fp2, '%s', buffer);
            fgets(fp);
            buffer = makeoutput(par.(cmpbuffer));
            fprintf(fp2, '%s \n', buffer);
            buffer = fgets(fp);
         elseif (strcmp(cmpbuffer, {'arrayelemts','ninterleave'}))
            fprintf(fp2, '%s', buffer);
            fgets(fp);
            buffer = ['1' sprintf(' %d', par.(cmpbuffer)(1))];
            fprintf(fp2, '%s \n', buffer);
            buffer = fgets(fp);
         end
      elseif (lenb == 12)
         if any(strcmp(cmpbuffer, {'spiral_gamma','spiral_alpha'}))
            fprintf(fp2, '%s', buffer);
            fgets(fp);
            buffer = makeoutput(par.(cmpbuffer));
            fprintf(fp2, '%s \n', buffer);
            buffer = fgets(fp);
         end 
      elseif (lenb == 13)
         if any(strcmp(cmpbuffer, {'spiral_filter'}))
            fprintf(fp2, '%s', buffer);
            fgets(fp);
            buffer = makeoutput(par.(cmpbuffer));
            fprintf(fp2, '%s \n', buffer);
            buffer = fgets(fp);
         end 
      elseif (lenb == 14)
         if any(strcmp(cmpbuffer, {'spiral_version','spiral_density','offlineAverage',...
               'referenceUnits','referenceLines'}))
            fprintf(fp2, '%s', buffer);
            fgets(fp);
            buffer = makeoutput(par.(cmpbuffer));
            fprintf(fp2, '%s \n', buffer);
            buffer = fgets(fp);
         end 
      elseif (lenb == 15)
         if any(strcmp(cmpbuffer, {'offlineAverages','reductionFactor'}))
            fprintf(fp2, '%s', buffer);
            fgets(fp);
            buffer = makeoutput(par.(cmpbuffer));
            fprintf(fp2, '%s \n', buffer);
            buffer = fgets(fp);
         end 
      elseif (lenb == 16)
         if any(strcmp(cmpbuffer, {'interleave_order'}))
            fprintf(fp2, '%s', buffer);
            fgets(fp);
            fprintf(fp2, '1 "%s" \n', par.(cmpbuffer));
            buffer = fgets(fp);
         end 
      end   
   end
   
   % Write out the current line buffer
   fprintf(fp2, '%s', buffer);
   
   % Read in the next line
   buffer = fgets(fp);
end  % End of while loop

% Close file
fclose(fp);
fclose(fp2);

% Move the procpar file to its final location
%-------------------------------------------------------
% Comment out next line to remove backup of procpar
if ~cleanup && (exist(procparfile, 'file') == 2)
   if unix(['mv -f ' procparfile ' ' procparfile '.orig']);
      error('CFMM:writeprocpar','Unable to save original %s', procparfile);
   end
end
%-------------------------------------------------------
if unix(['mv -f ' procparSwap ' '  procparfile])
   error('CFMM:writeprocpar', 'Unable to replace original %s/procpar', procparfile);
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

%-------------------------------------------------------
