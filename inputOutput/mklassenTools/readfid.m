% readfid allows one to read in the Varian 4T FID data acquired with:
%
% [K, HDR, BLOCK_HDR] = READFID(FPATH, PAR, DC_CORRECT) reads k-space data 
% from fid file if the FID directory FPATH. DC correction based on level
% and tilt can be can be suppressed by setting DC_CORRECT to false. If DC
% correction is applied, the level and tilt are set to zero. PAR is the 
% parameters as read by getparam. The main header is returned in HDR and 
% the block headers are returned in BLOCK_HDR
%
% see also WRITEFID, READPROCPAR, WRITEPROCPAR

function [k, hdr, block_hdr, par] = readfid(fpath, par, dc_correct, tfisp, ext)

if nargin < 1
   error('cfmm:fid:name','READFID requires fid file name');
end

if nargin < 5
   ext = [];
end

% Read procpar if not already provided
if nargin < 2 || isempty(par)
   par = readprocpar(fpath, ext);
end

% Turn on DC correction by default
if nargin < 3 || isempty(dc_correct)
   dc_correct = true;
end
% Turn off tfisp by default
if nargin < 4 || isempty(tfisp)
   tfisp = false;
end

% Check for validity of fid file
if exist([fpath '/fid' ext], 'file') ~= 2
   if exist([fpath '.fid/fid'], 'file') ~= 2
      error('cfmm:fid:name','Unable to open FID file: %s', fpath);
   else
      fpath = [fpath '.fid'];
   end
end

infid = fopen([fpath '/fid' ext],'r', 'ieee-be');
if infid == -1
   error('cfmm:fid:name','Unable to open FID file: %s', fpath);
end

% Read in all the information from the FID file
try
   hdr.nblocks    = fread(infid, 1, 'int32');
   hdr.ntraces    = fread(infid, 1, 'int32');
   hdr.np         = fread(infid, 1, 'int32');
   hdr.ebytes     = fread(infid, 1, 'int32');
   hdr.tbytes     = fread(infid, 1, 'int32');
   hdr.bbytes     = fread(infid, 1, 'int32'); 
   hdr.vers_id    = fread(infid, 1, 'int16');
   hdr.status     = fread(infid, 1, 'int16');
   hdr.nbheaders  = fread(infid, 1, 'int32');
   
   m = hdr.bbytes/hdr.ebytes;
   n = hdr.nblocks;
   
   % Read in the block headers
   % Imaging data all has only one block header and VERY BAD people use
   % nbheaders as their own spare field. Therefore we ignore the nbheaders
   % and treat everything as having one block header.
   if nargout > 2
      % Create the block_hdr structure
      block_hdr.scale   = 0;
      block_hdr.status  = 0;
      block_hdr.index   = 0;
      block_hdr.mode    = 0;
      block_hdr.ctcount = 0;
      block_hdr.lpval   = 0;
      block_hdr.rpval   = 0;
      block_hdr.lvl     = 0;
      block_hdr.tlt     = 0;
      block_hdr = block_hdr(ones(n, 1));

      % Read in the block header values
      for i = 1:n
         block_hdr(i).scale   = fread(infid, 1, 'short');
         block_hdr(i).status  = fread(infid, 1, 'short');
         block_hdr(i).index   = fread(infid, 1, 'short');
         block_hdr(i).mode    = fread(infid, 1, 'short');
         block_hdr(i).ctcount = fread(infid, 1, 'long');
         block_hdr(i).lpval   = fread(infid, 1, 'float');
         block_hdr(i).rpval   = fread(infid, 1, 'float');
         block_hdr(i).lvl     = fread(infid, 1, 'float');
         block_hdr(i).tlt     = fread(infid, 1, 'float');
         fseek(infid, hdr.bbytes-28*hdr.nbheaders, 'cof');
      end
   end
   
   % Get the data format from the header
   if bitand(hdr.status, 8)
      % Floating
      format = '*float32';
      if hdr.ebytes ~= 4
         error('cfmm:fid:readError','Inconsistent Format: %d-bytes and %s.', hdr.ebytes, format);
      end
   else
      % Integer
      if bitand(hdr.status, 4)
         % 32-bit integer
         format = '*int32';
         if hdr.ebytes ~= 4
            error('cfmm:fid:readError','Inconsistent Format: %d-bytes and %s.', hdr.ebytes, format);
         end
      else
         % 16-bit integer
         format = '*int16';
         if hdr.ebytes ~= 2
            error('cfmm:fid:readError','Inconsistent Format: %d-bytes and %s.', hdr.ebytes, format);
         end
      end
   end
   
   % Read in all the data
   fseek(infid, 32, 'bof');
   k = fread(infid, [m,n], format);
   
   % Strip block header and reorganize into complex values
   k = complex(k(28/hdr.ebytes+1:2:m,:), k(28/hdr.ebytes+2:2:m,:));
   
   % Apply the lvl and tlt correction
   if dc_correct && rem(par.nt(1),2)
      % Read in the lvl and tlt correction
      fseek(infid, 52, 'bof');
      cor = fread(infid, [2, n], '2*float32', (hdr.bbytes) - 8);
      
      % Set the block header values to zero after correcting
      if (nargout > 1)
         for i = 1:n
            block_hdr(i).lvl     = 0;
            block_hdr(i).tlt     = 0;
         end
      end
      
      % nt = 3,7,11,... have DC rotated by 90 degrees
      if rem(par.nt,4) == 3
         cor = complex(-cor(2,:), cor(1,:));
      else
         cor = complex(cor(1,:), cor(2,:));
      end
      if tfisp
         % tfisp switches detection direction on each phase encode
         % acquisition is always compressed phase encodes
         [m, n] = size(k);
         nro = par.np / 2;
         if par.nv == 0
            nv = 1;
         else
            nv = par.nv;
         end
         r = m / (nro * nv);
         k = reshape(k, [nro nv r n]);
         if any(cor ~= cor(1))
            cor = reshape(cor, [1 1 1 n]); 
            nv_1 = ceil(nv / 2);
            nv_2 = floor(nv / 2);
            if rem(par.ssc, 2) == 0
               k(:,1:2:end,:,:) = k(:,1:2:end,:,:) + cor(ones(1,nro), ones(1,nv_1), ones(1,r), :); 
               k(:,2:2:end,:,:) = k(:,2:2:end,:,:) - cor(ones(1,nro), ones(1,nv_2), ones(1,r), :); 
            else
               k(:,2:2:end,:,:) = k(:,2:2:end,:,:) + cor(ones(1,nro), ones(1,nv_2), ones(1,r), :); 
               k(:,1:2:end,:,:) = k(:,1:2:end,:,:) - cor(ones(1,nro), ones(1,nv_1), ones(1,r), :); 
            end
            % cor = reshape(cor, [1 n]);
         else
            if rem(par.ssc, 2) == 0
               k(:,1:2:end,:) = k(:,1:2:end,:) + cor(1);
               k(:,2:2:end,:) = k(:,2:2:end,:) - cor(1);
            else
               k(:,2:2:end,:) = k(:,2:2:end,:) + cor(1);
               k(:,1:2:end,:) = k(:,1:2:end,:) - cor(1);
            end
         end
         k = reshape(k, [m n]);
      else
         if any(cor ~= cor(1))
            k = k - cor(ones(1,size(k,1)),:);
         else
            k = k - cor(1);
         end
      end
   end
   
catch err
    fclose(infid);
    rethrow(err);
end
fclose(infid);
