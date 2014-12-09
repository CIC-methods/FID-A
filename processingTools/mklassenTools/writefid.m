% WRITEFID write k-space data to Varain FID fromat file
%
% WRITEFID(FPATH, K, PAR, HDR, BLOCK_HDR, CLEANUP) writes
% the k-space data K into the FID directory FPATH. PAR
% contains the acquisition parameters. HDR is the main 
% header which is automatically updated for size and datatype.
% BHDR is one block header. CLEANUP is a flag to remove
% original data from the write location (default true).
%
% see also READFID, WRITEPROCPAR, READPROCPAR

function writefid(fpath, k, par, hdr, bh, cleanup)

if nargin < 5
   error('cfmm:fid:arguments','WRITEFID requires five input arguments');
end

if nargin < 6
   cleanup = true;
end

% k must be 2 dimensional with blocks as the 2nd dimension
if ~ismatrix(k)
   error('cfmm:fid:arguments','WRITEFID requires k to be two dimensional');
end

[m,n] = size(k);

% Convert complex data to pairs of real numbers
k = k(:);
k = permute([real(k) imag(k)], [2,1]);

% Scale by a power of 2 such that data just fits into requested format
% abs is okay here because the real and imaginary channels are seperated
[f,p] = log2(max(abs(k(:)))); %#ok<ASGLU>

% Update the header based on the k-space and parameters
hdr.nblocks = n;
hdr.ntraces = 2 * m / par.np;
hdr.np = par.np;
% Get the data format from the header
if bitand(hdr.status, 8)
   % Floating
   hdr.ebytes = 4;
   writetype = [num2str(2*m) '*float32'];
   scale = 1;
else
   % Integer
   if bitand(hdr.status, 4)
      % 32-bit integer
      hdr.ebytes = 4;
      par.dp = 'y';
      writetype = [num2str(2*m) '*int32'];
      scale = 2^(31-p);
   else
      % 16-bit integer
      par.dp = 'n';
      hdr.ebytes = 2;
      writetype = [num2str(2*m) '*int16'];
      scale = 2^(15-p);
   end
end
hdr.tbytes = hdr.ebytes * hdr.np;
hdr.bbytes = 2 * m * hdr.ebytes + 28;
% hdr.vers_id = 
% hdr.status =
hdr.nbheaders = 1;
bh = bh(1);
clear f p;

% Scale the data ONLY if it is too big, i.e. data needs scaling down
if (scale < 1)
   k = k.*scale;
   % If we scale data, lvl and tlt should also be scaled
   bh.lvl = bh.lvl.*scale;
   bh.tlt = bh.tlt.*scale;
end
   
% Write out the new fid file 
fileid = fopen([fpath '/fid.out'], 'w+', 'ieee-be');
if fileid == -1
   error('cfmm:fid:output',['WRITEFID: Unable to create swap file ' fpath '/fid.out']);
end

try
   % Speed requires minimum possible fwrite statements and stream pointer
   % relocations. The compromise achieved below seems to offer the best 
   % performance
	% Write out the header file
	if fwrite(fileid, [hdr.nblocks hdr.ntraces hdr.np hdr.ebytes hdr.tbytes hdr.bbytes], 'int32') ~= 6
      error('cfmm:fid:output','WRITEFID: Error writing main header structure information.');
	end
	if fwrite(fileid, [hdr.vers_id hdr.status], 'int16') ~= 2
      error('cfmm:fid:output','WRITEFID: Error writing main header version information.');
	end
	if fwrite(fileid, hdr.nbheaders, 'int32') ~= 1
      error('cfmm:fid:output','WRITEFID: Error writing nbheaders field in main header.');
	end
   
   % Combine items of like type in block header
   bdata = zeros(4, hdr.nblocks);
   bdata(1,:) = bh.scale;
   bdata(2,:) = bh.status;
   bdata(3,:) = 1:hdr.nblocks;
   bdata(4,:) = bh.mode;
   if fwrite(fileid, bdata(:,1), 'int16') ~= 4
      error('cfmm:fid:output','WRITEFID: Error writing first block header integer values.');
   end
   if fwrite(fileid, bdata(:,2:end), '4*int16', hdr.bbytes-8) ~= (hdr.nblocks-1)*4
      error('cfmm:fid:output','WRITEFID: Error writing block header integer values.');
   end
   
	fseek(fileid, 40, 'bof');
	if fwrite(fileid, bh.ctcount, 'int32') ~= 1
      error('cfmm:fid:output','WRITEFID: Error writing ctcount field in first block header.');
	end
	if fwrite(fileid, bh.ctcount(ones(hdr.nblocks-1,1)), 'int32', hdr.bbytes-4) ~= hdr.nblocks-1
      error('cfmm:fid:output','WRITEFID: Error writing ctcount field in block headers.');
	end
   
   fseek(fileid, 44, 'bof');
   bdata(1,:) = bh.lpval;
   bdata(2,:) = bh.rpval;
   bdata(3,:) = bh.lvl;
   bdata(4,:) = bh.tlt;
   if fwrite(fileid, bdata(:,1), 'float32') ~= 4
      error('cfmm:fid:output','WRITEFID: Error writing first block header float values.');
   end
   if fwrite(fileid, bdata(:,2:end), '4*float32', hdr.bbytes-16) ~= (hdr.nblocks-1)*4
      error('cfmm:fid:output','WRITEFID: Error writing block header float values.');
   end
   clear bdata;
   
	% Goto the beginning of the data
	if fseek(fileid, 32, 'bof') == -1
      error('cfmm:fid:output','WRITEFID: Failure to move to start of data location.');
	end
	
	% Write out the data
	count = fwrite(fileid, k, writetype, 28);
catch err
   fclose(fileid);
   rethrow(err);
end

% Shut the file down 
fclose(fileid);

% Check that all the data was actually written
if numel(k) ~= count
   error('cfmm:fid:output','WRITEFID: Failure to write data.');
end

% Move the data file to its final location
%-----------------------------------------------------
% Comment out next line to remove backup fid file
if ~cleanup && (exist([fpath '/fid'], 'file') == 2)
   if unix(['mv -f ' fpath '/fid ' fpath '/fid.orig'])
      error('cfmm:fid:output','Unable to save original %s/fid', fpath);
   end
end
%-----------------------------------------------------
if unix(['mv -f ' fpath '/fid.out ' fpath '/fid'])
   error('cfmm:fid:output','Unable to replace original %s/fid', fpath);
end

return
