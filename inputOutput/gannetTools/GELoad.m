% GELoad.m
%
% USAGE:
% [FullData,WaterData,hdr]=GELoad(fname);
%
% DESCRIPTION:
% Reads in GE P-files (.7 file) using code adapted from GERead.m, provided
% as part of the Gannet software package by Richard Edden (gabamrs.com).
% Developed by Ralph Noeske and Mark Mikkelsen.
%
% INPUTS:
% filename   = filename of GE P file to be loaded.

function [FullData,WaterData,hdr]=GELoad(fname)

fid = fopen(fname, 'r', 'ieee-be');
fseek(fid, 0, 'bof');
rdbm_rev_num = fread(fid, 1, 'real*4');
if rdbm_rev_num == 7.0
    pfile_header_size = 39984; % LX
elseif rdbm_rev_num == 8.0
    pfile_header_size = 60464; % Cardiac / MGD
elseif rdbm_rev_num > 5.0 && rdbm_rev_num < 6.0
    pfile_header_size = 39940; % Signa 5.5
else
    % In 11.0 and later the header and data are stored as little-endian
    fclose(fid);
    fid = fopen(fname, 'r', 'ieee-le');
    fseek(fid, 0, 'bof');
    rdbm_rev_num = fread(fid, 1, 'real*4');
    if rdbm_rev_num == 9.0 % 11.0 product release
        pfile_header_size = 61464;
    elseif rdbm_rev_num == 11.0 % 12.0 product release
        pfile_header_size = 66072;
    end
end

% RTN 2018
% Added flexible P-file revision support
% Values are read from rdb_hdr and image sub-headers
% Position can be found in rdbm.h (RDB_HEADER_REC) and imagedb.h (MRIMAGEDATATYPE)

% RTN 2018
% unsigned int rdb_hdr_ps_mps_freq
% float rdb_hdr_user0
% float rdb_hdr_user4
% float rdb_hdr_user19
% short rdb_hdr_nechoes
% short rdb_hdr_navs
% short rdb_hdr_nframes
% short rdb_hdr_point_size
% unsigned short rdb_hdr_da_xres
% short rdb_hdr_da_yres
% short rdb_hdr_dab[0].start_rcv
% short rdb_hdr_dab[0].stop_rcv
% int rdb_hdr_off_image
% int rdb_hdr_off_data
%
% image sub-header
% int te
% int tr
% float user8-10    voxel dimensions
% float user19      rf waveform
% float user20-21   offset frequencies
% float user22      pulse width (-1 default)

switch num2str(rdbm_rev_num)

    case '14.3'

        % int
        rdb_hdr_off_image   = 377;
        rdb_hdr_off_data    = 368;
        rdb_hdr_ps_mps_freq = 107;

        % float
        rdb_hdr_user0  = 55;
        rdb_hdr_user4  = 59;
        rdb_hdr_user19 = 74;

        % short
        rdb_hdr_nechoes       = 36;
        rdb_hdr_navs          = 37;
        rdb_hdr_nframes       = 38;
        rdb_hdr_point_size    = 42;
        rdb_hdr_da_xres       = 52;
        rdb_hdr_da_yres       = 53;
        rdb_hdr_dab_start_rcv = 101;
        rdb_hdr_dab_stop_rcv  = 102;

        % int
        image_te = 181;
        image_tr = 179;

        % float
        image_user8  = 38;
        image_user19 = 49;
        image_user20 = 50;
        image_user22 = 52;

    case '16'

        % int
        rdb_hdr_off_image   = 377;
        rdb_hdr_off_data    = 368;
        rdb_hdr_ps_mps_freq = 107;

        % float
        rdb_hdr_user0  = 55;
        rdb_hdr_user4  = 59;
        rdb_hdr_user19 = 74;

        % short
        rdb_hdr_nechoes       = 36;
        rdb_hdr_navs          = 37;
        rdb_hdr_nframes       = 38;
        rdb_hdr_point_size    = 42;
        rdb_hdr_da_xres       = 52;
        rdb_hdr_da_yres       = 53;
        rdb_hdr_dab_start_rcv = 101;
        rdb_hdr_dab_stop_rcv  = 102;

        % int
        image_te = 193;
        image_tr = 191;

        % float
        image_user8  = 50;
        image_user19 = 61;
        image_user20 = 62;
        image_user22 = 64;

    case {'20.006','20.007','24'}

        % int
        rdb_hdr_off_image   = 377;
        rdb_hdr_off_data    = 368;
        rdb_hdr_ps_mps_freq = 107;

        % float
        rdb_hdr_user0  = 55;
        rdb_hdr_user4  = 59;
        rdb_hdr_user19 = 74;

        % short
        rdb_hdr_nechoes       = 36;
        rdb_hdr_navs          = 37;
        rdb_hdr_nframes       = 38;
        rdb_hdr_point_size    = 42;
        rdb_hdr_da_xres       = 52;
        rdb_hdr_da_yres       = 53;
        rdb_hdr_dab_start_rcv = 101;
        rdb_hdr_dab_stop_rcv  = 102;

        % int
        image_te = 267;
        image_tr = 265;

        % float
        image_user8  = 98;
        image_user19 = 109;
        image_user20 = 110;
        image_user22 = 112;

    case {'26.002','27','27.001','28.002','28.003','30','30.1'}

        % int
        rdb_hdr_off_image   = 11;
        rdb_hdr_off_data    = 2;
        rdb_hdr_ps_mps_freq = 123;

        % float
        rdb_hdr_user0  = 71;
        rdb_hdr_user4  = 75;
        rdb_hdr_user19 = 90;

        % short
        rdb_hdr_nechoes       = 74;
        rdb_hdr_navs          = 75;
        rdb_hdr_nframes       = 76;
        rdb_hdr_point_size    = 80;
        rdb_hdr_da_xres       = 90;
        rdb_hdr_da_yres       = 91;
        rdb_hdr_dab_start_rcv = 133;
        rdb_hdr_dab_stop_rcv  = 134;

        % int
        image_te = 267;
        image_tr = 265;

        % float
        image_user8  = 98;
        image_user19 = 109;
        image_user20 = 110;
        image_user22 = 112;

end


% Read rdb header as short, int and float
fseek(fid, 0, 'bof');
hdr_value = fread(fid, rdb_hdr_dab_stop_rcv, 'integer*2');
fseek(fid, 0, 'bof');
f_hdr_value = fread(fid, rdb_hdr_user19, 'real*4');
fseek(fid, 0, 'bof');
i_hdr_value = fread(fid, max(rdb_hdr_off_image, rdb_hdr_ps_mps_freq), 'integer*4');

if rdbm_rev_num > 11.0
    pfile_header_size = i_hdr_value(rdb_hdr_off_data);
end

%FINDING THE LARMOR FREQUENCY / B0:
% There is a 'probe-s' sequence in which the original code did not
% successfully retrieve the larmor frequency.  So here I will use a 
% different method if the seuqence is the 'probe-s' version.  This section
% is long because it takes some time to find the pulse sequence name (psd_nam):
%**************************************
SIZE_CHAR   = 1;
SIZE_IMTX   = 2;
SIZE_SHRT   = 2;
SIZE_ATMC   = 4;
SIZE_INT    = 4;
SIZE_FLT    = 4;
SIZE_DTYP   = 4;
SIZE_PTYP   = 4;
SIZE_RSPT   = 4;
SIZE_DBL    = 8;
SIZE_ASUB   = 48;

BYTE_ORDER = 'ieee-le';

G_to_T      = 1e-4;
MM_to_CM    = 0.1;
US_to_S     = 1e-6;
US_to_MS    = 1e-3;

%OFFSETS AND VERSION
fseek(fid,0,'eof');
off_end=ftell(fid);

fseek(fid,0,'bof');
verF=fread(fid,1,'float');
fseek(fid,0,'bof');
verI=fread(fid,1,'int');

off_data=fread(fid,1,'int');

off_pass=fread(fid,1,'int');

off = 4 * SIZE_INT;
fseek(fid,off,'cof');
off_tool=fread(fid,1,'int');

off_exm=fread(fid,1,'int');

off_ser=fread(fid,1,'int');

off_img=fread(fid,1,'int');

%IMAGE INFO:
off = off_img + SIZE_ASUB + 32 * SIZE_DBL;
fseek(fid,off,'bof');
fov = fread(fid,1,'float') * MM_to_CM;

off = SIZE_FLT;
fseek(fid,off,'cof');
t_acq = fread(fid,1,'float') * US_to_S;

slthk = fread(fid,1,'float');

off = 5 * SIZE_FLT;
fseek(fid,off,'cof');
sar_avg = fread(fid,1,'float');

sar_pk = fread(fid,1,'float');

off = 50 * SIZE_FLT;
fseek(fid,off,'cof');
ptl_dt = fread(fid,1,'float');

ptl_ns = fread(fid,1,'float');

n_x = fread(fid,1,'float');

n_ptl_tr = fread(fid,1,'float');

n_ptl_fov = fread(fid,1,'float');

n_avg = fread(fid,1,'float');

ltavg = fread(fid,1,'float');

exc_fa = fread(fid,1,'float');

off = 19 * SIZE_FLT + 24 * SIZE_RSPT + 2 * SIZE_DTYP + 2 * SIZE_PTYP ...
    + 30 * SIZE_FLT + 2 * SIZE_ATMC + 40 * SIZE_INT;
fseek(fid,off,'cof');
t_rep = fread(fid,1,'int') * US_to_MS;

off = SIZE_INT;
fseek(fid,off,'cof');
t_ech = fread(fid,1,'int') * US_to_MS;

off = 75 * SIZE_INT + 2 * SIZE_IMTX + 3 * SIZE_SHRT;
fseek(fid,off,'cof');
n_ech = fread(fid,1,'uint16');

off = 126 * SIZE_SHRT;
fseek(fid,off,'cof');
psd_nam = fread(fid,[1 33]*SIZE_CHAR,'char');
psd_nam = char(psd_nam(psd_nam>0));
%**************************************

%OK now we have the pulse sequence name (psd_nam).  If psd_nam is
%'probe-s', we will get the Larmor frequeny and B0 values a new way.
%Otherwise, we will just do it the old way.

if strcmp(psd_nam,'probe-s')
    off = off_exm + 32 * SIZE_DBL + 36 * SIZE_FLT + 11 * SIZE_ATMC + 32 * SIZE_INT;
    fseek(fid,off,'bof');
    Bo = fread(fid,1,'int') * G_to_T; %Field strength
    hdr.Larmor = Bo*42.577;
else
    hdr.Larmor = i_hdr_value(rdb_hdr_ps_mps_freq)/1e7;
end

hdr.sw = f_hdr_value(rdb_hdr_user0);

nechoes = hdr_value(rdb_hdr_nechoes);
nex = hdr_value(rdb_hdr_navs);
nframes = hdr_value(rdb_hdr_nframes);
point_size = hdr_value(rdb_hdr_point_size);
npoints = hdr_value(rdb_hdr_da_xres);
nrows = hdr_value(rdb_hdr_da_yres);

start_recv = hdr_value(rdb_hdr_dab_start_rcv);
stop_recv = hdr_value(rdb_hdr_dab_stop_rcv);
nreceivers = (stop_recv - start_recv) + 1;

% RTN 2018
dataframes = f_hdr_value(rdb_hdr_user4)/nex;
refframes = f_hdr_value(rdb_hdr_user19);

% Read image header as int and float
% MM (170118): Find TE/TR
fseek(fid, i_hdr_value(rdb_hdr_off_image), 'bof');
t_hdr_value = fread(fid, image_te, 'integer*4');
hdr.TE = t_hdr_value(image_te)/1e3;
hdr.TR = t_hdr_value(image_tr)/1e3;

% Spectro prescan pfiles
if npoints == 1 && nrows == 1
    npoints = 2048;
end

% Compute size (in bytes) of data
data_elements = npoints * 2;
totalframes = nrows * nechoes; % RTN nechoes mulitply
data_elements = data_elements * totalframes * nreceivers;

fseek(fid, pfile_header_size, 'bof');
% Read data: point_size = 2 means 16-bit data, point_size = 4 means EDR
if point_size == 2
    raw_data = fread(fid, data_elements, 'integer*2');
else
    raw_data = fread(fid, data_elements, 'integer*4');
end
fclose(fid);

% 110303 CJE
% Calculate Navg from nframes, 8 water frames, 2 phase cycles
% Needs to be specific to single experiment - for frame rejection
% RTN edits to accommodate Noeske version raee 20141007
% MM (160916): Incorporating more edits from RTN to handle dual-echo data
%              acquired with one of four possible encoding schemes:
%              NEX=2/noadd=0, NEX=2/noadd=1, NEX=8/noadd=0, NEX=8/noadd=1
% MM (171120): RTN edits to accomodate HERMES aquisitions; better looping
%              over phase cycles
% MM (200713): RTN edits for proper handling of no_add if nechoes == 1
if nechoes == 1
    
    if (dataframes + refframes) ~= nframes
        mult = 1;
        dataframes = dataframes * nex;
        refframes = nframes - dataframes;
    else
        mult = 1/nex;
    end
    
    ShapeData = reshape(raw_data, [2 npoints totalframes nreceivers]);
    WaterData = ShapeData(:,:,2:refframes+1,:) * mult;
    FullData = ShapeData(:,:,refframes+2:end,:) * mult;
    
    totalframes = totalframes - (refframes+1);
    waterframes = refframes;
    
else
    
    if (dataframes + refframes) ~= nframes
        mult = nex/2; % RTN 2016
        multw = nex; % RTN 2016
        %mult = 1; % RTN 2017
        %multw = 1; % RTN 2017
        noadd = 1;
        dataframes = dataframes * nex;
        refframes = nframes - dataframes;
    else
        mult = nex/2; % RTN 2016
        multw = 1; % RTN 2016
        %mult = 1; % RTN 2017
        %multw = 1/nex; % RTN 2017
        noadd = 0;
    end
    
    if totalframes ~= (dataframes + refframes + 1) * nechoes % RTN 2017
        error('# of totalframes not same as (dataframes + refframes + 1) * nechoes');
    end
    
    ShapeData = reshape(raw_data, [2 npoints totalframes nreceivers]);
    
    % MM (180404)
    [X1,X2] = ndgrid(1:refframes, 1:nechoes);
    X1 = X1'; X1 = X1(:);
    X2 = X2'; X2 = X2(:);
    Y1 = (-1).^(noadd * (X1-1));
    Y1 = permute(repmat(Y1, [1 npoints 2 nreceivers]), [3 2 1 4]);
    Y2 = 1 + (totalframes/nechoes) * (X2-1) + X1;
    WaterData = Y1 .* ShapeData(:,:,Y2,:) * multw;
    
    [X1,X2] = ndgrid(1:dataframes, 1:nechoes);
    X1 = X1'; X1 = X1(:);
    X2 = X2'; X2 = X2(:);
    Y1 = (-1).^(noadd * (X1-1));
    Y1 = permute(repmat(Y1, [1 npoints 2 nreceivers]), [3 2 1 4]);
    Y2 = 1 + refframes + (totalframes/nechoes) * (X2-1) + X1;
    FullData = Y1 .* ShapeData(:,:,Y2,:) * mult;
    
    totalframes = totalframes - (refframes + 1) * nechoes; % RTN 2017
    waterframes = refframes * nechoes; % RTN 2017
    
end

FullData = FullData .* repmat([1; 1i], [1 npoints totalframes nreceivers]);
FullData = squeeze(sum(FullData,1));
WaterData = WaterData .* repmat([1; 1i], [1 npoints waterframes nreceivers]);
WaterData = squeeze(sum(WaterData,1));



