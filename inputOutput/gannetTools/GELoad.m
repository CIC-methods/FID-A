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
        
    case '24'
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
        
    case '20.007'
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
        
     case '20.006'
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
        
    case '26.002'
        
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

hdr.Larmor = i_hdr_value(rdb_hdr_ps_mps_freq)/1e7;
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
if nechoes == 1
    
    ShapeData = reshape(raw_data, [2 npoints totalframes nreceivers]);
    WaterData = ShapeData(:,:,2:refframes+1,:);
    FullData = ShapeData(:,:,refframes+2:end,:);
    
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



