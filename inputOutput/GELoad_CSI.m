% GELoad_CSI.m
%
% USAGE:
% [FullData,hdr]=GELoad_CSI(fname);
%
% DESCRIPTION:
% Reads in GE P-files (.7 file) using code adapted from nGetData_dv26.py, 
% provided by Louis Lauzon.

% INPUTS:
% fname     = filename of GE P file to be loaded.

% OUTPUTS:
% FullData  = MRSI data matrix.
% hdr       = Struct of header information.

function [data,hdr]=GELoad_CSI(fname)

FROM_BEG = 0;
FROM_CUR = 1;

OFFSET_TOOLS    = 188424;
OFFSET_EXAM     = 193660;
OFFSET_SERIES   = 195620;
OFFSET_IMAGE    = 198180;
HDR_LEN         = 213684;

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

off = OFFSET_EXAM + 32 * SIZE_DBL + 36 * SIZE_FLT ...
    + 11 * SIZE_ATMC + 32 * SIZE_INT;

fid=fopen(fname,'r',BYTE_ORDER);
fseek(fid,off,'bof');

%Get a bunch of params for the header:
Bo=fread(fid,1,'int')*G_to_T; %Field strength

off = 30 * SIZE_INT;
fseek(fid,off,0);
exm_no = fread(fid,1,'uint16');

off = 51 * SIZE_SHRT + (257 + 3 * 65) * SIZE_CHAR;
fseek(fid,off,0);
exm_dsc = fread(fid,[1 65]*SIZE_CHAR,'char');
exm_dsc = char(exm_dsc(exm_dsc>0));

off = (3 + 17 + 13 + 34 + 8 + 32 + 4 + 3*32) * SIZE_CHAR;
fseek(fid,off,0);
pat_nam = fread(fid,[1 65]*SIZE_CHAR,'char');
pat_nam = char(pat_nam(pat_nam>0));

pat_id = fread(fid,[1 65]*SIZE_CHAR,'char');
pat_id = char(pat_id(pat_id>0));

%SERIES INFO:
off = OFFSET_SERIES + 32 * SIZE_DBL + 54 * SIZE_FLT + 7 * SIZE_ATMC ...
    + 59 * SIZE_INT;
fseek(fid,off,'bof');
ser_no = fread(fid,1,'int');

off = 17 * SIZE_INT + 76 * SIZE_SHRT + 22 * SIZE_CHAR;
fseek(fid,off,0);
ser_dsc = fread(fid,[1 65]*SIZE_CHAR,'char');
ser_dsc = char(ser_dsc(ser_dsc>0));

%TOOLS INFO:
off = OFFSET_TOOLS + 4 * SIZE_SHRT + 4 * SIZE_FLT + 49 * SIZE_INT;
fseek(fid,off,'bof');
coil = fread(fid,[1 16]*SIZE_CHAR,'char');
coil = char(coil(coil>0));

off = 4 * SIZE_CHAR + 2 * SIZE_INT + 2 * SIZE_SHRT + 4 * SIZE_FLT;
fseek(fid,off,0);
n_chn = fread(fid,1,'int');

%IMAGE INFO:
off = OFFSET_IMAGE + SIZE_ASUB + 32 * SIZE_DBL;
fseek(fid,off,'bof');
fov = fread(fid,1,'float') * MM_to_CM;

off = SIZE_FLT;
fseek(fid,off,0);
t_acq = fread(fid,1,'float') * US_to_S;

slthk = fread(fid,1,'float');

off = 5 * SIZE_FLT;
fseek(fid,off,0);
sar_avg = fread(fid,1,'float');

sar_pk = fread(fid,1,'float');

off = 50 * SIZE_FLT;
fseek(fid,off,0);
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
fseek(fid,off,0);
t_rep = fread(fid,1,'int') * US_to_MS;

off = SIZE_INT;
fseek(fid,off,0);
t_ech = fread(fid,1,'int') * US_to_MS;

off = 75 * SIZE_INT + 2 * SIZE_IMTX + 3 * SIZE_SHRT;
fseek(fid,off,0);
n_ech = fread(fid,1,'uint16');

%RECON_INFO
off = SIZE_FLT + 20 * SIZE_INT + SIZE_SHRT + 6 * SIZE_CHAR;
fseek(fid,off,'bof');
scn_date = fread(fid,[1 10]*SIZE_CHAR,'char');
scn_date = char(scn_date(scn_date>0));

scn_time = fread(fid,[1 8]*SIZE_CHAR,'char');
scn_time = char(scn_time(scn_time>0));

off = 10 * SIZE_CHAR + 12 * SIZE_SHRT;
fseek(fid,off,0);
n_slc = fread(fid,1,'uint16');

off = 2 * SIZE_SHRT;
fseek(fid,off,0);
adc_y = fread(fid,1,'int16');

off = 2 * SIZE_SHRT;
fseek(fid,off,0);
adc_x = fread(fid,1,'int16');

pt_siz = fread(fid,1,'int16');

hdr.pt_siz      = pt_siz;
hdr.adc_x       = adc_x;
hdr.adc_y       = adc_y;
hdr.n_slc       = n_slc;
hdr.scn_time    = scn_time;
hdr.scn_date    = scn_date;
hdr.n_ech       = n_ech;
hdr.t_ech       = t_ech;
hdr.t_rep       = t_rep;
hdr.exc_fa      = exc_fa;
hdr.ltavg       = ltavg;
hdr.n_avg       = n_avg;
hdr.n_ptl_fov   = n_ptl_fov;
hdr.n_ptl_tr    = n_ptl_tr;
hdr.n_x         = n_x;
hdr.ptl_ns      = ptl_ns;
hdr.ptl_dt      = ptl_dt;
hdr.sar_pk      = sar_pk;
hdr.sar_avg     = sar_avg;
hdr.slthk       = slthk;
hdr.t_acq       = t_acq;
hdr.fov         = fov;
hdr.n_chn       = n_chn;
hdr.coil        = coil;
hdr.ser_dsc     = ser_dsc;
hdr.ser_no      = ser_no;
hdr.pat_id      = pat_id;
hdr.pat_nam     = pat_nam;
hdr.exm_dsc     = exm_dsc;
hdr.exm_no      = exm_no;
hdr.Bo          = Bo;


%**********FINISHED READING HEADER INFO***********

%NOW TIME TO START READING THE DATA:
%ORDERING
n_1 = 2 * adc_x; %t
n_2 = 1 + adc_y; %PE 
n_3 = n_ech;     %ECHO
n_4 = n_slc;     %SLICE
n_5 = n_chn;     %COIL

n_tot = n_1 * n_2 * n_3 * n_4 * n_5;

%Start reading the data:
fseek(fid,HDR_LEN,'bof');
data = fread(fid,n_tot,'int32');

%reshape the data matrix:
data = reshape(data,[n_1,n_2,n_3,n_4,n_5]);

%Squeeze out the singleton dims
data=squeeze(data);

%remove the baseline:
data = data(:,2:end,:,:);

%Correct for RF chopping - negate every line along PE
data(:,2:2:end,:,:) = -1 * data(:,2:2:end,:,:);

%Combine the real and imaginary parts into a complex vector
data = data(1:2:end,:,:,:) + 1i * data(2:2:end,:,:,:);

%Concatenate the echoes onto eachother.  But before doing that, add one
%point to the end of each adc chunk to account for the fact that there is a
%dead time of one petal rotation (ptl_dt) between the three ADC chunks.
%The value of that point will be the average of the last point in the
%previous chunk and the first point in the subsequent chunk:
data(adc_x+1:adc_x+ptl_ns,:,1,:) = mean(cat(3,data(adc_x-(ptl_ns-1):adc_x,:,1,:),data(1:ptl_ns,:,2,:)),3);
data(adc_x+1:adc_x+ptl_ns,:,2,:) = mean(cat(3,data(adc_x-(ptl_ns-1):adc_x,:,2,:),data(1:ptl_ns,:,3,:)),3); 
data(adc_x+1:adc_x+ptl_ns,:,3,:) = zeros(ptl_ns,adc_y,1,n_chn);
data = reshape(permute(data,[1,3,2,4]),[(adc_x+ptl_ns)*n_ech,adc_y,n_chn]);

%update the number of loops per petal, and the number of ADC points:
n_ptl_tr        = n_ptl_tr+3;
hdr.n_ptl_tr    = hdr.n_ptl_tr+3;
hdr.adc_x       = hdr.adc_x + hdr.ptl_ns;
adc_x           = adc_x + ptl_ns;

%Separate the averages from the rosette petals:
if ltavg
    data=reshape(permute(data,[1,3,2]),[adc_x*n_ech,n_chn,n_ptl_fov,n_avg]);
elseif ~ltavg
    data=reshape(permute(data,[1,3,2]),[adc_x*n_ech,n_chn,n_avg,n_ptl_fov]);
end

%Trim off the last point in the time dimension, which are all
%zeros anyway:
%FullData=dat(1:end-1,:,:,:);

%Finally, write a struct to fill in the dims order.  Using the Siemens notation;
hdr.dims={'Col','Cha','Ave','Lin'};






