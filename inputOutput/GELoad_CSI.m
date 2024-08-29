% GELoad_CSI.m
% Jamie Near, Sunnybrook Research Institute 2024
% Louis Lauzon, University of Calgary 2024
%
% USAGE:
% [FullData,hdr]=GELoad_CSI(fname);
%
% DESCRIPTION:
% Reads in GE P-files (.7 file) using code adapted from nGetData_dv26.py
% and aux_get_hdr.py, both provided by Louis Lauzon.
%
% Alo reads in GE .h5 files (from MR30 version systems) using Orchestra (GE), 
% as well as code adapted from nGetData_mr30.py, also provided by Louis 
% Lauzon.

% INPUTS:
% fname     = filename of GE P file or .h5 file to be loaded.

% OUTPUTS:
% FullData  = MRSI data matrix.
% hdr       = Struct of header information.

function [data,hdr]=GELoad_CSI(fname)

if strcmp(fname(1),'P') && strcmp(fname(end-1:end),'.7')
    %This is a p-file.  Process accordingly
    HDR_LEN     = 213684;

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

    fid=fopen(fname,'r',BYTE_ORDER);

    %Retrieving parameters from the header:

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

    %RECON_INFO
    off = SIZE_FLT + 20 * SIZE_INT + SIZE_SHRT + 6 * SIZE_CHAR;
    fseek(fid,off,'bof');
    scn_date = fread(fid,[1 10]*SIZE_CHAR,'char');
    scn_date = char(scn_date(scn_date>0));

    scn_time = fread(fid,[1 8]*SIZE_CHAR,'char');
    scn_time = char(scn_time(scn_time>0));

    off = 10 * SIZE_CHAR + 12 * SIZE_SHRT;
    fseek(fid,off,'cof');
    n_slc = fread(fid,1,'uint16');

    off = 2 * SIZE_SHRT;
    fseek(fid,off,'cof');
    adc_y = fread(fid,1,'int16');

    off = 2 * SIZE_SHRT;
    fseek(fid,off,'cof');
    adc_x = fread(fid,1,'int16');

    pt_siz = fread(fid,1,'int16');

    %Re-populate n_slc for MR30
    if verF==30
        off = off_pass - 398 * SIZE_SHRT - 7 * SIZE_INT - 18 * SIZE_FLT;
        fseek(fid,off,'bof');
        n_slc=fread(fid,1,'uint16');
    end

    %TOOLS INFO:
    off = off_tool + 4 * SIZE_SHRT + 4 * SIZE_FLT + 49 * SIZE_INT;
    fseek(fid,off,'bof');
    coil = fread(fid,[1 16]*SIZE_CHAR,'char');
    coil = char(coil(coil>0));

    off = 4 * SIZE_CHAR + 2 * SIZE_INT + 2 * SIZE_SHRT + 4 * SIZE_FLT;
    fseek(fid,off,'cof');
    n_chn = fread(fid,1,'int');

    %EXAM INFO:
    off = off_exm + 32 * SIZE_DBL + 36 * SIZE_FLT + 11 * SIZE_ATMC + 32 * SIZE_INT;
    fseek(fid,off,'bof');
    Bo = fread(fid,1,'int') * G_to_T; %Field strength

    off = 30 * SIZE_INT;
    fseek(fid,off,'cof');
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
    off = off_ser + 32 * SIZE_DBL + 54 * SIZE_FLT + 7 * SIZE_ATMC ...
        + 59 * SIZE_INT;
    fseek(fid,off,'bof');
    ser_no = fread(fid,1,'int');

    off = 17 * SIZE_INT + 76 * SIZE_SHRT + 22 * SIZE_CHAR;
    fseek(fid,off,'cof');
    ser_dsc = fread(fid,[1 65]*SIZE_CHAR,'char');
    ser_dsc = char(ser_dsc(ser_dsc>0));

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

    %Calculate ACTUAL n_chn based on file size:
    val = off_end;
    val = val - off_data;
    val = val / (ptl_ns * 2 * pt_siz);
    val = val / n_ptl_tr;
    val = val / (adc_y +1);
    val = val / n_slc;
    n_chnA = val;


    %Populate the "hdr" sctructure:
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
    hdr.psd_nam     = psd_nam;
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
    n_5 = n_chnA;    %COIL


    n_tot = n_1 * n_2 * n_3 * n_4 * n_5;

    %Start reading the data:
    fseek(fid,off_data,'bof');
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
    data(adc_x+1:adc_x+ptl_ns,:,3,:) = zeros(ptl_ns,adc_y,1,n_chnA);
    data = reshape(permute(data,[1,3,2,4]),[(adc_x+ptl_ns)*n_ech,adc_y,n_chnA]);

    %update the number of loops per petal, and the number of ADC points:
    n_ptl_tr        = n_ptl_tr+3;
    hdr.n_ptl_tr    = hdr.n_ptl_tr+3;
    hdr.adc_x       = hdr.adc_x + hdr.ptl_ns;
    adc_x           = adc_x + ptl_ns;

    %Separate the averages from the rosette petals:
    if ltavg
        data=reshape(permute(data,[1,3,2]),[adc_x*n_ech,n_chnA,n_ptl_fov,n_avg]);
    elseif ~ltavg
        data=reshape(permute(data,[1,3,2]),[adc_x*n_ech,n_chnA,n_avg,n_ptl_fov]);
    end

    %Trim off the last point in the time dimension, which are all
    %zeros anyway:
    %FullData=dat(1:end-1,:,:,:);

    %Finally, write a struct to fill in the dims order.  Using the Siemens notation;
    hdr.dims={'Col','Cha','Ave','Lin'};

elseif strcmp(fname(end-2:end),'.h5')
    %This is a .h5 file.  Process accordingly:

    %Load the archive:
    handle = GERecon('Archive.Load',fname);

    %Now get the header information:
    %------------------------------

    %Date and Time
    

    %Store most header information in temporary variable hdrTmp:
    hdrTmp1 = handle.DownloadData.rdb_hdr_rec;
    hdrTmp2 = handle.DownloadData.rdb_hdr_image;
    hdrTmp3 = handle.DownloadData.rdb_hdr_exam;
    hdrTmp4 = handle.DownloadData.rdb_hdr_series;

    hdr.scn_date = hdrTmp1.rdb_hdr_scan_date;   %Scan date
    hdr.scn_time = hdrTmp1.rdb_hdr_scan_time;   %Scan time
    hdr.pt_siz = hdrTmp1.rdb_hdr_frame_size;    %point size
    hdr.adc_x = hdrTmp1.rdb_hdr_rc_xres;        %adc points
    hdr.adc_y = hdrTmp1.rdb_hdr_rc_yres;        %
    hdr.n_slc = hdrTmp1.rdb_hdr_nslices;        %number of slices
    hdr.n_ech = hdrTmp1.rdb_hdr_nechoes;        %number of echoes
    hdr.t_ech = hdrTmp2.te/1000;                %echo time [ms]
    hdr.t_rep = hdrTmp2.tr/1000;                %repetition time [ms]
    hdr.exc_fa = hdrTmp2.mr_flip;               %flip angle [degrees]
    hdr.ltavg = 0;                              %JN - Long term averaging.  Hard coded becuase not sure how to find this param in Orchestra
    hdr.n_avg = hdrTmp2.user36;                 %Number of averages
    hdr.n_ptl_fov = hdrTmp2.user35;             %Number of petal angles acquired (i.e. Nsh) 
    hdr.n_ptl_tr = hdrTmp2.user34;              %Number of petal traversals per TR
    hdr.n_x = hdrTmp2.user33;                   %Matrix size in x direction (Nx)
    hdr.ptl_ns = hdrTmp2.user32;                %Number of samples per petal traversal
    hdr.ptl_dt = hdrTmp2.user31;                %ADC dwell time [us]
    hdr.sar_pk = hdrTmp2.sarpeak;               %Peak SAR [W/kg]
    hdr.sar_avg = hdrTmp2.saravg;               %Average SAR [W/kg]
    hdr.slthk = hdrTmp2.slthick;                %Slice thickness [mm]
    hdr.t_acq = hdrTmp2.sctime/1e6;             %Exam duration [s]
    hdr.fov = hdrTmp2.dfov/10;                  %Field of view [cm]
    hdr.n_chn = 44;                             %JN - Number of channels.  Hard coded because not sure how to find this param in Orchestra.  will find this later
    hdr.coil = hdrTmp2.usedCoilData;            %Name of RF coil
    hdr.ser_dsc = hdrTmp4.se_desc;              %Series description (Same as pulse sequence name)
    hdr.ser_no = handle.SeriesNumber;           %Series number
    hdr.pat_id = hdrTmp3.patidff;               %Patient ID
    hdr.pat_nam = hdrTmp3.patnameff;            %Patient Name
    hdr.psd_nam = hdrTmp2.psdname;              %Pulse sequence name (Same as series description)
    hdr.exm_dsc = hdrTmp3.ex_desc;              %Exam description
    hdr.ex_no = handle.ExamNumber;              %Exam number
    hdr.Bo = hdrTmp3.magstrength/10000;         %Magnetic field strength [T]


    %Now get the raw data:
    %------------------------------
    
    %Load the first "slice":
    control = GERecon('Archive.Next',handle);
    
    %Now loop through the remaining "slices" and keep loading them until there are
    %none left:
    n=1;
    while ~isfield(control,'isAcqDone') || ~isfield(control,'isScanDone')
        data(:,:,n) = control.Data;
        control = GERecon('Archive.Next', handle);
        n=n+1;
    end

    %Correct for RF chopping - negate every 3rd and 4th line along PE
    data(:,:,3:4:end) = -1 * data(:,:,3:4:end);
    data(:,:,4:4:end) = -1 * data(:,:,4:4:end);
    
    %Now organize the data into the following shape: [nx, nptl*nave, nech, nch]
    data=reshape(data,[hdr.adc_x,hdr.n_chn,hdr.n_ech,hdr.n_avg*hdr.n_ptl_fov]);
    data=permute(data,[1,4,3,2]);

    %Concatenate the echoes onto eachother.  But before doing that, add one
    %point to the end of each adc chunk to account for the fact that there is a
    %dead time of one petal rotation (ptl_dt) between the three ADC chunks.
    %The value of that point will be the average of the last point in the
    %previous chunk and the first point in the subsequent chunk:
    data(hdr.adc_x+1:hdr.adc_x+hdr.ptl_ns,:,1,:) = mean(cat(3,data(hdr.adc_x-(hdr.ptl_ns-1):hdr.adc_x,:,1,:),data(1:hdr.ptl_ns,:,2,:)),3);
    data(hdr.adc_x+1:hdr.adc_x+hdr.ptl_ns,:,2,:) = zeros(hdr.ptl_ns,hdr.adc_y,1,hdr.n_chn);
    data = reshape(permute(data,[1,3,2,4]),[(hdr.adc_x+hdr.ptl_ns)*hdr.n_ech,hdr.adc_y,hdr.n_chn]);

    %update the number of loops per petal, and the number of ADC points:
    hdr.n_ptl_tr    = hdr.n_ptl_tr+2;
    hdr.adc_x       = hdr.adc_x + hdr.ptl_ns;

    %Separate the averages from the rosette petals:
    if hdr.ltavg
        data=reshape(permute(data,[1,3,2]),[hdr.adc_x*hdr.n_ech,hdr.n_chn,hdr.n_ptl_fov,hdr.n_avg]);
    elseif ~hdr.ltavg
        data=reshape(permute(data,[1,3,2]),[hdr.adc_x*hdr.n_ech,hdr.n_chn,hdr.n_avg,hdr.n_ptl_fov]);
    end

    %Trim off the last point in the time dimension, which are all
    %zeros anyway:
    %FullData=dat(1:end-1,:,:,:);

    %Finally, write a struct to fill in the dims order.  Using the Siemens notation;
    hdr.dims={'Col','Cha','Ave','Lin'};

end






