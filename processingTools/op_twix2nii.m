function mrs_struct=op_twix2nii(mrs_struct)
%Parsing FID-A format to NIfTI-MRS format as described in:
%`Clarke WT, Bell TK, Emir UE, Mikkelsen M, Oeltzschner G, Shamaei A, Soher
%BJ, Wilson M. NIfTI-MRS: A standard data format for magnetic resonance
%spectroscopy. Magn Reson Med. 2022. doi: 10.1002/mrm.29418.` & https://www.nitrc.org/forum/message.php?msg_id=3738


% C Header Struct
% /*! \struct nifti_2_header
%     \brief Data structure defining the fields in the nifti2 header.
%            This binary header should be found at the beginning of a valid
%            NIFTI-2 header file.
%  */
%                         /*************************/  /************************/ /************/
% struct nifti_2_header { /* NIFTI-2 usage         */  /* NIFTI-1 usage        */ /*  offset  */
%                         /*************************/  /************************/ /************/
% int   sizeof_hdr;     /*!< MUST be 540           */  /* int sizeof_hdr; (348) */  /*   0 */
% char  magic[8] ;      /*!< MUST be valid signature. */  /* char magic[4];     */  /*   4 */
% short datatype;       /*!< Defines data type!    */  /* short datatype;       */  /*  12 */
% short bitpix;         /*!< Number bits/voxel.    */  /* short bitpix;         */  /*  14 */
% int64_t dim[8];       /*!< Data array dimensions.*/  /* short dim[8];         */  /*  16 */
% double intent_p1 ;    /*!< 1st intent parameter. */  /* float intent_p1;      */  /*  80 */
% double intent_p2 ;    /*!< 2nd intent parameter. */  /* float intent_p2;      */  /*  88 */
% double intent_p3 ;    /*!< 3rd intent parameter. */  /* float intent_p3;      */  /*  96 */
% double pixdim[8];     /*!< Grid spacings.        */  /* float pixdim[8];      */  /* 104 */
% int64_t vox_offset;   /*!< Offset into .nii file */  /* float vox_offset;     */  /* 168 */
% double scl_slope ;    /*!< Data scaling: slope.  */  /* float scl_slope;      */  /* 176 */
% double scl_inter ;    /*!< Data scaling: offset. */  /* float scl_inter;      */  /* 184 */
% double cal_max;       /*!< Max display intensity */  /* float cal_max;        */  /* 192 */
% double cal_min;       /*!< Min display intensity */  /* float cal_min;        */  /* 200 */
% double slice_duration;/*!< Time for 1 slice.     */  /* float slice_duration; */  /* 208 */
% double toffset;       /*!< Time axis shift.      */  /* float toffset;        */  /* 216 */
% int64_t slice_start;  /*!< First slice index.    */  /* short slice_start;    */  /* 224 */
% int64_t slice_end;    /*!< Last slice index.     */  /* short slice_end;      */  /* 232 */
% char  descrip[80];    /*!< any text you like.    */  /* char descrip[80];     */  /* 240 */
% char  aux_file[24];   /*!< auxiliary filename.   */  /* char aux_file[24];    */  /* 320 */
% int qform_code ;      /*!< NIFTI_XFORM_* code.   */ /* short qform_code;      */  /* 344 */
% int sform_code ;      /*!< NIFTI_XFORM_* code.   */ /* short sform_code;      */  /* 348 */
% double quatern_b ;    /*!< Quaternion b param.   */ /* float quatern_b;       */  /* 352 */
% double quatern_c ;    /*!< Quaternion c param.   */ /* float quatern_c;       */  /* 360 */
% double quatern_d ;    /*!< Quaternion d param.   */ /* float quatern_d;       */  /* 368 */
% double qoffset_x ;    /*!< Quaternion x shift.   */ /* float qoffset_x;       */  /* 376 */
% double qoffset_y ;    /*!< Quaternion y shift.   */ /* float qoffset_y;       */  /* 384 */
% double qoffset_z ;    /*!< Quaternion z shift.   */ /* float qoffset_z;       */  /* 392 */
% double srow_x[4] ;    /*!< 1st row affine transform. */  /* float srow_x[4];  */  /* 400 */
% double srow_y[4] ;    /*!< 2nd row affine transform. */  /* float srow_y[4];  */  /* 432 */
% double srow_z[4] ;    /*!< 3rd row affine transform. */  /* float srow_z[4];  */  /* 464 */
% int slice_code ;      /*!< Slice timing order.   */  /* char slice_code;      */  /* 496 */
% int xyzt_units ;      /*!< Units of pixdim[1..4] */  /* char xyzt_units;      */  /* 500 */
% int intent_code ;     /*!< NIFTI_INTENT_* code.  */  /* short intent_code;    */  /* 504 */
% char intent_name[16]; /*!< 'name' or meaning of data. */ /* char intent_name[16]; */  /* 508 */
% char dim_info;        /*!< MRI slice ordering.   */      /* char dim_info;        */  /* 524 */
% char unused_str[15];  /*!< unused, filled with \0 */                                  /* 525 */
% } ;                   /**** 540 bytes total ****/
% typedef struct nifti_2_header nifti_2_header ;
 
%         sizeof_hdr: 540
%              magic: 'n+2 ←↵'
%           datatype: 32
%             bitpix: 64
%                dim: [6 1 1 1 4096 16 128 1]
%          intent_p1: 0
%          intent_p2: 0
%          intent_p3: 0
%             pixdim: [1 30 40 30 2.5000e-04 1 1 1]
%         vox_offset: 1680
%          scl_slope: 1
%          scl_inter: 0
%            cal_max: 0
%            cal_min: 0
%     slice_duration: 0
%            toffset: 0
%        slice_start: 0
%          slice_end: 0
%            descrip: ''
%           aux_file: ''
%         qform_code: 0
%         sform_code: 2
%          quatern_b: 1
%          quatern_c: 0
%          quatern_d: 0
%          qoffset_x: -1.1622
%          qoffset_y: -11.5254
%          qoffset_z: 16.9492
%             srow_x: [30 0 0 -1.1622]
%             srow_y: [0 -40 0 -11.5254]
%             srow_z: [0 0 -30 16.9492]
%         slice_code: 0
%         xyzt_units: 0
%        intent_code: 0
%        intent_name: 'mrs_v0_5'
%           dim_info: 0
%         unused_str: ''
%          extension: [1 0 0 0]
%            version: 2
%        swap_endian: 0
%          file_name: '/Users/competer/Documents/data_MR_scans/2022-04-20_JNEA1U_HERMESINVIVO/MEGA-GABA/peter_test_spec2nii.nii.gz'

          
nii_mrs.hdr.sizeof_hdr=540;
nii_mrs.hdr.magic=magic(8);
nii_mrs.hdr.datatype=32;
nii_mrs.hdr.bitpix=64;
%nii_mrs.hdr.dim=[? vox_x vox_y vox_z num_pts coils averages ?]
nii_mrs.hdr.intent_p1=0;
nii_mrs.hdr.intent_p2=0;
nii_mrs.hdr.intent_p3=0;
%nii_mrs.hdr.pixdim=[? sz_x sz_y sz_z dwell_time ? ? ?]
nii_mrs.hdr.vox_offset=1680;
nii_mrs.hdr.scl_slope=1; 
nii_mrs.hdr.scl_inter=0;
nii_mrs.hdr.cal_max=0;
nii_mrs.hdr.cal_min=0;
nii_mrs.hdr.slice_duration=0;
nii_mrs.hdr.toffset=0;
nii_mrs.hdr.slice_start=0;
nii_mrs.hdr.slice_end=0;
nii_mrs.hdr.descrip='';
nii_mrs.hdr.aux_file='';
nii_mrs.hdr.qform_code=0;
nii_mrs.hdr.sform_code=2;
nii_mrs.hdr.quatern_b=1;
nii_mrs.hdr.quatern_c=0;
nii_mrs.hdr.quatern_d=0;
nii_mrs.hdr.qoffset_x=-mrs_struct.hdr.MeasYaps.sSpecPara.sVoI.sPosition.dSag;
nii_mrs.hdr.qoffset_y=-mrs_struct.hdr.MeasYaps.sSpecPara.sVoI.sPosition.dCor;
nii_mrs.hdr.qoffset_z=mrs_struct.hdr.MeasYaps.sSpecPara.sVoI.sPosition.dTra;
nii_mrs.hdr.srow_x=[mrs_struct.hdr.MeasYaps.sSpecPara.sVoI.dThickness 0 0 nii_mrs.hdr.qoffset_x];
nii_mrs.hdr.srow_y=[0 mrs_struct.hdr.MeasYaps.sSpecPara.sVoI.dPhaseFOV 0 nii_mrs.hdr.qoffset_y];
nii_mrs.hdr.srow_z=[0 0 mrs_struct.hdr.MeasYaps.sSpecPara.sVoI.dReadoutFOV nii_mrs.hdr.qoffset_z];
nii_mrs.hdr.slice_code=0;
nii_mrs.hdr.xyzt_units=0;
nii_mrs.hdr.intent_code=0;
nii_mrs.hdr.intent_name='mrs_v0_5';
nii_mrs.hdr.dim_info=0;
nii_mrs.hdr.unused_str='';
nii_mrs.hdr.extension=[1 0 0 0];
nii_mrs.hdr.version=2;
nii_mrs.hdr.swap_endian=0;
nii_mrs.hdr.filename=file_name;


nii_mrs.hdr_ext.SpectrometerFrequency=mrs_struct.hdr.Meas.lFrequency/1e6;
nii_mrs.hdr_ext.ResonantNucleus=mrs_struct.hdr.Meas.Nucleus;
nii_mrs.hdr_ext.dim5='DIM_COIL';
nii_mrs.hdr_ext.dim6='DIM_DYN';
nii_mrs.hdr_ext.EchoTime=mrs_struct.hdr.Meas.TE/1e6;
nii_mrs.hdr_ext.RepetitionTime=mrs_struct.hdr.Meas.TR/1e6;
nii_mrs.hdr_ext.ExcitationFlipAngle=mrs_struct.hdr.Meas.FlipAngle;
nii_mrs.hdr_ext.TxOffset=mrs_struct.hdr.Meas.dDeltaFrequency;
nii_mrs.hdr_ext.Manufacturer=mrs_struct.hdr.Meas.Manufacturer;
nii_mrs.hdr_ext.ManufacturersModelName=mrs_struct.hdr.Meas.ManufacturersModelName;
nii_mrs.hdr_ext.DeviceSerialNumber=mrs_struct.hdr.Meas.DeviceSerialNumber;
nii_mrs.hdr_ext.SoftwareVersions=mrs_struct.hdr.Meas.SoftwareVersions;
nii_mrs.hdr_ext.InstitutionName=mrs_struct.hdr.Meas.InstitutionName;
nii_mrs.hdr_ext.InstitutionAddress=mrs_struct.hdr.Meas.InstitutionAddress;
nii_mrs.hdr_ext.RxCoil=mrs_struct.hdr.MeasYaps.sCoilSelectMeas.aRxCoilSelectData{1}.asList{1}.sCoilElementID.tCoilID
nii_mrs.hdr_ext.SequenceName=mrs_struct.hdr.Meas.SequenceString;
nii_mrs.hdr_ext.ProtocolName=mrs_struct.hdr.Meas.ProtocolName;
nii_mrs.hdr_ext.PatientPosition=mrs_struct.hdr.Meas.PatientPosition;
nii_mrs.hdr_ext.PatientName='';%mrs_struct.hdr.Meas.PatientName;
nii_mrs.hdr_ext.PatientWeight='';%mrs_struct.hdr.Meas.PatientWeight;
nii_mrs.hdr_ext.PatientDoB='';%mrs_struct.hdr.Meas.PatientBirthday;
nii_mrs.hdr_ext.PatientSex='';%mrs_struct.hdr.Meas.PatientSex;
nii_mrs.hdr_ext.ConversionMethod='FID-A spec2nii v0';
nii_mrs.hdr_ext.ConvertionTime=datetime(now,'ConvertFrom','datenum');
nii_mrs.hdr_ext.OriginalFile=mrs_struct.filename;
nii_mrs.hdr_ext.kSpace=[0,0,0];
nii_mrs.hdr_ext.PulseSequenceFile.Value=mrs_struct.hdr.Meas.tSequenceFileName;
nii_mrs.hdr_ext.PulseSequenceFile.Description='Sequence binary path.';
nii_mrs.hdr_ext.IceProgramFile.Value=mrs_struct.hdr.Meas.tICEProgramName;
nii_mrs.hdr_ext.IceProgramFile.Description='Reconstruction binary path.';

%      SpectrometerFrequency: 123.2500
%            ResonantNucleus: {'1H'}
%                      dim_5: 'DIM_COIL'
%                      dim_6: 'DIM_DYN'
%                   EchoTime: 0.0680
%             RepetitionTime: 3
%        ExcitationFlipAngle: 90
%                   TxOffset: -1.7000
%               Manufacturer: 'SIEMENS'
%     ManufacturersModelName: 'Prisma'
%         DeviceSerialNumber: '66049.0'
%           SoftwareVersions: 'syngo MR E11'
%            InstitutionName: 'Sunnybrook Research Institute'
%         InstitutionAddress: 'Bayview Avenue 2075,Toronto,Ontario,CA,M4N 3M5'
%                     RxCoil: '"HeadNeck_20"'
%               SequenceName: 'svs_se'
%               ProtocolName: 'megapress_GABA'
%            PatientPosition: 'HFS'
%                PatientName: 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
%              PatientWeight: 79.3787
%                 PatientDoB: 'xxxxxxxx'
%                 PatientSex: 'F'
%           ConversionMethod: 'spec2nii v0.4.9'
%             ConversionTime: '2023-01-10T14:21:06.251'
%               OriginalFile: {'meas_MID01209_FID16258_megapress_GABA.dat'}
%                     kSpace: [3×1 logical]
%          PulseSequenceFile: [1×1 struct]
%             IceProgramFile: [1×1 struct]

mrs_struct.nii_mrs=nii_mrs;


           
mrs_struct.vox_offset = typecast(header(169 :176 ), 'int64');
mrs_struct.scl_slope = typecast(header(177 :184 ), 'double');
mrs_struct.scl_inter = typecast(header(185 :192 ), 'double');
mrs_struct.cal_max = typecast(header(193 :200 ), 'double');
mrs_struct.cal_min = typecast(header(201 :208 ), 'double');
mrs_struct.slice_duration = typecast(header(209: 216), 'double');
mrs_struct.tofset = typecast(header(217 :224 ), 'double');
mrs_struct.slice_start = typecast(header(225 :232 ), 'int64');
mrs_struct.slice_end = typecast(header(233 :240 ), 'int64');
mrs_struct.descrip = cast(header(241 :320 ), 'char');
mrs_struct.aux_file = cast(header(321 :344 ), 'char');
mrs_struct.qform_cod = typecast(header(345 :348 ), 'int32');
mrs_struct.sform_cod = typecast(header(349 :352 ), 'int32');
mrs_struct.quatern_b = typecast(header(353 :360 ), 'double');
mrs_struct.quatern_c = typecast(header(361 :368 ), 'double');
mrs_struct.quatern_d = typecast(header(369 :376 ), 'double');
mrs_struct.qoffset_x = typecast(header(377 :384 ), 'double');
mrs_struct.qoffset_y = typecast(header(385 :392 ), 'double');
mrs_struct.qoffset_z = typecast(header(393 :400 ), 'double');
mrs_struct.srow_x = typecast(header(401 :432 ), 'double');
mrs_struct.srow_y = typecast(header(433 :464 ), 'double');
mrs_struct.srow_z = typecast(header(465 :496 ), 'double');
mrs_struct.slice_code = typecast(header(497 :500 ), 'int32');
mrs_struct.xyzt_units = typecast(header(501 :504 ), 'int32');
mrs_struct.intent_cod3 = typecast(header(505 :508 ), 'int32');
mrs_struct.intent_name = cast(header(509: 524), 'char');
mrs_struct.dim_info = cast(header(525 : 525), 'char');
end