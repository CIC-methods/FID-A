function [svs_av,svs_cs,svs_lb,svs_cc,svs_raw]=peter_analysis_svs(in)

lw=0;

if nargin<1
    [x,y]=uigetfile('*.dat');
    svs_raw=io_loadspec_twix([y x]);
else
    svs_raw=in;
end
svs_cc=op_combineRcvrs(svs_raw,svs_raw);
svs_lb=op_filter(svs_cc,lw);
svs_aa=op_alignISIS(op_alignAverages(op_alignISIS(op_alignAverages(svs_lb,0.4),0.4),0.4),0.4);
svs_cs=op_combinesubspecs(svs_aa,'diff');
svs_av=op_averaging(svs_cs);

op_plotspec(svs_av,-30,30);


% csi_zf=op_CSIZeroFillTime(csi_apd,4096);
% csi_fft=op_CSIFourierTransform(csi_zf);
% csi_fft.flags.coilCombined=1;
% csi_fft.flags.addedrcvrs=1;
% % csi_b0=op_CSIB0Correction(csi_fft);
% csi_av=op_CSIAverage(csi_fft);


