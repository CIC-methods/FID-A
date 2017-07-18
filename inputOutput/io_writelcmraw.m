%io_writelcmraw.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% RF=io_writelcmraw(data_struct,outfile,metab);
% 
% DESCRIPTION:
% Take a simulated metabolite basis spectrum in matlab structure format, and
% output it into LCModel RAW format to be used in an LCModel basis spectrum.
% 
% INPUTS:
% data_struct    = simulated metabolite basis spectrum in matlab structure
%                   format.
% outfile        = name of the output .RAW file.
% metab          = Abbreviated name of the metabolite (ie. 'Cr', 'Glu',
%                   etc.)
%
% OUTPUTS:
% RF            = Same as input.  Not used.  The primary output of this
%                   function is a text file in LCModel raw format.  
    
function RF=io_writelcmraw(data_struct,outfile,metab);
%function RF=writelcmraw(data_struct,outfile);

datsets=1;
zop=0;
t0=0;
Bo=2.9;
Nuc=0;
PatName='No Name';
scanner='TrioTim';
addinfo='jnear';
expname=[data_struct.seq num2str(data_struct.te) '_' data_struct.sim 'Sim'];
investigator='jnear';
comment=['Simulated LCModel basis set for ' data_struct.sim ' ' data_struct.seq...
    ' experiment with echo time ' num2str(data_struct.te) ' ms'];

%Take the complex conjugate becuase the sense of rotation in LCModel seems to
%be opposite to that used in FID-A.
data_struct=op_complexConj(data_struct);

RF=zeros(length(data_struct.fids),2);
RF(:,1)=real(data_struct.fids);
RF(:,2)=imag(data_struct.fids);



%write to txt file for jmrui
fid=fopen(outfile,'w+');
fprintf(fid,' This RAW file was created using the MATLAB spectral');
fprintf(fid,'\n simulation program developed by Robin Simpson');
fprintf(fid,'\n and Jamie Near, FMRIB 2010');
fprintf(fid,'\n\n Experiment Name : %s' ,expname);
fprintf(fid,'\n Investigator    : %s',investigator);
fprintf(fid,'\n Comment : %s',comment);
fprintf(fid,'\n\n User defined parameters:');
fprintf(fid,'\n\n Sweep Width = %1.4f Hz',data_struct.spectralwidth);
fprintf(fid,'\n Vector Size = %i points',data_struct.n);
fprintf(fid,'\n Apodization = %1.4f Hz',data_struct.linewidth);
fprintf(fid,'\n B0 Field    = %1.4f T',data_struct.Bo);
fprintf(fid,'\n\n $NMID ID=%s',metab);
fprintf(fid,'\n FMTDAT=(2E16.6)');
fprintf(fid,'\n VOLUME=1.00000');
fprintf(fid,'\n TRAMP=1.00000');
fprintf(fid,'\n $END');
fprintf(fid,'\n  % 1.6E  % 1.6E',RF'); %space after "%" sign forces hanging negative sign
fprintf(fid,'\n');
fclose(fid);
