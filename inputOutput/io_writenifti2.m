function nifti=io_writenifti2(mrs_struct,filename)
arguments
    mrs_struct (1, 1) struct
    filename (1, :) char
end
matlabError=[];
nifti=struct();

if isfield(mrs_struct,'nii_mrs')
    printf('Data has not been parsed for nii_mrs format, please run "op_spec2nii()" first\n');
    return;
end

try
    %open the file
    header=mrs_struct.nii_mrs.hdr;
    header_ext=mrs_struct.nii_mrs.hdr_ext;
    
    fileID=fopen(filename,'w');
    fwrite(fileID, 540, 'int32');
    fwrite(fileID, header.magic, 'char');
    fwrite(fileID, 32, 'int16');
    fwrite(fileID, 128, 'int16');
    fwrite(fileID,header.dim, 'int64');
    fwrite(fileID,header.intent_p1, 'double');
    fwrite(fileID,header.intent_p2, 'double');
    fwrite(fileID,header.intent_p3, 'double');
    fwrite(fileID,header.pixdim, 'double');
    fwrite(fileID,header.vox_offset,'int64');
    fwrite(fileID,header.scl_slope,'double');
    fwrite(fileID,header.scl_inter,'double');
    fwrite(fileID,header.cal_max,'double');
    fwrite(fileID,header.cal_min,'double');
    fwrite(fileID,header.slice_duration,'double');
    fwrite(fileID,header.tofset,'double');
    fwrite(fileID,header.slice_start,'int64');
    fwrite(fileID,header.slice_end,'int64');
    fwrite(fileID,header.descrip,'char');
    fwrite(fileID,header.aux_file,'char');
    fwrite(fileID,header.qform_cod,'int32');
    fwrite(fileID,header.sform_cod,'int32');
    fwrite(fileID,header.quatern_b,'double');
    fwrite(fileID,header.quatern_c,'double');
    fwrite(fileID,header.quatern_d,'double');
    fwrite(fileID,header.quatern_x,'double');
    fwrite(fileID,header.quatern_y,'double');
    fwrite(fileID,header.quatern_z,'double');
    fwrite(fileID,header.srow_x,'double');
    fwrite(fileID,header.srow_y,'double');
    fwrite(fileID,header.srow_z,'double');
    fwrite(fileID,header.slice_code,'int32');
    fwrite(fileID,header.xyzt_units,'int32');
    fwrite(fileID,header.intent_code,'int32');
    fwrite(fileID,header.intent_name,'char');
    fwrite(fileID,header.dim_info,'char');
    fwrite(fileID,header.unused_str,'char');
    fwrite(fileID,header.extension,'uint8');
    
    ecode=44;
    edata=uint8(jsonencode(header_ext))';
    esize=size(edata,2)+ecode;
    fwrite(fileID,esize,'int32');
    fwrite(fileID,ecode,'int32');
    fwrite(fileID,edata,'uint8');
    
%     
%     
%     if(contains(nifti.datatype, 'complex', 'IgnoreCase', true))
%         if(contains(nifti.datatype, '64'))
%             dataType = 'float';
%         elseif(contains(nifti.datatype, '128'))
%             dataType = 'double';
%         else
%             fprintf('Datatype is %s bits long! Matlab does not contain this datatype. Using chars', bits{1});
%             dataType = 'char';
%         end
%         data = fread(fileID, Inf, dataType);
%         complexData = reshape(data, 2, []);
%         data = complex(complexData(1, :), complexData(2, :));
%     else
%         data = fread(fileID, Inf, nifti.datatype);
%     end
    data = reshape(data, nifti.dim(2:end));
    nifti.data = data;
catch matlabError
end
fclose(fileID);
if ~isempty(matlabError)
    rethrow(matlabError);
end
end


% function nifti = readExtensions(nifti, fileID)
% extension = fread(fileID, 4, 'int8=>int8');
% extensionLength = fread(fileID, 1, 'int32');
% nifti.extensionCode = fread(fileID, 1, 'int32');
% extensionJson = fread(fileID, extensionLength - 8, '*char')';
% nifti.extension = jsondecode(deblank(extensionJson));
% end