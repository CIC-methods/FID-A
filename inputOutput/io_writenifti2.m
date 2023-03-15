function nifti=io_writenifti2(in,filename)
arguments
    in (1, 1) struct
    filename (1, :) char
end
matlabError=[];
nifti=struct();

if ~isfield(in,'nii_mrs')
    printf('Data has not been parsed for nii_mrs format, please run "op_spec2nii()" first\n');
    return;
end

try
    %open the file
    header=in.nii_mrs.hdr;
    header_ext=in.nii_mrs.hdr_ext;
    
    fileID=fopen([filename '.nii'],'w');
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
    fwrite(fileID,header.toffset,'double');
    fwrite(fileID,header.slice_start,'int64');
    fwrite(fileID,header.slice_end,'int64');
    
    descrip=zeros(80,1);
    if numel(header.descrip) > 80 %if the description is larger than 80 elements
        descrip=header.descrip(1:80);
    elseif numel(header.descrip) < 80 %if the description is smaller than 80 elements
        descrip(numel(header.descrip)+1:end)=0;
    end
    fwrite(fileID,descrip,'char');
    
    aux_file=zeros(24,1);
    if numel(header.aux_file) > 24 %if the description is larger than 80 elements
        aux_file=header.aux_file(1:24);
    elseif numel(header.aux_file) < 24 %if the description is smaller than 80 elements
        aux_file(numel(header.aux_file)+1:end)=0;
    end
    fwrite(fileID,aux_file,'char');
    
    fwrite(fileID,header.qform_code,'int32');
    fwrite(fileID,header.sform_code,'int32');
    fwrite(fileID,header.quatern_b,'double');
    fwrite(fileID,header.quatern_c,'double');
    fwrite(fileID,header.quatern_d,'double');
    fwrite(fileID,header.qoffset_x,'double');
    fwrite(fileID,header.qoffset_y,'double');
    fwrite(fileID,header.qoffset_z,'double');
    fwrite(fileID,header.srow_x,'double');
    fwrite(fileID,header.srow_y,'double');
    fwrite(fileID,header.srow_z,'double');
    fwrite(fileID,header.slice_code,'int32');
    fwrite(fileID,header.xyzt_units,'int32');
    fwrite(fileID,header.intent_code,'int32');
    
    intent_name=zeros(16,1);
    if numel(header.intent_name) > 16 %if the description is larger than 80 elements
        intent_name=header.intent_name(1:16);
    elseif numel(header.intent_name) < 16 %if the description is smaller than 80 elements
        intent_name(numel(header.intent_name)+1:end)=0;
    end
    fwrite(fileID,intent_name,'char');
    
    fwrite(fileID,header.dim_info,'uint8');
    
    unused_str=zeros(15,1);
    if numel(header.unused_str) > 15 %if the description is larger than 80 elements
        unused_str=header.unused_str(1:15);
    elseif numel(header.unused_str) < 15 %if the description is smaller than 80 elements
        unused_str(numel(header.unused_str)+1:end)=0;
    end
    fwrite(fileID,unused_str,'char');
    
    fwrite(fileID,header.extension,'uint8');
    
    ecode=44;
    edata=uint8(jsonencode(header_ext))';
    esize= header.vox_offset - 540 - 4;
    fwrite(fileID,esize,'int32');
    fwrite(fileID,ecode,'int32');
    fwrite(fileID,edata,'uint8');
    
    if in.dims.kx||in.dims.ky||in.dims.kz
        img_tmp=reshape(in.data,1,[]);
    else
        img_tmp=reshape(in.fids,1,[]);
    end
    img=zeros([size(img_tmp,2)*2 1]);
    img(1:2:end)=real(img_tmp);
    img(2:2:end)=imag(img_tmp);
    
    n = header.vox_offset - ftell(fileID);
    if n<0 % seen n=-1 for unknown reason
        fseek(fileID, n, 'cof');
    elseif n>0
        fwrite(fileID, zeros(n,1), 'uint8');
    end
    
    fwrite(fileID, img, 'single');
    fclose(fileID); % all written
    gzip([filename '.nii']);
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
%     data = reshape(data, nifti.dim(2:end));
%     nifti.data = data;
catch matlabError
end
% fclose(fileID);
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