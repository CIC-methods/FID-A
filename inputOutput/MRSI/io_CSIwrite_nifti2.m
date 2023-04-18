function nifti = io_CSIwrite_nifti2(MRSIStruct, filename)
    arguments
        MRSIStruct (1, 1) struct
        filename (1, :) char
    end
    matlabError = [];
    nifti = struct();
    try
        %open the file
        fileID  = fopen(filename, 'w');

%         %sizeof_hdr: Needs to be 540
%         header = fread(fileID, 540, '*int8')';

        nifti = writeHeader(MRSIStruct, fileID);
        nifti = readExtensions(nifti, fileID);

        if(contains(nifti.datatype, 'complex', 'IgnoreCase', true))
            if(contains(nifti.datatype, '64'))
                dataType = 'float';
            elseif(contains(nifti.datatype, '128'))
                dataType = 'double';
            else
                fprintf('Datatype is %s bits long! Matlab does not contain this datatype. Using chars', bits{1});
                dataType = 'char';
            end
            data = fread(fileID, Inf, dataType);
            complexData = reshape(data, 2, []);
            data = complex(complexData(1, :), complexData(2, :));
        else
            data = fread(fileID, Inf, nifti.datatype);
        end
        data = reshape(data, nifti.dim(2:end));
        nifti.data = data;
    catch matlabError
    end
    fclose(fileID);
    if ~isempty(matlabError)
        rethrow(matlabError);
    end
end


function dataType = getDataTypeFromInt(dataTypeInt)

    intTypes = [2, 4, 8, 16, 32, 64, 128, 256, 512, 768, 1024, 1280, 1536, 1792, 2048];
    dataLabels = ["UINT8", "INT16", "INT32", ...
        "FLOAT32", "COMPLEX64", "FLOAT64", ...
        "RGB24", "INT8", "UINT16", "UINT32", ...
        "INT64", "UINT64", "FLOAT128", ...
        "COMPLEX128", "COMPLEX256"];
    map = containers.Map(intTypes, dataLabels);
    dataType = map(dataTypeInt);
end

function MRSIStruct = writeHeader(MRSIStruct, fileID)
    % header length
    fwrite(fileID, 540, 'int32');
    % magic string
    fwrite(fileID, 'n+2', 'char')
    % data type
    fwrite(fileID, 1792, 'int16') %fwrite(fileID, 1792, 5, 'int16')
    % data size
    fwrite(fileID, 128, 'int16')
    % dimensions
    dimensions = getSizeFromDimensions(MRSIStruct, {'x', 'y', 'z', 't', 'averages', 'coils', 'extras'});
    dimensions = [7, dimensions];
    fwrite(fileID, dimensions, 'int64');
    % intent code
    fwrite(fileID, [0, 0, 0], 'double')
    % pixel dimensions
    pixelDimensions = [1, getVoxSize(MRSIStruct, 'x'), getVoxSize(MRSIStruct, 'y'), ...
                        getVoxSize(MRSIStruct, 'z'), getSpectralWidth(MRSIStruct), 1, 1, 1];

    fwrite(fileID, pixelDimensions, 'double')
    
    MRSIStruct.vox_offset = typecast(header(169 :176 ), 'int64');
    MRSIStruct.scl_slope = typecast(header(177 :184 ), 'double');
    MRSIStruct.scl_inter = typecast(header(185 :192 ), 'double');
    MRSIStruct.cal_max = typecast(header(193 :200 ), 'double');
    MRSIStruct.cal_min = typecast(header(201 :208 ), 'double');
    MRSIStruct.slice_duration = typecast(header(209: 216), 'double');
    MRSIStruct.tofset = typecast(header(217 :224 ), 'double');
    MRSIStruct.slice_start = typecast(header(225 :232 ), 'int64');
    MRSIStruct.slice_end = typecast(header(233 :240 ), 'int64');
    MRSIStruct.descrip = cast(header(241 :320 ), 'char');
    MRSIStruct.aux_file = cast(header(321 :344 ), 'char');
    MRSIStruct.qform_cod = typecast(header(345 :348 ), 'int32');
    MRSIStruct.sform_cod = typecast(header(349 :352 ), 'int32');
    MRSIStruct.quatern_b = typecast(header(353 :360 ), 'double');
    MRSIStruct.quatern_c = typecast(header(361 :368 ), 'double');
    MRSIStruct.quatern_d = typecast(header(369 :376 ), 'double');
    MRSIStruct.qoffset_x = typecast(header(377 :384 ), 'double');
    MRSIStruct.qoffset_y = typecast(header(385 :392 ), 'double');
    MRSIStruct.qoffset_z = typecast(header(393 :400 ), 'double');
    MRSIStruct.srow_x = typecast(header(401 :432 ), 'double');
    MRSIStruct.srow_y = typecast(header(433 :464 ), 'double');
    MRSIStruct.srow_z = typecast(header(465 :496 ), 'double');
    MRSIStruct.slice_code = typecast(header(497 :500 ), 'int32');
    MRSIStruct.xyzt_units = typecast(header(501 :504 ), 'int32');
    MRSIStruct.intent_cod3 = typecast(header(505 :508 ), 'int32');
    MRSIStruct.intent_name = cast(header(509: 524), 'char');
    MRSIStruct.dim_info = cast(header(525 : 525), 'char');
end

function nifti = readExtensions(nifti, fileID)
    extension = fread(fileID, 4, 'int8=>int8');
    extensionLength = fread(fileID, 1, 'int32');
    nifti.extensionCode = fread(fileID, 1, 'int32');
    extensionJson = fread(fileID, extensionLength - 8, '*char')';
    nifti.extension = jsondecode(deblank(extensionJson));
end