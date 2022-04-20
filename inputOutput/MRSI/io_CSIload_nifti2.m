function nifti = MRSI_load_nifti_2(filename)
    arguments
        filename (1, :) char {mustBeFile}
    end
    matlabError = [];
    nifti = struct();
    try
        %open the file
        fileID  = fopen(filename, 'r');

        %sizeof_hdr: Needs to be 540
        header = fread(fileID, 540, '*int8')';

        nifti = readHeader(nifti, header);
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

function nifti = readHeader(nifti, header)
    nifti.sizeof_hdr = typecast(header(1:4), 'int32');
    nifti.magic = cast(header(5:12), 'char');
    nifti.datatype =  getDataTypeFromInt(typecast(header(13:14)', 'int16'));
    nifti.bitpix =  typecast(header(15:16), 'int16');
    nifti.dim =  typecast(header(17:80), 'int64');
    nifti.intent_p1 =  typecast(header(81:88), 'double');
    nifti.intent_p2 =  typecast(header(89 :96 ), 'double');
    nifti.intent_p3 =  typecast(header(97 :104 ), 'double');
    nifti.pixdim = typecast(header(105 :168 ), 'double');
    nifti.vox_offset = typecast(header(169 :176 ), 'int64');
    nifti.scl_slope = typecast(header(177 :184 ), 'double');
    nifti.scl_inter = typecast(header(185 :192 ), 'double');
    nifti.cal_max = typecast(header(193 :200 ), 'double');
    nifti.cal_min = typecast(header(201 :208 ), 'double');
    nifti.slice_duration = typecast(header(209: 216), 'double');
    nifti.tofset = typecast(header(217 :224 ), 'double');
    nifti.slice_start = typecast(header(225 :232 ), 'int64');
    nifti.slice_end = typecast(header(233 :240 ), 'int64');
    nifti.descrip = cast(header(241 :320 ), 'char');
    nifti.aux_file = cast(header(321 :344 ), 'char');
    nifti.qform_cod = typecast(header(345 :348 ), 'int32');
    nifti.sform_cod = typecast(header(349 :352 ), 'int32');
    nifti.quatern_b = typecast(header(353 :360 ), 'double');
    nifti.quatern_c = typecast(header(361 :368 ), 'double');
    nifti.quatern_d = typecast(header(369 :376 ), 'double');
    nifti.qoffset_x = typecast(header(377 :384 ), 'double');
    nifti.qoffset_y = typecast(header(385 :392 ), 'double');
    nifti.qoffset_z = typecast(header(393 :400 ), 'double');
    nifti.srow_x = typecast(header(401 :432 ), 'double');
    nifti.srow_y = typecast(header(433 :464 ), 'double');
    nifti.srow_z = typecast(header(465 :496 ), 'double');
    nifti.slice_code = typecast(header(497 :500 ), 'int32');
    nifti.xyzt_units = typecast(header(501 :504 ), 'int32');
    nifti.intent_cod3 = typecast(header(505 :508 ), 'int32');
    nifti.intent_name = cast(header(509: 524), 'char');
    nifti.dim_info = cast(header(525 : 525), 'char');
end

function nifti = readExtensions(nifti, fileID)
    extension = fread(fileID, 4, 'int8=>int8');
    extensionLength = fread(fileID, 1, 'int32');
    nifti.extensionCode = fread(fileID, 1, 'int32');
    extensionJson = fread(fileID, extensionLength - 8, '*char')';
    nifti.extension = jsondecode(deblank(extensionJson));
end