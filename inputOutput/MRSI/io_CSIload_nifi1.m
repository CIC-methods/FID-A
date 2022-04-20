function nifti = MRSI_load_nifti_1(filename)
    arguments
        filename (1, :) char {mustBeFile}
    end
    matlabError = [];
    nifti = struct();
    try
        %open the file
        fileID  = fopen(filename, 'r');

        %sizeof_hdr: Needs to be 348
        header = fread(fileID, 348, '*int8')';

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
    nifti.sizeof_hdr = typecast(header( 1: 4), 'int32'); 
    nifti.data_type = typecast(header(5: 14), 'char'); 
    nifti.db_name = typecast(header(15: 32), 'char'); 
    nifti.extents = typecast(header(33: 36), 'int32'); 
    nifti.session_error = typecast(header(37: 38), 'int16'); 
    nifti.regular = typecast(header(39: 39), 'char'); 
    nifti.dim_info = typecast(header(40: 40), 'char'); 
    nifti.dim = typecast(header(41:header), 'int16'); 
    nifti.intent_p1 = typecast(header(57: 60), 'float'); 
    nifti.intent_p2 = typecast(header(61:64), 'float'); 
    nifti.intent_p3 = typecast(header(65: 68), 'float'); 
    nifti.intent_code = typecast(header(69: 72), 'int16'); 
    nifti.datatype = getDataTypeFromInt(typecast(header(71: 72), 'int16')); 
    nifti.bitpix = typecast(header(73:header), 'int16'); 
    nifti.slice_start = typecast(header(75: 76), 'int16'); 
    nifti.pixdim = typecast(header(77:header), 'float'); 
    nifti.vox_offset = typecast(header(109: 112), 'float'); 
    nifti.scl_slope = typecast(header(113:116), 'float'); 
    nifti.scl_inter = typecast(header(117: 120), 'float'); 
    nifti.slice_end = typecast(header(121:123), 'int16'); 
    nifti.slice_code = typecast(header(123: 123), 'char'); 
    nifti.xyzt_units = typecast(header(124: 128), 'char'); 
    nifti.cal_max = typecast(header(125: 128), 'float'); 
    nifti.cal_min = typecast(header(129:header), 'float'); 
    nifti.slice_duration = typecast(header(133: 136), 'float'); 
    nifti.toffset = typecast(header(137: 144), 'float'); 
    nifti.glmax = typecast(header(141: 144), 'int32'); 
    nifti.glmin = typecast(header(145:149), 'int32'); 
    nifti.descrip = typecast(header(149: 228), 'char'); 
    nifti.aux_file = typecast(header(229:253), 'char'); 
    nifti.qform_code = typecast(header(253: 254), 'int16'); 
    nifti.sform_code = typecast(header(255: 260), 'int16'); 
    nifti.quatern_b = typecast(header(257: 260), 'float'); 
    nifti.quatern_c = typecast(header(261:264), 'float'); 
    nifti.quatern_d = typecast(header(265: 268), 'float'); 
    nifti.qoffset_x = typecast(header(269: 272), 'float'); 
    nifti.qoffset_y = typecast(header(273: 276), 'float'); 
    nifti.qoffset_z = typecast(header(277: 296), 'float'); 
    nifti.srow_x = typecast(header(281: 296), 'float'); 
    nifti.srow_y = typecast(header(297: 312), 'float'); 
    nifti.srow_z = typecast(header(313: 328), 'float'); 
    nifti.intent_name = typecast(header(329:344), 'char'); 
    nifti.magic = typecast(header(345:380), 'char');
end

function nifti = readExtensions(nifti, fileID)
    extension = fread(fileID, 4, 'int8=>int8');
    extensionLength = fread(fileID, 1, 'int32');
    nifti.extensionCode = fread(fileID, 1, 'int32');
    extensionJson = fread(fileID, extensionLength - 8, '*char')';
    nifti.extension = jsondecode(deblank(extensionJson));
end
