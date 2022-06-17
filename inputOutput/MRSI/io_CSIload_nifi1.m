function nifti = io_CSIload_nifi1(filename)
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
                dataType = 'single';
            elseif(contains(nifti.datatype, '128'))
                dataType = 'double';
            else
                fprintf('Datatype is %s bits long! Matlab does not contain this datatype. Using chars', nifti.datatype{1});
                dataType = 'char';
            end
            data = fread(fileID, Inf, dataType);
            complexData = reshape(data, 2, []);
            data = complex(complexData(1, :), complexData(2, :));
        else
            try 
                data = fread(fileID, Inf, nifti.datatype);
            catch matlabException
                if(isequal(matlabException.identifier, 'MATLAB:badprecision_mx'))
                    fprintf('Datatype %s is not identified in MATLAB. Reading in chars instead', nifti.datatype);
                    data = fread(fileID, Inf, nifti.datatype);
                else
                    rethrow(matlabException)
                end
            end
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
    dataLabels = ["unit8", "int16", "int32", ...
        "single", "complex64", "double", ...
        "RGB24", "int8", "uint16", "uint32", ...
        "int64", "uint64", "float128", ...
        "complex128", "complex256"];
    map = containers.Map(intTypes, dataLabels);
    dataType = map(dataTypeInt);
end

function nifti = readHeader(nifti, header)
    nifti.sizeof_hdr = typecast(header( 1: 4), 'int32');
    nifti.dim_info = cast(header(40: 40), 'char');
    nifti.dim = typecast(header(41:56), 'int16');
    nifti.intent_p1 = typecast(header(57: 60), 'single');
    nifti.intent_p2 = typecast(header(61:64), 'single');
    nifti.intent_p3 = typecast(header(65: 68), 'single');
    nifti.intent_code = typecast(header(69: 72), 'int16');
    nifti.datatype = getDataTypeFromInt(typecast(header(71: 72), 'int16'));
    nifti.bitpix = typecast(header(73:74), 'int16');
    nifti.slice_start = typecast(header(75: 76), 'int16');
    nifti.pixdim = typecast(header(77:108), 'single');
    nifti.vox_offset = typecast(header(109: 112), 'single');
    nifti.scl_slope = typecast(header(113:116), 'single');
    nifti.scl_inter = typecast(header(117: 120), 'single');
    nifti.slice_end = typecast(header(121:122), 'int16');
    nifti.slice_code = cast(header(123: 123), 'char');
    nifti.xyzt_units = dec2bin(typecast(header(124: 124), 'int8'));
    nifti.cal_max = typecast(header(125: 128), 'single');
    nifti.cal_min = typecast(header(129:132), 'single');
    nifti.slice_duration = typecast(header(133: 136), 'single');
    nifti.toffset = typecast(header(137: 144), 'single');
    nifti.glmax = typecast(header(141: 144), 'int32');
    nifti.glmin = typecast(header(145:148), 'int32');
    nifti.descrip = cast(header(149: 228), 'char');
    nifti.aux_file = cast(header(229:252), 'char');
    nifti.qform_code = typecast(header(253: 254), 'int16');
    nifti.sform_code = typecast(header(255: 256), 'int16');
    nifti.quatern_b = typecast(header(257: 260), 'single');
    nifti.quatern_c = typecast(header(261: 264), 'single');
    nifti.quatern_d = typecast(header(265: 268), 'single');
    nifti.qoffset_x = typecast(header(269: 272), 'single');
    nifti.qoffset_y = typecast(header(273: 276), 'single');
    nifti.qoffset_z = typecast(header(277: 280), 'single');
    nifti.srow_x = typecast(header(281: 296), 'single');
    nifti.srow_y = typecast(header(297: 312), 'single');
    nifti.srow_z = typecast(header(313: 328), 'single');
    nifti.intent_name = cast(header(329:344), 'char');
    nifti.magic = cast(header(345:348), 'char');
end

function nifti = readExtensions(nifti, fileID)
    extension = fread(fileID, 4, 'int8=>int8');
    if(extension(1) ~= 0)
        disp('Extension present. Reading...')
        extensionLength = fread(fileID, 1, 'int32');
        nifti.extensionCode = fread(fileID, 1, 'int32');
        extensionJson = fread(fileID, extensionLength - 8, '*char')';
        nifti.extension = jsondecode(deblank(extensionJson));
    end
end
