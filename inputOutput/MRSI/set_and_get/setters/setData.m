function MRSIStruct = setData(MRSIStruct, data)
    MRSIStruct.data = data;
    MRSIStruct = setSize(MRSIStruct, size(data));
end