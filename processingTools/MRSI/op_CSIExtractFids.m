function [fids, positions] = op_CSIExtractFids(MRSIStruct, kSpaceFile)
    table = readTable(kSpaceFile);
    TR = table.TR;
    Kx = table.Kx;
    Ky = table.Ky;
    t = table.t;
    [MRSIStruct, prevPermute, prevShape] = reshapeDimensions(MRSIStruct, {'t', 'ky', 'kx'});



        
end
