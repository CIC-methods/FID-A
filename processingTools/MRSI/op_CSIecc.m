function [MRSISupressed, MRSIUnsupressed] = op_CSIecc(MRSISupressed, MRSIUnsupressed)
    [MRSISupressed, dimOrderSup, shapeSup] = reshapeDimensions(MRSISupressed, {'t', 'y', 'x'});
    [MRSIUnsupressed, dimOrderUnsup, shapeUnsup] = reshapeDimensions(MRSIUnsupressed, {'t', 'y', 'x'});
    
    dataUnsupressed = getData(MRSIUnsupressed);
    dataSupressed = getData(MRSIUnsupressed);
    for e = 1:getSizeFromDimensions(MRSISupressed, {'extras'})
        for y = 1:getSizeFromDimensions(MRSISupressed, {'y'})
            for x = 1:getSizeFromDimensions(MRSISupressed, {'x'})
                mrsUnsupressed = op_CSItoMRS(MRSIUnsupressed, x, y, 'extraIndex', e);
                mrsSupressed = op_CSItoMRS(MRSISupressed, x, y, 'extraIndex', e);
                [mrsUnsupressed, mrsSupressed] = op_ecc(mrsUnsupressed, mrsSupressed);
                if(getFlags(MRSISupressed, 'spectralFT') == 1)
                    dataUnsupressed(:, y, x, e) = mrsUnsupressed.specs;
                    dataSupressed(:, y, x, e) = mrsSupressed.specs;
                else
                    dataUnsupressed(:, y, x, e) = mrsUnsupressed.fids;
                    dataSupressed(:, y, x, e) = mrsSupressed.fids;
                end
            end
        end
    end
    MRSISupressed = setData(MRSISupressed, dataSupressed);
    MRSIUnsupressed = setData(MRSIUnsupressed, dataUnsupressed);

    MRSISupressed = reshapeBack(MRSISupressed, dimOrderSup, shapeSup);
    MRSIUnsupressed = reshapeBack(MRSIUnsupressed, dimOrderUnsup, shapeUnsup);

end