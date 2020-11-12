function out = CSI(scanID)
    twix = io_loadCSI_twix(scanID);
    combined = op_combineCSICoils(twix);
    transformed = op_CSIFourierTransform(combined);
    out = transformed;
    sim_plotCSI(transformed);

end