%op_CSIoverlayMRI_LCModel.m
%Brenden Kadota, Jamie Near, Sunnybrook 2021.
%
% USAGE:
% RF=op_CSIoverlayMRI_LCModel(in,nifti,lcTable);
%
% DESCRIPTION:
% Takes csv output of lcmodel and maps the concentrations of each voxel
% onto an mri. Uses spm12 to display mri.
%
% INPUTS:
% in        = input data in FID-A CSI structure format.
% nifti     = mri in nifti format.
% lcTable   = lcmodel output file in csv format.
%
% NOTE: Lcmodel's output in csv has no information on position other than
% row and column numbers. Therefore, this code assumes neurological
% convention (right -> x+).

function op_CSIoverlayMRI_LCModel(MRSIStruct, mriFileNameInNifti, lcModelCSV)
    arguments
        MRSIStruct (1, 1) struct
        mriFileNameInNifti (1, :) char {mustBeFile}
        lcModelCSV (1, :) char {mustBeFile}
    end
    % call spm global variable
    global st

    checkArguments(mriFileNameInNifti, lcModelCSV);

    buildMRIPlotUsingSPM(mriFileNameInNifti);
    st = addMRSIOverlayToSPMCallback(st, 'overlayConcMapOverMRI');
    st = buildConcentrationMapAxes(st);

    voxels = buildArrayOfVoxels(MRSIStruct);
    
    spmFigure = st.fig;
    buildMetaboliteDropDownSelector(spmFigure, readtable(lcModelCSV));
    initalizeFigureUserData(spmFigure, voxels)
    
    overlayConcMapOverMRI;
end

function sendMetaboliteMapToPlot(src, ~, lcModelTable)
    selectedMetaboliteIndex = src.Value;
    spmFigure = src.Parent;
    
    [metabolitePlotArray, xIndexBounds, yIndexBounds] = makeLcModelHeatMap(lcModelTable, selectedMetaboliteIndex);
    spmFigure.UserData.metabolitePlot = metabolitePlotArray;
    spmFigure.UserData.xIndexBounds = xIndexBounds;
    spmFigure.UserData.yIndexBounds = yIndexBounds; 
    overlayConcMapOverMRI;
end

function checkArguments(mriFileNameInNifti, lcModelCSV)
    %Some error checks
    if ~exist('spm_image', 'file')
        error('please add spm12 to your path before continuing');
    end
    if ~isa(mriFileNameInNifti, 'char') && ~isa(mriFileNameInNifti, 'string')
        error('please enter a string or char array for the nifti variable')
    end
    if ~isa(lcModelCSV, 'char') && ~isa(lcModelCSV, 'string')
        error('please enter a string')
    end
end

function st = buildConcentrationMapAxes(st)
    transversePlaneAxes = st.vols{1}.ax{1}.ax;
    concentrationMapAxes = axes('Position', transversePlaneAxes.Position, ...
        'XLim', transversePlaneAxes.XLim, ...
        'YLim', transversePlaneAxes.YLim, ...
        'XTick', [], ...
        'XTickLabel', [], ...
        'ytick', [], ...
        'yticklabel', [], ...
        'Color', 'none', ...
        'PickableParts','none', ...
        'Colormap', 'default');
    colormap(concentrationMapAxes, "default");

    % add colormap
    st.vols{1}.ax{4}.ax = concentrationMapAxes;
end

function buildMetaboliteDropDownSelector(spmFigure, lcModelTable)
    metaboliteLabels = buildMetaboliteLabelsFromTable(lcModelTable);

    % Dropdown ui for selecting metabolite heat map to plot
    metaboliteDropDownSelector = uicontrol('Parent', spmFigure, ...
        'Units', 'Pixels', ...
        'Position', [320 405 100 20], ...
        'Style', 'Popupmenu', ...
        'String', metaboliteLabels);
    setCallbackForMetaboliteSelector(metaboliteDropDownSelector, lcModelTable);
end

function setCallbackForMetaboliteSelector(metaboliteDropDownSelector, lcModelTable)
    metaboliteDropDownSelector.Callback = {@sendMetaboliteMapToPlot, lcModelTable};
end

function buildMRIPlotUsingSPM(mriFileNameInNifti)
    spm_image('display', mriFileNameInNifti)
end

function initalizeFigureUserData(spmFigure, voxels)
    spmFigure.UserData.metabolitePlot = [];
    spmFigure.UserData.xIndexBounds = [];
    spmFigure.UserData.yIndexBounds = [];
    spmFigure.UserData.voxels = voxels;
end

function metaboliteLabels = buildMetaboliteLabelsFromTable(lcModelTable)
    %labels for the dropdown menu. Remove first 2 rows because they are row and column indecies
    metaboliteLabels = lcModelTable.Properties.VariableNames(3:end);
    metaboliteLabels = ['None', metaboliteLabels];
end