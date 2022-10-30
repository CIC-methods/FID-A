% op_CSIoverlayMRI.m
% Brenden Kadota, Sunnybrook 2022
%
% DESCRIPTION:
% Overlay MRSI structure on spm_image graphic. Displays the MRSI specs onto
% the axial plane of the brain.
%
% USAGE:
% op_CSIoverlayMRI('/path/to/nifit/mri.nii', MRSIStruct, 1);
%
%
% INPUT:
% mriFileName           = string path to nifti file
% in                    = MRSI structure to overlay. In FID-A format
% coilIndex             = coilIndex to plot if coils have not been combined (default value = 1)
% averageIndex          = average index to plot
%
% OUTPUT:
%


function op_CSIoverlayMRI(MRSIStruct, mriFileName, coilIndex, averageIndex)
    arguments
        MRSIStruct (1, 1) struct
        mriFileName (:, 1) string
        coilIndex (1, 1) double = 1
        averageIndex (1, 1) double = 1
    end

    checkArguments(MRSIStruct);
    %call spm global variable
    global st
    global voxels
    global showVoxelOutlineFlag

    showVoxelOutlineFlag = false;

    % ensure filename is a character array not string
    mriFileName = char(mriFileName);
    if(contains(mriFileName, '~'))        
        mriFileName = fullfile(getuserdir, mriFileName(2:end));
    end
    % display image using spm
    spm_image('display', mriFileName)

    % adds MRSI overlay function to SPM callback function
    st = addMRSIOverlayToSPMCallback(st, 'overlaySpectraOnMRI()');

    spmFigure = st.fig;
    buildMouseClickDropdown(spmFigure);
    buildPPMBoundTextBox(spmFigure);
    buildPlotTypeDropdown(spmFigure);
    buildShowVoxelButton(spmFigure);

    voxels = buildArrayOfVoxels(MRSIStruct);
    %plot the MRSI onto the spm mri graphic along the axial dimension.
    overlaySpectraOnMRI('real', min(MRSIStruct.ppm), max(MRSIStruct.ppm), MRSIStruct, coilIndex, averageIndex);
end

function change_plot(src,~)
    overlaySpectraOnMRI(src.String{src.Value});
end

function set_min(src, ~)
    max_obj = findobj('Tag', 'max');
    max_obj.String;
    ppm_min = str2double(src.String);
    ppm_max = str2double(max_obj.String);
    if(~isempty(ppm_min) && ~isempty(ppm_max))
        if(ppm_min <= ppm_max)
            plot_obj = findobj('Tag', 'plot_type');
            cur_type = plot_obj.String{plot_obj.Value};
            overlaySpectraOnMRI(cur_type, ppm_min, ppm_max);
        end
    end
end

function set_max(src, ~)
    min_obj = findobj('Tag', 'min');
    min_obj.String;
    ppm_max = str2double(src.String);
    ppm_min = str2double(min_obj.String);
    if(~isempty(ppm_min) && ~isempty(ppm_max))
        if(ppm_min <= ppm_max)
            plot_obj = findobj('Tag', 'plot_type');
            cur_type = plot_obj.String{plot_obj.Value};
            overlaySpectraOnMRI(cur_type, ppm_min, ppm_max);
        end
    end
end

function setMouseClickOption(src, ~)
    global st

    %get the axis object for the axial plane
    axis = st.vols{1}.ax{1}.ax;
    
    if(src.Value == 2)
        %If the popup menu is set to select_spec
        %set the callback function of click on axis to the method
        %select_spec.m
        set(axis, 'ButtonDownFcn', @select_spec);
        %position new axis here with properties
        axes('Visible','off', ...
             'Parent', st.fig, ...
             'YDir','normal', ...
             'Units','Pixels', ...
             'Box','on', ...
             'Position', [380 445 250 220], ...
             'Units','normalized', ...
             'Visible','on', ...
             'Tag', 'csi_axis');
    else
        %If the popup menu is set to reposition
        %set click on axial plane to reposition
        set(axis, 'ButtonDownFcn', @repos_start);
        %find the axis used to plot MRSI specs
        ax = findobj(gcf, 'Tag', 'csi_axis');
        %delete the specs
        delete(ax);
    end
end

%callback method used for repositioning
function repos_start(~, ~)
    currentFigure = gcbf;
    %checking if the click is a right click
    if ~strcmpi(get(currentFigure, 'SelectionType'),'alt')
        %set callback functions for movement of mouse and finished
        %clicking. spm_orthviews('reposition') is for movement of mouse
        % and repos_end for end of click.
        set(currentFigure, 'windowbuttonmotionfcn',@()(spm_orthviews('reposition')), 'windowbuttonupfcn',@repos_end);
        %call reposition in spm_orthviews
        spm_orthviews('reposition');
    end
end

function repos_end(~, ~)
    %remove callback functions for movement of mouse and wiindows
    %button up.
    set(gcbf,'windowbuttonmotionfcn','', 'windowbuttonupfcn','');
end

function showVoxelCallback(src, ~)
    global show_vox
    show_vox = src.Value;
    overlaySpectraOnMRI();
end

function checkArguments(in)
    if ~exist('spm_image', 'file')
        error('please add spm12 to your path before continuing');
    end
    checkSpatialFT(in)
end

function buildMouseClickDropdown(spmFigure)
    %Labels for dropdown
    labels={'Reposition', 'Select Spec'};
    %putting the dropdown menu in the figure to reposition or select spec
    mouseDropdown = uicontrol('Parent', spmFigure, ...
        'Units', 'Pixels', ...
        'Position', [320 400 100 20], ...
        'Style', 'Popupmenu', ...
        'String', labels);
    %selector callback function
    mouseDropdown.Callback = @setMouseClickOption;
end

function [maxPPMTextBox, minPPMTextBox] = buildPPMBoundTextBox(spmFigure)
    maxPPMTextBox = uicontrol('Parent', spmFigure, ...
        'Units', 'Pixels', ...
        'Position', [100, 405, 50, 20],...
        'Style', 'edit', ...
        'String', 'ppmmax', ...
        'Tag', 'max', ...
        'Callback', @set_max);
    
    minPPMTextBox = uicontrol('Parent', spmFigure, ...
        'Units', 'Pixels', ...
        'Position', [50, 405, 50, 20],...
        'Style', 'edit', ...
        'String', 'ppmmin', ...
        'Tag', 'min', ...
        'Callback', @set_min);
end


%
function buildPlotTypeDropdown(spmFigure)
    %labels for the dropdown menu
    plotLabels = {'real', 'imaginary', 'magnitude'};
    %creates a dropdown menu in the figure to decide to plot real, imaginary,
    %or magnitud
    plot_type_selector = uicontrol('Parent', spmFigure, ...
        'Units', 'Pixels', ...
        'Position', [200 400 100 20], ...
        'Style', 'Popupmenu', ...
        'String', plotLabels, ...
        'Tag', 'plot_type', ...
        'Callback', @change_plot);
end

function buildShowVoxelButton(spmFigure)
    uicontrol('Parent', spmFigure, ...
        'Units', 'Pixels', ...
        'Position', [500 405 60 20],...
        'Style', 'togglebutton', ...
        'String', 'Show Voxels', ...
        'Tag', 'show_voxels', ...
        'Callback', @showVoxelCallback);
end