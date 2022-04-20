% sim_CSIoverlayMRI.m
%
% Overlay MRSI structure on spm_image graphic. Displays the MRSI specs onto
% the axial plane of the brain.
%
% USAGE:
% sim_CSIoverlayMRI(mriFileName, in, (optinal) coilNum);
%
%
% INPUT:
% mriFileName           = char array or string representing the file path
%                         to mriFile in nifti format
% in                    = MRSI structure to overlay on MRI
% coilNum               = coilNum to plot if coils have not been combined
%                       (default value = 1)


function op_CSIoverlayMRI(mriFileName, in, coilNum, average_num)

    checkArguments(in);
    %call spm global variable
    global st
    global voxels
    global show_vox

    show_vox = 0;
    if ~exist('coliNum', 'var')
        coilNum = 1;
    end
    if ~exist('average_num', 'var')
        average_num = 1;
    end

    %ensure filename is a character array not string
    mriFileName = char(mriFileName);
    %display image using spm
    spm_image('display', mriFileName)

    %check if st.callback is a cell array of strings or a string
    if(ischar(st.callback))
        %create cell array and add overaly() callback onto it
        callback = st.callback;
        st.callback = cell(1,2);
        st.callback{1} = callback;
        st.callback{2} = 'overlay_specs()';
    else
        %append the overlay callback to the last spot in the cell array
        st.callback{size(st.callback, 2) + 1} = 'overlay_specs()';
    end

    %Labels for dropdown
    labels={'Reposition', 'Select Spec'};
    %putting the dropdown menu in the figure to reposition or select spec
    selector = uicontrol('Parent',st.fig,'Units','Pixels','Position',[320 400 100 20],...
        'Style','Popupmenu','String',labels);
    %selector callback function
    selector.Callback = @selection;

    ppm_max_control = uicontrol('Parent', st.fig, 'Units', 'Pixels', 'Position', [100, 405, 50, 20],...
        'Style', 'edit', 'String', 'ppmmax', 'Tag', 'max');
    ppm_max_control.Callback = @set_max;

    ppm_min_control = uicontrol('Parent', st.fig, 'Units', 'Pixels', 'Position', [50, 405, 50, 20],...
        'Style', 'edit', 'String', 'ppmmin', 'Tag', 'min');
    ppm_min_control.Callback = @set_min;

    %labels for the dropdown menu
    plot_labels = {'real', 'imaginary', 'magnitude'};
    %creates a dropdown menu in the figure to decide to plot real, imaginary,
    %or magnitud
    plot_type_selector = uicontrol('Parent',st.fig,'Units','Pixels','Position',[200 400 100 20],...
        'Style','Popupmenu','String',plot_labels, 'Tag', 'plot_type');

    %assigning callback method used to change the plotting type
    plot_type_selector.Callback = @change_plot;

    uicontrol('Parent',st.fig,'Units','Pixels','Position',[500 405 60 20],...
        'Style','togglebutton','String', 'Show Voxels', 'Tag', 'show_voxels', ...
        'Callback', @show_voxels);


    voxels = createVoxels(in);
    %plot the MRSI onto the spm mri graphic along the axial dimension.
    
    overlay_specs('real', min(in.ppm), max(in.ppm), in, coilNum, average_num);
end

function change_plot(src,~)
    overlay_specs(src.String{src.Value});
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
            overlay_specs(cur_type, ppm_min, ppm_max);
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
            overlay_specs(cur_type, ppm_min, ppm_max);
        end
    end
end

function selection(src, ~)
    global st
    %get the axis object for the axial plane
    axis = st.vols{1}.ax{1}.ax;

    if(src.Value == 2)
        %If the popup menu is set to select_spec
        %set the callback function of click on axis to the method
        %select_spec.m
        set(axis, 'ButtonDownFcn', @select_spec);
        %position new axis here with properties
        axes('Visible','off', 'Parent',st.fig, ...
            'YDir','normal', 'Units','Pixels', 'Box','on',...
            'Position',[380 445 250 220],...
            'Units','normalized',...
            'Visible','on', 'Tag', 'csi_axis');
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

    %checking if the click is a right click
    if ~strcmpi(get(gcbf,'SelectionType'),'alt')
        %set callback functions for movement of mouse and finished
        %clicking. spm_orthviews('reposition') is for movement of mouse
        % and repos_end for end of click.
        set(gcbf,'windowbuttonmotionfcn',@()(spm_orthviews('reposition')), 'windowbuttonupfcn',@repos_end);
        %call reposition in spm_orthviews
        spm_orthviews('reposition');
    end
end

function repos_end(~, ~)
    %remove callback functions for movement of mouse and wiindows
    %button up.
    set(gcbf,'windowbuttonmotionfcn','', 'windowbuttonupfcn','');
end

function show_voxels(src, ~)
    global show_vox
    show_vox = src.Value;
    overlay_specs();
end

function checkArguments(in)
    if ~exist('spm_image', 'file')
        error('please add spm12 to your path before continuing');
    end
    if ~getFlags(in, 'spatialFT')
        error('please fourier tranform along the spatial dimension before using')
    end
end