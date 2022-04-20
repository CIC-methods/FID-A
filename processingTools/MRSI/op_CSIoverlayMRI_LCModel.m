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

function op_CSIoverlayMRI_LCModel(in, nifti, lcTable)
    arguments
        in (1, 1) struct
        nifti (1, :) {mustBeFile}
        lcTable (1, :) {mustBeFile}
    end
    %call spm global variable
    global st
    global voxels

    %Some error checks
    if ~exist('spm_image', 'file')
        error('please add spm12 to your path before continuing');
    end
    if ~isa(nifti, 'char') && ~isa(nifti, 'string')
        error('please enter a string or char array for the nifti variable')
    end
    if ~isa(lcTable, 'char') && ~isa(lcTable, 'string')
        error('please enter a string')
    end

    met_data = readtable(lcTable);

    %ensure filename is a character array not string
    nifti = char(nifti);
    %display image using spm
    spm_image('display', nifti)

    %check if st.callback is a cell array of strings or a string
    if(ischar(st.callback))
        %create cell array and add overaly() callback onto it
        callback = st.callback;
        st.callback = cell(1,2);
        st.callback{1} = callback;
        st.callback{2} = 'overlay_conc()';
    else
        %append the overlay callback to the last spot in the cell array
        st.callback{size(st.callback, 2) + 1} = 'overlay_conc()';
    end

    conc = axes('Position', st.vols{1}.ax{1}.ax.Position, 'xtick', [], ...
        'xticklabel', [], 'ytick', [], 'yticklabel',[], 'Color', 'none', ...
        'XLim', st.vols{1}.ax{1}.ax.XLim, 'YLim', st.vols{1}.ax{1}.ax.YLim, 'PickableParts','none');
    st.vols{1}.ax{4}.ax = conc;
    colormap(conc, 'default');


    %labels for the dropdown menu. Remove first 2 rows because they are row and column indecies
    plot_labels = met_data.Properties.VariableNames(3:end);
    plot_labels = ['None', plot_labels];
    %putting the dropdown menu in the figure to reposition or select spec
    selector = uicontrol('Parent',st.fig,'Units','Pixels','Position',[320 405 100 20],...
        'Style','Popupmenu','String',plot_labels);

    %selector callback function
    selector.Callback = @display_met;
    function display_met(src, ~)
        met = src.Value;
        y_range(1) = min(met_data.Row);
        y_range(2) = max(met_data.Row);
        x_range(1) = min(met_data.Col);
        x_range(2) = max(met_data.Col);
        met_arr = zeros(x_range(2)-x_range(1)+1,y_range(2)-y_range(1)+1);
        if(met ~= 0)
            for i = 1:length(met_data.(met+1))
                x = met_data.Col(i)-x_range(1)+1;
                y = met_data.Row(i)-y_range(1)+1;
                met_arr(x,y) = met_data.(met+1)(i);
            end
        else
            met_arr = [];
        end
        overlay_conc(in, met_arr, x_range, y_range);
    end

    voxels = createVoxels(in);
    overlay_conc(in, []);
end


