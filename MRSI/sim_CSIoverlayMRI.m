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
% mriFileName           = char array or string representing the file path to mriFile
% in                    = MRSI structure to overlay on MRI
% coilNum               = coilNum to plot if coils have not been combined
%                       (default value = 1)


function sim_CSIoverlayMRI(mriFileName, in, coilNum)

%call spm global variable
global st
if ~exist('coliNum', 'var')
    coilNum = 1;
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
    st.callback{2} = 'overlay()';
else
    %append the overlay callback to the last spot in the cell array
    st.callback{size(st.callback, 2) + 1} = 'overlay()';
end

%labels for the dropdown menu
labels={'Reposition', 'Select Spec'};
plot_labels = {'real', 'imaginary', 'magnitude'};
%putting the dropdown menu in the figure to reposition or select spec
selector = uicontrol('Parent',st.fig,'Units','Pixels','Position',[320 405 100 20],...
    'Style','Popupmenu','String',labels);

%selector callback function
selector.Callback = @selection; 
    
%creates a dropdown menu in the figure to decide to plot real, imaginary,
%or magnitud 
plot_type = uicontrol('Parent',st.fig,'Units','Pixels','Position',[200 405 100 20],...
    'Style','Popupmenu','String',plot_labels);
%assigning callback method used to change the plotting type
plot_type.Callback = @change_plot;
    %the callback method used to change the plotting type
    function change_plot(src,~)
        overlay(src.String{src.Value})
    end
%plot the MRSI onto the spm mri graphic along the axial dimension.    
overlay('real', in, coilNum);
end


    