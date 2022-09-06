%select_spec.m 
%
% callback function used for uicontrol popup menu for choosing to select
% MRSI specs to plot. Plots the closest clicked MRSI spec in spm_image GUI.
% must be used with sim_CSIoverlayMRI and spm_image. See
% https://www.mathworks.com/help/matlab/creating_plots/callback-definition.html
% for more info on callback funcitons.
%
% USAGE: can't be used by itself, must be called back
% 
% INPUT: 
% src       = modified uicontrol object
%

function select_spec(src,~) 
    %get global MRSI object
    global CSI_OBJ;
    global ppm_min;
    global ppm_max;
    global type;
    global voxels;
    global st;

    %only start if the users clicks with the right mouse button
    if ~strcmpi(get(gcbf,'SelectionType'),'alt')
        %get where the user clicked on the plot
        mouseSelectCoordinates = get(src,'CurrentPoint');
        %get x and y coordinates
        mouseSelectCoordinates = mouseSelectCoordinates(1:2:3) + st.bb(1,1:2);
        %get all the MRSI line plots
        
        voxelCenterCoordinates = [voxels.center];
        distance = mouseSelectCoordinates' - voxelCenterCoordinates(1:2,:);
        [~, indexOfClosestVoxel] = min(sum(distance.^2, 1));
        %display the linear index of the closest MRSI spec
        %get the corresponding x,y indecies from linear indexing
        
        %delete the current MRSI shown plot if there is one
        ax = findobj(gcf, 'Tag', 'csi_axis');
        if(numel(ax.Children) ~= 0)
            delete(ax.Children)
        end
        %hold the axis
        hold(ax, 'on')
        %plot the selected MRSI spec.4
        range_bool = CSI_OBJ.ppm > ppm_min & CSI_OBJ.ppm < ppm_max;
        ppm = CSI_OBJ.ppm(range_bool);
        set(ax,'XDir','reverse')
        xlim(ax, [ppm_min ppm_max]);
        specs = permute(getData(CSI_OBJ), nonzeros([CSI_OBJ.dims.t, CSI_OBJ.dims.x, CSI_OBJ.dims.y,...
                            CSI_OBJ.dims.z, CSI_OBJ.dims.coils, CSI_OBJ.dims.averages]));
        specs = specs(range_bool, voxels(indexOfClosestVoxel).index(1), voxels(indexOfClosestVoxel).index(2), voxels(indexOfClosestVoxel).index(3));
        switch(type)
            case 'real'
                specs = real(specs);
            case 'imaginary'
                specs = imag(specs);
            case 'magnitude'
                specs = abs(specs);
            otherwise
                error("please enter a valid plot_type");
        end
        
        plot(ax, ppm, specs)
        hold(ax, 'off')
        
    end 

end
