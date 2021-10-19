%Helper function to be used to plot voxels onto coronal plane

function plot_plane(ax, vox, dims, center)
%Get the CSI object
global st
bb = st.bb;

if(dims == 1)
    tag = 'sagital';
elseif(dims == 2)
    tag = 'coronal';
end
offset = [bb(1,1), -bb(2,2), bb(1,3)]';
%Get plots if they exis

%see if we are in the correct position to plot
all_dims = [1,2,3];
diff = setdiff(all_dims, dims);

plane = findobj(ax, 'Tag', tag);
delete(plane);
voxel_positions = zeros(2, 5, numel(vox));
for i = 1:numel(vox)
    if(dims == 1)
        coordinates = vox(i).find_intersection(dims, center(dims));
        coordinates(1,:) = -coordinates(1,:);
    else
        coordinates = vox(i).find_intersection(dims, center(dims));
    end
    coordinates = coordinates - offset(diff);
    voxel_positions(:,:, i) = coordinates;
 

end
hold(ax, 'on');
    patch(ax, squeeze(voxel_positions(1,:, :)), squeeze(voxel_positions(2,:, :)), ...
        'w', 'FaceAlpha', 0, 'LineWidth', 1, 'Tag', tag, 'EdgeColor', 'g', 'HitTes', 'off');
hold(ax, 'off');
