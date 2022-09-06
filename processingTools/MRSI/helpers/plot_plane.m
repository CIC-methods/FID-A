%Helper function to be used to plot voxels onto coronal plane

function plot_plane(axis, intersectingVoxels, planeLabel, cursorPosition)
    
    % Get the spm plotting object
    global st

    % bouding box for the MRI
    boundingBox = st.bb;
    sagitalBoundingBox = boundingBox(:, 1);
    coronalBoundingBox = boundingBox(:, 2);
    axialBoundingBox = boundingBox(:, 3);
    
    % Cursor position
    cursorSagitalPosition = cursorPosition(1);
    cursorCoronalPosition = cursorPosition(2);
    cursorAxialPosition = cursorPosition(3);

    % get the proper coordinate dimension. Dimensions are [sagital, coronal,
    % axial]
    if(strcmp(planeLabel, 'sagital'))
        plottingDimensions = [2, 3];
    elseif(strcmp(planeLabel, 'coronal'))
        plottingDimensions = [1, 3];
    end

    % offset to center world coordinates with plotting coordinates
    coordinateOffset = [sagitalBoundingBox(1), coronalBoundingBox(1), axialBoundingBox(1)]';

    % Delete previous plot
    planePlot = findobj(axis, 'Tag', planeLabel);
    delete(planePlot);

    % go through all voxels and find where they intersect with the current plane
    voxelIntersectionPositions = cell(1, numel(intersectingVoxels));
    for iVoxel = 1:numel(intersectingVoxels)
        if(strcmp(planeLabel, 'sagital'))
            coordinates = intersectingVoxels(iVoxel).find_intersection(planeLabel, cursorSagitalPosition);
            % sagital plane is reversed so we need to swap the y coordinates
            coordinates(2, :) = -coordinates(2,:);
        else
            coordinates = intersectingVoxels(iVoxel).find_intersection(planeLabel, cursorCoronalPosition);
        end
        coordinates = coordinates(plottingDimensions, :) - coordinateOffset(plottingDimensions);
        voxelIntersectionPositions{iVoxel} = coordinates;
    end

    %plot voxel bounding box along the plane
    hold(axis, 'on');
    cellfun(@(voxel) patch(axis, voxel(1, :), voxel(2, :), ...
        'w', 'FaceAlpha', 0, 'LineWidth', 1, 'Tag', planeLabel, ...
        'EdgeColor', 'g', 'HitTes', 'off'), voxelIntersectionPositions)
    hold(axis, 'off');
end