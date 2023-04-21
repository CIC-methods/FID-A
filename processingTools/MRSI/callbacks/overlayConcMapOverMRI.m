% overlay_conc
%
% callback on mouse down event on spm_image figure
% should display concentration maps on mri
%
% INPUTS:
% in: CSI object
% met: metabolite array
% X: (1,2) vector with the first and last index of x coordinates
% Y: (1,2) vector with the first and last of y coordinates

function overlayConcMapOverMRI
    %call global variable from spm
    global st;
    %default value of 1 for coilNum
    spmFigure = st.fig;

    [metabolitePlot, xIndexBounds, yIndexBounds, voxels] = getCurrentUserDataValues(spmFigure);

    if(isempty(metabolitePlot))
        return
    end
    % cursorPosition in the dimensions (l-r,a-p,s-i)
    cursorPosition = st.centre;
    %dimension size of the MRI
    dimensionSize = st.vols{1}.dim;
    %resolution of MRI (in mm)(not used, might be useful)
    resolution = (st.bb(2,:) - st.bb(1,:) + 1)/dimensionSize(:)';

    [sagital_voxels, coronal_voxels, transverse_voxels] = findIntersectingVoxels(voxels, cursorPosition);

    %3D bounding box of the MRI scan in mm. ie the coordinates where the MRI is
    %plotted onto.
    mriBoundingBox = st.bb;

    deleteCurrentConcentrationMap(st);
    %(needs to be modified if 3D MRSI is to be done)
    voxelIntersectionPositions = {};
    counter = 1;
    metabolitesToPlotVectorized = zeros(1, numel(transverse_voxels));
    for i = 1:numel(transverse_voxels)
        xIndex = transverse_voxels(i).fid_aIndex(1) - xIndexBounds(1) + 1;
        yIndex = transverse_voxels(i).fid_aIndex(2) - yIndexBounds(1) + 1;
        if (xIndex > 0 && yIndex > 0 && xIndex <= size(metabolitePlot, 1) && yIndex <= size(metabolitePlot, 2))
            transverseIntersection = transverse_voxels(i).findIntersection('axial', cursorPosition(3));
            voxelIntersectionPositions{counter} = transverseIntersection(1:2, :) - mriBoundingBox(1,1:2)';
            metabolitesToPlotVectorized(counter) = metabolitePlot(xIndex, yIndex);
            counter = counter + 1;
        end
    end

    
    for i = 1:length(voxelIntersectionPositions)
        currentCoordiante = voxelIntersectionPositions{i};
        patch(st.vols{1}.ax{4}.ax, currentCoordiante(1, :), ...
            currentCoordiante(2, :), metabolitesToPlotVectorized(i), ...
            'Tag', 'trans_plot');
    end
    %PLOTTING IN THE SAGITAL PLANE
    plotVoxelsOnMRIPlane(st.vols{1}.ax{3}.ax, sagital_voxels, 'sagital', cursorPosition);

    %PLOTTING ON THE CORONAL PLANE
    plotVoxelsOnMRIPlane(st.vols{1}.ax{2}.ax, coronal_voxels, 'coronal', cursorPosition);
end

function [metabolitePlot, xIndexBounds, yIndexBounds, voxels] = getCurrentUserDataValues(spmFigure)
    figureUserData = spmFigure.UserData;
    metabolitePlot = figureUserData.metabolitePlot;
    xIndexBounds = figureUserData.xIndexBounds;
    yIndexBounds = figureUserData.yIndexBounds;
    voxels = figureUserData.voxels;
end

function deleteCurrentConcentrationMap(st)
    concentrationMap = findobj(st.vols{1}.ax{4}.ax,'Tag','trans_plot');
    delete(concentrationMap);
end
