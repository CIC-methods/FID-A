%overlay.m
% Used as a callback function in sim_CSIoverlayMRI to overaly CSI data onto
% MRI. The overlay is done with spm_image first, so spm needs to be installed.
% Coordinates are in RAS or neurological coordinates (ie, right = x, anterior = y, superior = z.
%
% USAGE:
% overlay((optional) plot_type, (optional) in, (optional) coilNum)
% (must call spm_image before using overlay)
%
% INPUT:
% plot_type         = a char array of 'real', 'magnitude', and 'imaginary
%                   to be plotted
% in                = MRSI object
% coilNum           = coilNum to plot if coils have not been combined.
%                   (default value)  1

function overlay_specs(plot_type, ppmmin, ppmmax, in, coilNum, average_num)
    %call global variable from spm
    global st;

    %these are global variables to define plotting on callbacks
    global type;
    global show_vox;
    global CSI_OBJ;
    global ppm_min;
    global ppm_max;
    global coil;
    global voxels;
    global average

    %update global variables
    if(exist('plot_type', 'var'))
        type = plot_type;
    end
    if(exist('in','var'))
        CSI_OBJ = in;
    end
    if(exist('ppmmin','var'))
        ppm_min = ppmmin;
    end
    if(exist('ppmmax','var'))
        ppm_max = ppmmax;
    end
    if(exist('coilNum','var'))
        coil = coilNum;
    end
    if(exist('average_num','var'))
        average = average_num;
    end

    %center in the dimensions (l-r,a-p,s-i)
    cursorPosition = st.centre;

    %dimension size of the MRI
    dimensionSizes_mm = st.vols{1}.dim;
    %3D bounding box of the MRI scan in mm. ie the coordinates where the MRI is
    %plotted onto.
    boundingBox_mm = st.bb;

    [sagitalVoxels, coronalVoxels, axialVoxels] = getInersectingVoxels(cursorPosition);

    %get a logical array of points to plot for the ppm range
    timeRange = CSI_OBJ.ppm >= ppm_min & CSI_OBJ.ppm <= ppm_max;
    ppm = CSI_OBJ.ppm(timeRange);
    ppm = ppm - min(ppm);

    %permute specs
    spectralData = permute(getData(CSI_OBJ), nonzeros([CSI_OBJ.dims.t, CSI_OBJ.dims.x, CSI_OBJ.dims.y,...
        CSI_OBJ.dims.z, CSI_OBJ.dims.coils, CSI_OBJ.dims.averages]));
    switch(type)
        case 'real'
            spectralData = real(spectralData);
        case 'imaginary'
            spectralData = imag(spectralData);
        case 'magnitude'
            spectralData = abs(spectralData);
        otherwise
            error("please enter a valid plot_type");
    end

    %temp variables
    min_amp = realmax;
    max_amp = realmin;
    spectrumToPlot = zeros(length(ppm), numel(axialVoxels));
    for iVoxel = 1:numel(axialVoxels)

        currentVoxel = axialVoxels(iVoxel);
        voxelIndex = currentVoxel.index;
        spectrumToPlot(:, iVoxel) = spectralData(timeRange, voxelIndex(1), ...
        voxelIndex(2), voxelIndex(3), 1);
        spectrumToPlot(:, iVoxel) = flip(spectrumToPlot(:, iVoxel));

        spectrumToPlot(:, iVoxel) = spectrumToPlot(:,iVoxel) - spectrumToPlot(end, iVoxel);


        %find min and maximum amplitudes
        if(min_amp > min(spectrumToPlot(:,iVoxel)))
            min_amp = min(spectrumToPlot(:,iVoxel));
        end
        if(max_amp < max(spectrumToPlot(:,iVoxel)))
            max_amp = max(spectrumToPlot(:,iVoxel));
        end

    end
    %change to plotting type (ie. imaginary, real, or absolute)
    yrange = max_amp - min_amp;

    %get the range of x and y values
    xrange = ppm_max - ppm_min;

    %scale factors to fit the spectral dimension at each x and y coordinates
    minSagital = voxels(1).getVoxelMinimumCoordinate('sagital');
    maxSagital = voxels(1).getVoxelMaximumCoordinate('sagital');
    minCoronal = voxels(1).getVoxelMinimumCoordinate('coronal');
    maxCoronal = voxels(1).getVoxelMaximumCoordinate('coronal');
    sagitalVoxelWidth = maxSagital - minSagital;
    coronalVoxelWidth = maxCoronal - minCoronal;
    scalefactorX = (0.8 * sagitalVoxelWidth) / xrange;
    scalefactorY = (0.8 * coronalVoxelWidth) / yrange;

    ppmAxisToPlot = repmat(ppm', [1, numel(axialVoxels)]);

    %scale axes by x and y scale factors and shift by offsets
    spectrumToPlot = spectrumToPlot .* scalefactorY + coronalVoxelWidth - boundingBox_mm(1,2);
    ppmAxisToPlot = ppmAxisToPlot .* scalefactorX - sagitalVoxelWidth/2 - boundingBox_mm(1,1);

    if (numel(axialVoxels) > 0)
        voxelCenterCoordinates = [axialVoxels.center];
        sagitalCenters = voxelCenterCoordinates(1, :);
        coronalCenters = voxelCenterCoordinates(2, :);
        spectrumToPlot = spectrumToPlot + coronalCenters;
        ppmAxisToPlot = ppmAxisToPlot + sagitalCenters;
    end


    h = findobj(st.vols{1}.ax{1}.ax,'Tag','csi_plot');
    delete(h);

    if(show_vox)
        for iVoxel = 1:numel(axialVoxels)
            box = axialVoxels(iVoxel).find_intersection('axial', cursorPosition(3));
            box = box(1:2, :);
            coordinates = box - boundingBox_mm(1,1:2)';
            patch(st.vols{1}.ax{1}.ax, coordinates(1, :), coordinates(2,:), ...
                'w', 'FaceAlpha', 0, 'LineWidth', 1, 'Tag', 'csi_plot', 'EdgeColor', 'g', 'HitTes', 'off');
        end
    end
    hold(st.vols{1}.ax{1}.ax, 'on');
    line(st.vols{1}.ax{1}.ax, ppmAxisToPlot, spectrumToPlot, 'Tag', 'csi_plot', 'HitTes', 'off');
    hold(st.vols{1}.ax{1}.ax, 'off');
    %stop holding

    %PLOTTING IN THE SAGITAL PLANE
    plot_plane(st.vols{1}.ax{3}.ax, sagitalVoxels, 'sagital', cursorPosition);

    %PLOTTING ON THE CORONAL PLANE
    plot_plane(st.vols{1}.ax{2}.ax, coronalVoxels, 'coronal', cursorPosition);
end

% see if plane coordinate intersects with lower and upper coord bounds in 1
% dimension


function isIntersect = intersect(lowerCoords, upperCoords, plane)

    isIntersect = upperCoords >= plane && lowerCoords <= plane;
end

function [sagitalVoxels, coronalVoxels, axialVoxels] = getInersectingVoxels(center)
    currentSagitalposition = center(1);
    currentCoronalPosition = center(2);
    currentAxialPosition = center(3);
    global voxels
    %counters
    sagitalVoxelCounters = 1;
    coronalVoxelCounters = 1;
    axialVoxelCounters = 1;
    %initalize variables 
    sagitalVoxels(numel(voxels)) = Voxel();
    coronalVoxels(numel(voxels)) = Voxel();
    axialVoxels(numel(voxels)) = Voxel();
    for iVoxel = 1:numel(voxels)
        currentVoxel = voxels(iVoxel);
        %check which voxels intersect with the current sagital plane.
        if(intersect(currentVoxel.getVoxelMinimumCoordinate('sagital'), ...
                    currentVoxel.getVoxelMaximumCoordinate('sagital'), ...
                    currentSagitalposition))
            sagitalVoxels(sagitalVoxelCounters) = currentVoxel;
            sagitalVoxelCounters = sagitalVoxelCounters + 1;  
        end
        %check which voxels intersectwith the sagital plane
        if(intersect(currentVoxel.getVoxelMinimumCoordinate('coronal'), ...
                    currentVoxel.getVoxelMaximumCoordinate('coronal'), ...
                    currentCoronalPosition))
            coronalVoxels(coronalVoxelCounters) = currentVoxel;
            coronalVoxelCounters = coronalVoxelCounters + 1;
        end
        %check which voxels intersect with the axial plane
        if(intersect(currentVoxel.getVoxelMinimumCoordinate('axial'), ...
                    currentVoxel.getVoxelMaximumCoordinate('axial'), ...
                    currentAxialPosition))
            axialVoxels(axialVoxelCounters) = currentVoxel;
            axialVoxelCounters = axialVoxelCounters + 1;
        end
    end
    if(sagitalVoxelCounters < numel(voxels))
        sagitalVoxels(sagitalVoxelCounters:numel(voxels)) = [];
    end
    if(coronalVoxelCounters < numel(voxels))
        coronalVoxels(coronalVoxelCounters:numel(voxels)) = [];
    end
    if(axialVoxelCounters < numel(voxels))
        axialVoxels(axialVoxelCounters:numel(voxels)) = [];
    end
end