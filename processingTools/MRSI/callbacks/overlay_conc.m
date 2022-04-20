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

function overlay_conc(in, met, x_range, y_range)
    %call global variable from spm
    global st;
    %default value of 1 for coilNum

    % Need to be global since callback happens ion mouse down and st updates every time
    global CSI_OBJ
    global METABOLITE
    global X
    global Y
    global voxels

    changed = false;

    if(exist('met', 'var'))
        METABOLITE = met;
        changed = true;
    end
    if(exist('in','var'))
        CSI_OBJ = in;
    end
    if(exist('x_range', 'var'))
        X = x_range;
    end
    if(exist('y_range','var'))
        Y = y_range;
    end
    if(isempty(METABOLITE))
        return
    end

    %center in the dimensions (l-r,a-p,s-i)
    center = st.centre;
    %dimension size of the MRI
    dims = st.vols{1}.dim;
    %resolution of MRI (in mm)(not used, might be useful)
    resolution = (st.bb(2,:) - st.bb(1,:) + 1)/dims(:)';

    s = 1;
    c = 1;
    t = 1;
    sagital_voxels(numel(voxels)) = Voxel();
    coronal_voxels(numel(voxels)) = Voxel();
    transverse_voxels(numel(voxels)) = Voxel();
    for i = 1:numel(voxels)
        if(voxels(i).minimumCoordinates(1) <= center(1) && voxels(i).maximumCoordinates(1) >= center(1))
            sagital_voxels(s) = voxels(i);
            s = s + 1;
        end
        if(voxels(i).minimumCoordinates(2) <= center(2) && voxels(i).maximumCoordinates(2) >= center(2))
            coronal_voxels(c) = voxels(i);
            c = c + 1;
        end
        if(voxels(i).minimumCoordinates(3) <= center(3) && voxels(i).maximumCoordinates(3) >= center(3))
            transverse_voxels(t) = voxels(i);
            t = t + 1;
        end
    end
    if(s < numel(voxels))
        sagital_voxels(s:end) = [];
    end
    if(c < numel(voxels))
        coronal_voxels(c:end) = [];
    end
    if(t < numel(voxels))
        transverse_voxels(t:end) = [];
    end


    %3D bounding box of the MRI scan in mm. ie the coordinates where the MRI is
    %plotted onto.
    bb = st.bb;


    %PLOTING IN THE TRANSVERSE PLANE
    %set current axis object to be of the transverse image. (st.vols{1}.ax{2} and
    %st.vols{1}.ax{3} are objects for the sagital and coronal)

    %don't plot unless in the correct z position
    trans_plot = findobj(st.vols{1}.ax{4}.ax,'Tag','trans_plot');
    delete(trans_plot);
    %(needs to be modified if 3D MRSI is to be done)
    coordinates = zeros(2, 5, numel(transverse_voxels));
    counter = 1;
    mets_to_plot = zeros(1,numel(transverse_voxels));
    for i = 1:numel(transverse_voxels)
        if(transverse_voxels(i).index(1) <= X(2) && transverse_voxels(i).index(1)  >= X(1))
            if(transverse_voxels(i).index(2) <= Y(2) && transverse_voxels(i).index(2) >= Y(1))
                coordinates(:,:,counter) = transverse_voxels(i).find_intersection('axial', center(3)) - bb(1,1:2)';
                x_index = transverse_voxels(i).index(1) - X(1) + 1;
                y_index = transverse_voxels(i).index(2) - Y(1) + 1;
                mets_to_plot(counter) = METABOLITE(x_index, y_index);
                counter = counter + 1;
            end
        end
    end
    if(counter < numel(transverse_voxels))
        coordinates(:, :, counter:end) = [];
        mets_to_plot(counter:end) = [];
    end
    patch(st.vols{1}.ax{4}.ax, squeeze(coordinates(1,:,:)), squeeze(coordinates(2,:,:)), mets_to_plot, 'Tag', 'trans_plot');

    %PLOTTING IN THE SAGITAL PLANE
    plot_plane(st.vols{1}.ax{3}.ax, sagital_voxels, 'sagital', center);

    %PLOTTING ON THE CORONAL PLANE
    plot_plane(st.vols{1}.ax{2}.ax, coronal_voxels, 'coronal', center);
end

