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
center = st.centre;
%dimension size of the MRI
dims = st.vols{1}.dim;
%resolution of MRI (in mm)(not used might be useful)
resolution = (st.bb(2,:) - st.bb(1,:) + 1)/dims(:)';
%3D bounding box of the MRI scan in mm. ie the coordinates where the MRI is
%plotted onto.
bb = st.bb;

%counters
s = 1;
c = 1;
t = 1;
%initalize variables (TODO: Figure out how to preallocate size... Does it
%save time?)
sagital_voxels(numel(voxels)) = Voxel();
coronal_voxels(numel(voxels)) = Voxel();
transverse_voxels(numel(voxels)) = Voxel();
for i = 1:numel(voxels)
    %check which voxels intersect with the current sagital plane. 
    if(voxels(i).min_coords(1) <= center(1) && voxels(i).max_coords(1) >= center(1))
        sagital_voxels(s) = voxels(i);
        s = s + 1;
    end
    %check which voxels intersectwith the axial plane
    if(voxels(i).min_coords(2) <= center(2) && voxels(i).max_coords(2) >= center(2))
        coronal_voxels(c) = voxels(i);
        c = c + 1;
    end
    %check which voxels intersect with the coronal plane
    if(voxels(i).min_coords(3) <= center(3) && voxels(i).max_coords(3) >= center(3))
        transverse_voxels(t) = voxels(i);
        t = t + 1;
    end
end
if(s < numel(voxels))
    sagital_voxels(s:numel(voxels)) = [];
end
if(c < numel(voxels))
    coronal_voxels(c:numel(voxels)) = [];
end
if(t < numel(voxels))
    transverse_voxels(t:numel(voxels)) = [];
end

%get a logical array of points to plot for the ppm range
range_bool = CSI_OBJ.ppm >= ppm_min & CSI_OBJ.ppm <= ppm_max;
ppm = CSI_OBJ.ppm(range_bool);
ppm = ppm - min(ppm);

%permute specs
specs = permute(getData(CSI_OBJ), nonzeros([CSI_OBJ.dims.t, CSI_OBJ.dims.x, CSI_OBJ.dims.y,...
                            CSI_OBJ.dims.z, CSI_OBJ.dims.coils, CSI_OBJ.dims.averages]));
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
%temp variables
min_amp = realmax;
max_amp = realmin;
plot_specs = zeros(sum(range_bool), numel(transverse_voxels));
for i = 1:numel(transverse_voxels)
    %TODO: This will throw an error if there is another dimension bigger than one other
    %than x, y, and z.
    plot_specs(:,i) = flip(specs(range_bool, transverse_voxels(i).index(1),...
        transverse_voxels(i).index(2), transverse_voxels(i).index(3), 1), 1);
    
    plot_specs(:,i) = plot_specs(:,i) - plot_specs(end,i);

    
    %find min and maximum amplitudes
    if(min_amp > min(plot_specs(:,i)))
        min_amp = min(plot_specs(:,i));
    end
    if(max_amp < max(plot_specs(:,i)))
        max_amp = max(plot_specs(:,i));
    end
    
end
%change to plotting type (ie. imaginary, real, or absolute)
yrange = max_amp - min_amp;

%get the range of x and y values
xrange = ppm_max - ppm_min;

%scale factors to fit the spectral dimension at each x and y coordinates
min_points = min(voxels(1).points,[], 2);
max_points = max(voxels(1).points,[], 2);
x_pix_width = max_points(1)-min_points(1);
y_pix_width = max_points(2)-min_points(2);
scalefactorX=(0.8*(x_pix_width))/xrange;
scalefactorY=(0.8*(y_pix_width))/yrange;

ppm_plot = repmat(ppm', [1, numel(transverse_voxels)]);

%scale axes by x and y scale factors and shift by offsets
plot_specs = plot_specs .* scalefactorY - y_pix_width/2 - bb(1,2);
ppm_plot = ppm_plot .* scalefactorX - x_pix_width/2 - bb(1,1);

if (numel(transverse_voxels) > 0)
    vox_centers = [transverse_voxels.center];
    plot_specs = plot_specs + vox_centers(1, :);
    ppm_plot = ppm_plot + vox_centers(2,:);
end


h = findobj(st.vols{1}.ax{1}.ax,'Tag','csi_plot');
delete(h);

if(show_vox)
    for i = 1:numel(transverse_voxels)
        coordinates = transverse_voxels(i).find_intersection(3, center(3)) - bb(1,1:2)';
        patch(st.vols{1}.ax{1}.ax, coordinates(1, :), coordinates(2,:), ...
        'w', 'FaceAlpha', 0, 'LineWidth', 1, 'Tag', 'csi_plot', 'EdgeColor', 'g', 'HitTes', 'off');
    end
end
hold(st.vols{1}.ax{1}.ax, 'on');
line(st.vols{1}.ax{1}.ax, ppm_plot, plot_specs, 'Tag', 'csi_plot', 'HitTes', 'off');
hold(st.vols{1}.ax{1}.ax, 'off');
%stop holding

%PLOTTING IN THE SAGITAL PLANE
plot_plane(st.vols{1}.ax{3}.ax, sagital_voxels, 1, center);

%PLOTTING ON THE CORONAL PLANE
plot_plane(st.vols{1}.ax{2}.ax, coronal_voxels, 2, center);
end
