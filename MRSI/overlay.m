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

function overlay(plot_type, in, coilNum)

%call global variable from spm
global st;
%call global variable for MRSI object
%(needed for callbacks) 
global csi_obj;

%replace global csi_obj if a new MRSI object is defined
if exist('in', 'var')
    csi_obj = in;
end

%default value of 1 for coilNum
if ~exist('coilNum', 'var')
    coilNum = 1;
end

if ~exist('plot_type', 'var')
    %default is to plot the real
    plot_type = 'real';
else
    %ensure plot_type is lowercase
    plot_type = char(plot_type);
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


%get tne yrange for scaling of the y axis to plot the specs 
if(csi_obj.dims.coils == 0)
    yrange=max(real(csi_obj.specs),[],'all') - min(real(csi_obj.specs),[],'all');
else
    yrange=max(real(csi_obj.specs(:,coilNum,:,:)),[],'all') - min(real(csi_obj.specs(:,coilNum,:,:)),[],'all');
end

%scale factors to fit the spectral dimension at each (x,y) coordinates
scalefactorX=(0.8*csi_obj.deltaX)/max(csi_obj.t);
scalefactorY=(0.8*csi_obj.deltaY)/yrange;

%set the x dimension and y dimension
xdim = csi_obj.dims.x;
ydim = csi_obj.dims.y;

%reduce the spec scales by scalefactorY
tempSpec=op_ampScale(csi_obj,scalefactorY);

%reverse the y dimension so plotting anterior = +y and posterior = -y
tempSpec.specs = flip(tempSpec.specs, ydim);

%rearange dimensions so that we plot time, x, y
dimsToPlot = [csi_obj.dims.t, xdim, ydim];
extraDims = setdiff(numel(size(csi_obj.sz)), dimsToPlot);
tempSpec = permute(tempSpec.specs, [dimsToPlot, extraDims]);
tempSpec = reshape(tempSpec, [csi_obj.sz(dimsToPlot), prod(csi_obj.sz(extraDims))]);



%set current axis object to be of the axial image. (st.vols{1}.ax{2} and
%st.vols{1}.ax{3} are objects for the sagital and coronal)
axes(st.vols{1}.ax{1}.ax)
%delete previous MRSI plot on MRI
h = findobj(gca,'Tag','csi_plot');
delete(h)
%hold axis
hold on;
%don't plot unless in the correct z position
%(needs to be modified if 3D MRSI is to be done) 
if(center(3) - csi_obj.deltaZ/2 < csi_obj.imageOrigin(3) && csi_obj.imageOrigin(3) < center(3) + csi_obj.deltaZ/2)
    for x = 1:size(csi_obj.specs,xdim)
        for y = 1:size(csi_obj.specs,ydim)
            %first scale the ppm scale so that the range is correct;
            time = csi_obj.t*scalefactorX;
            %now shift it to the correct x-position;
            time = time + (x-1)*csi_obj.deltaX - (0.8*csi_obj.deltaX)/2 + csi_obj.xCoordinates(1) - bb(1,1);
            %check plottype
            fids = tempSpec(:,x,y,coilNum);
            
            %plot whichever plot type is chosen
            switch(plot_type)
                case 'real'
                    fids = real(fids);
                case 'imaginary'
                    fids = imag(fids);
                case 'magnitude'
                    fids = abs(fids);
                otherwise
                    error("please enter a valid plot_type");
            end
            
            %Now start plotting
            p = plot(gca, time,fids + (y-1)*csi_obj.deltaY + csi_obj.yCoordinates(1)+ -bb(1,2));
            %set the label of the line object
            set(p, 'Tag', 'csi_plot', 'HitTes', 'off');
        end
    end
end
%stop holding
hold off;


end