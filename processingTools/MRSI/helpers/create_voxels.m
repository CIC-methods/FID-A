function vox = create_voxels(in)
    if(in.dims.z == 0)
        z_size = 1;
    else
        z_size = in.sz(in.dims.z);
    end
    
    [x_ind,y_ind,z_ind] = meshgrid(0:(in.sz(in.dims.x)), ...
                                   0:(in.sz(in.dims.y)), ...
                                   0.5:z_size+0.5);
    coord_vector = [x_ind(:), y_ind(:), z_ind(:) ones(numel(z_ind), 1)];
    affineMatrix = getAffineMatrix(in);
    world_coord = affineMatrix * coord_vector';
    
    world_coord = reshape(world_coord(1:3, :), 3, size(x_ind,1), size(x_ind,2), size(x_ind, 3));

    vox(in.sz(in.dims.x), in.sz(in.dims.y), z_size) = Voxel();
    for z_pix = z_size
        for y_pix = 1:in.sz(in.dims.y)
            for x_pix = 1:in.sz(in.dims.x)
                x_coordinates = [world_coord(1, y_pix, x_pix, z_pix) world_coord(1, y_pix, x_pix+1, z_pix)];
                y_coordinates = [world_coord(2, y_pix, x_pix, z_pix) world_coord(2, y_pix+1, x_pix, z_pix)];
                z_coordinates = [world_coord(3, y_pix, x_pix, z_pix) world_coord(3, y_pix, x_pix, z_pix+1)];
                vox(x_pix, y_pix, z_pix) = Voxel(x_coordinates, y_coordinates, z_coordinates, [x_pix, y_pix, z_pix]);
            end
        end
    end
end
