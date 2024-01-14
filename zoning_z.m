% This function divide the volumetric image into smaller subvolumes along
% the z axis
function [list_z, list_z_im] = zoning_z(z,dz,dz_im,n_sub,overlap2,z_im)
    % Thickness of the image
    l = length(z)*dz;
    % Half size of the overlap
    d = 0;
    for sub_id = 1:n_sub
        % Equally divide the image into subvolumes that are not overlapped.
        % Then, we will expand the boundary of the subvolumes to make them 
        % overlap. The inner subvolumes are expanded in both directions while 
        % the two outer subvolumes are expanded in only one direction. We 
        % will set the boundaries first. Here are the upper and lower
        % bounds of the zones when we run optimization and when
        % reconstructing the image
        low_bound = (sub_id-1)*l/n_sub+min(z,[],'all');
        up_bound = sub_id*l/n_sub+min(z,[],'all');
        low_bound_im = (sub_id-1)*l/n_sub+min(z_im,[],'all');
        up_bound_im = sub_id*l/n_sub+min(z_im,[],'all');
        % Then expand the bounds
        if sub_id > 1
            low_bound = low_bound-d;
            low_bound_im = low_bound_im-d;
        end
        if sub_id < n_sub
            up_bound = up_bound+d;
            up_bound_im = up_bound_im+d;
        end
        % Now, assign the coordinates
        list_z{sub_id,1} = low_bound:dz:up_bound-dz/2;
        list_z_im{sub_id,1} = low_bound_im:dz_im:up_bound_im-dz_im/2;
    end
end