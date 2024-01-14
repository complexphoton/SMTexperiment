% Function to reconstruct the fully corrected 2D SMT image
function [I] = reconstruct2D(list_x_im,list_y_im,x_im,y_im,overlap2,k,n_division,r_z,directory_save,rad_order_start,rad_order_inc)

    I = zeros(length(y_im),length(x_im));
    subvolume_id = 1;
    n_zone = 2^(n_division*2); % Number of zones
    dx_im = round(abs(x_im(2)-x_im(1)),1); % Pixel size
    % Input and output k:
    kout_x = k(:,1); kout_y = k(:,2);
    kin_x = k(:,1).'; kin_y = k(:,2).';
    % The Fourier components
    fx = single(kout_x-kin_x);
    fy = single(kout_y-kin_y);
    % Build the final image of each zone
    for zone_id = 1:n_zone
        % Zone coordinate
        x = list_x_im{zone_id,n_division+1};
        y = list_y_im{zone_id,n_division+1};
        Nx = length(x); % Number of pixels in x and y of the zone
        Ny = length(y);
        [X,Y] = meshgrid(x,y); % Make grid points
        nx = round((x+dx_im/2)/dx_im); % Zone coordinate in pixels
        ny = round((y+dx_im/2)/dx_im);
        % Find the final aberration phase of the zone by summing the phase
        % at each division step
        N = length(k); % Number of k
        phi_in = single(zeros(N,1)); % Initialize the input and output aberration phase
        phi_out = single(zeros(N,1));
        for division_step = n_division:-1:0
            if division_step == n_division
                zone = zone_id; 
            else
                zone = ceil(zone/4);
            end
            c_in_step = load(""+directory_save+"c_in_"+division_step+"_zone_"+zone+"_subvolume_"+subvolume_id+".mat").c_in;
            c_out_step = load(""+directory_save+"c_out_"+division_step+"_zone_"+zone+"_subvolume_"+subvolume_id+".mat").c_out;
            Z_step =  single(load(""+directory_save+"Z_"+division_step+"_re.mat").Z_re);
            phi_in_step = Z_step*c_in_step;
            phi_out_step = Z_step*c_out_step;
            phi_in = phi_in+phi_in_step;
            phi_out = phi_out+phi_out_step;
        end
        % Update the reflection matrix
        r_update = exp(1i*phi_out).*r_z.*exp(1i*phi_in.');
        % Build the zone image
        I_zone = abs(finufft2d3(fy(:),fx(:),r_update(:),1,1e-2,Y(:),X(:))).^2;
        I_zone = reshape(I_zone,Ny,Nx);
        % Stitching window in x and y
        nb = 2*round(overlap2/dx_im); % Number of pixels overlapped
        windowx = ones(size(I_zone)); % The window to stitch in x
        windowy = ones(size(I_zone)); % The window to stitch in x
        if ~ismember(min(x_im,[],'all'),x) % If the zone is not on the left of the image
            windowx(:,1:nb) = (0:1/nb:1-1/nb).*ones(height(I_zone),nb); % In the overlapping area, the weight of the window declines linearly from 1 to 0
        end
        if ~ismember(max(x_im,[],'all'),x) % If the zone is not on the right of the image
            windowx(:,end-nb+1:end) = fliplr((0:1/nb:1-1/nb).*ones(height(I_zone),nb));
        end
        if ~ismember(min(y_im,[],'all'),y) % If the zone is not at the top of the image
            windowy(1:nb,:) = transpose(0:1/nb:1-1/nb).*ones(nb,width(I_zone));
        end
        if ~ismember(max(y_im,[],'all'),y) % If the zone is not at the bottom of the image
            windowy(end-nb+1:end,:) = flipud(transpose(0:1/nb:1-1/nb).*ones(nb,width(I_zone)));
        end
        I_zone = I_zone.*windowx.*windowy; % Apply the stitching window
        I_zone_expand = sparse(length(y_im),length(x_im));
        I_zone_expand(ny,nx) = I_zone; % Put the zone image to the corresponding position in the big image
        I = I+I_zone_expand;
    end
    
    % Save the image
    save(""+directory_save+"./I_2D_division_"+n_division+"_single_"+rad_order_start+"_"+rad_order_inc+".mat",'I')
end

