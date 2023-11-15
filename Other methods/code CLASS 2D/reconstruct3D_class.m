% Function to reconstruct the fully corrected 2D SMT image
function [I] = reconstruct3D_class(list_x_im,list_y_im,list_z_class_im,x_im,y_im,z_class_im,overlap2,k,list_k0,phase_list,list_amp,n_division,directory_r,prefix,directory_save)
    
    n_sub = length(list_z_class_im); % Number of subvolumes
    dz_im = abs(z_class_im(2)-z_class_im(1)); % Pixel size in yz
    n_zone = 2^(n_division*2); % Number of zones
    dx_im = round(abs(x_im(2)-x_im(1)),1); % Pixel size in x,y
    % Input and output k:
    kout_x = k(:,1); kout_y = k(:,2);
    kin_x = k(:,1).'; kin_y = k(:,2).';
    N = length(k); % Number of k
    % The Fourier components
    fx = kout_x-kin_x;
    fy = kout_y-kin_y;
    % Build the image of each subvolume
    Psi = zeros(length(x_im),length(y_im),length(z_class_im)); % Initialize the volumetric complex image
    n_freq = length(list_k0); % Number of frequency    
    % Scan each frequency and build the volumetric image at each frequency
    for i_freq = 1:n_freq
        fprintf("Working with frequency"+i_freq+"\n")
        k0 = list_k0(i_freq);
        
        % Load the reflection matrix at this frequency
        r = single(load(""+directory_r+""+prefix+""+i_freq+".mat").r_pad);
        
        % Find kz and Fourier components fz 
        kz_in = sqrt(k0^2-kin_x.^2-kin_y.^2);
        kz_out = -sqrt(k0^2-kout_x.^2-kout_y.^2);
        fz0 = kz_out-kz_in;
        
        % Compensate dispersion and normalize with laser amplitude
        r = r.*exp(-1i*phase_list(i_freq))./list_amp(i_freq);
        
        % Update the reflection matrix with the phase correction of the
        % corresponding subvolume and zone
        psi_s = zeros(length(x_im),length(y_im),length(z_class_im)); % Initialize the single-frequency volumetric complex image
        for subvolume_id = 1:n_sub
            % z coordinate of the subvolume
            z = list_z_class_im{subvolume_id,1};
            Nz = length(z); % Number of pixels in z of the zone
            % Build the single-frequency image of each zone in the subvolume
            for zone_id = 1:n_zone
                % Zone coordinate
                x = list_x_im{zone_id,n_division+1};
                y = list_y_im{zone_id,n_division+1};
                Nx = length(x); % Number of pixels in x and y of the zone
                Ny = length(y);
                [X,Y,Z] = meshgrid(single(x),single(y),single(z)); % Make grid points
                nx = round((x+dx_im/2)/dx_im); % Zone coordinate in pixels
                ny = round((y+dx_im/2)/dx_im);
                nz = round((z-list_z_class_im{1,1}(1)+dz_im/2)/dz_im);
                phi_in = load(""+directory_save+"phi_in_"+n_division+"_zone_"+zone_id+"_subvolume_"+subvolume_id+"_class.mat").phi_in;
                phi_out = load(""+directory_save+"phi_out_"+n_division+"_zone_"+zone_id+"_subvolume_"+subvolume_id+"_class.mat").phi_out;
                % Update the reflection matrix
                r_zone = exp(1i*phi_out).*r.*exp(1i*phi_in.');
                % Build the single-frequency zone image
                psi_zone = finufft3d3(fy(:),fx(:),fz0(:),r_zone(:),1,1e-2,Y(:),X(:),Z(:));
                psi_zone = reshape(psi_zone,Ny,Nx,Nz);
                % Stitching window in x and y
                nb = 2*round(overlap2/dx_im); % Number of pixels overlapped
                windowx = ones(Ny,Nx,Nz); % The window to stitch in x
                windowy = ones(Ny,Nx,Nz); % The window to stitch in x
                if ~ismember(min(x_im,[],'all'),x) % If the zone is not on the left of the image
                    windowx(:,1:nb,:) = (0:1/nb:1-1/nb).*ones(Ny,nb,Nz); % In the overlapping area, the weight of the window declines linearly from 1 to 0
                end
                if ~ismember(max(x_im,[],'all'),x) % If the zone is not on the right of the image
                    windowx(:,end-nb+1:end,:) = fliplr((0:1/nb:1-1/nb).*ones(Ny,nb,Nz));
                end
                if ~ismember(min(y_im,[],'all'),y) % If the zone is not at the top of the image
                    windowy(1:nb,:,:) = transpose(0:1/nb:1-1/nb).*ones(nb,Nx,Nz);
                end
                if ~ismember(max(y_im,[],'all'),y) % If the zone is not at the bottom of the image
                    windowy(end-nb+1:end,:,:) = flipud(transpose(0:1/nb:1-1/nb).*ones(nb,Nx,Nz));
                end
                psi_zone = psi_zone.*windowx.*windowy; % Apply the stitching window
                % Stitching window in z
                windowz = ones(Ny,Nx,Nz);
                for jj = 1:nb % Scan each depth in the overlapping region and assign the corresponding weight to that depth
                    if ~ismember(min(z_class_im,[],'all'),z) % If it is not the top subvolumes (the volume with smallest z)
                        windowz(:,:,jj) = (ones(Ny,Nx)*(jj-1)*1/nb);
                    end
                    if ~ismember(max(z_class_im,[],'all'),z) % If it is not the bottom subvolumes
                        windowz(:,:,end-jj+1) = (ones(Ny,Nx)*(jj-1)*1/nb);
                    end
                end
                psi_zone = psi_zone.*windowz; % Apply the stitching window
                psi_expand = zeros(length(x_im),length(y_im),length(z_class_im));
                psi_expand(ny,nx,nz) = psi_zone; % Put the zone image to its corresponding position in the big image
                psi_s = psi_s+psi_expand;
            end
        end
        Psi = Psi+psi_s;
    end
    % Finally, the volumetric image is
    I = abs(Psi).^2;
    save(""+directory_save+"./I_3D_class.mat",'I')
end
