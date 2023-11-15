% Function to reconstruct the fully corrected 2D SMT image
function [I] = reconstruct3D_vrm(list_z_vrm_im,x_im,y_im,z_vrm_im,overlap2,k,list_k0,phase_list,list_amp,directory_r,prefix,directory_save)
    n_sub = length(list_z_vrm_im); % Number of subvolumes
    dz_im = round(abs(z_vrm_im(2)-z_vrm_im(1)),1); % Pixel size in yz
    % Input and output k:
    kout_x = k(:,1); kout_y = k(:,2);
    kin_x = k(:,1).'; kin_y = k(:,2).';
    N = length(k); % Number of k
    % The Fourier components
    fx = kout_x-kin_x;
    fy = kout_y-kin_y;
    % Build the image of each subvolume
    Psi = zeros(length(x_im),length(y_im),length(z_vrm_im)); % Initialize the volumetric complex image
    n_freq = length(list_k0); % Number of frequency    
    % Scan each frequency and build the volumetric image at each frequency
    for i_freq = 1:n_freq
        fprintf("Working with frequency"+i_freq+"\n")
        k0 = list_k0(i_freq);
        
        % Load the reflection matrix at this frequency
        r = load(""+directory_r+""+prefix+""+i_freq+".mat").r_pad;
        
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
            z = list_z_vrm_im{subvolume_id,1};
            Nz = length(z); % Number of pixels in z of the zone
            % Build the single-frequency image of the subvolume
            Nx = length(x); 
            Ny = length(y);
            [X,Y,Z] = meshgrid(x_im,y_im,z); % Make grid points
            nz = round((z-list_z_vrm_im{1,1}(1)+dz_im/2)/dz_im);
            phi_in = load(""+directory_save+"phi_in_"+division_step+"_zone_"+zone+"_subvolume_"+subvolume_id+"_class.mat").phi_in;
            phi_out = load(""+directory_save+"phi_out_"+division_step+"_zone_"+zone+"_subvolume_"+subvolume_id+"_class.mat").phi_out;
            % Update the reflection matrix
            r_update = exp(1i*phi_out).*r.*exp(1i*phi_in.');
            % Build the single-frequency zone image
            psi_subvolume = finufft3d3(fy(:),fx(:),fz0(:),r_update(:),1,1e-2,Y(:),X(:),Z(:));
            psi_subvolume = reshape(psi_subvolume,Ny,Nx,Nz);
            % Stitching window in z
            nb = 2*round(overlap2/dx_im); % Number of pixels overlapped
            windowz = ones(Ny,Nx,Nz);
            for jj = 1:nb % Scan each depth in the overlapping region and assign the corresponding weight to that depth
                if ~ismember(min(z_im,[],'all'),z) % If it is not the top subvolumes (the volume with smallest z)
                    windowz(:,:,jj) = (ones(Ny,Nx)*(jj-1)*1/nb);
                end
                if ~ismember(max(z_im,[],'all'),z) % If it is not the bottom subvolumes
                    windowz(:,:,end-jj+1) = (ones(Ny,Nx)*(jj-1)*1/nb);
                end
            end
            psi_subvolume = psi_subvolume.*windowz; % Apply the stitching window
            psi_expand = zeros(length(x_im),length(y_im),length(z_class_im));
            psi_expand(:,:,nz) = psi_subvolume; % Put the zone image to its corresponding position in the big image
            psi_s = psi_s+psi_expand;
        end
        Psi = Psi+psi_s;
    end
    % Finally, the volumetric image is
    I = abs(Psi).^2;
    save(""+directory_save+"./I_3D_vrm.mat",'I')
end