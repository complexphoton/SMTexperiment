% Function to reconstruct the fully corrected 2D SMT image
function [I] = reconstruct3D(list_x_im,list_y_im,list_z_im,x_im,y_im,z_im,overlap2xy,overlap2z,k,list_k0,phase_list,list_amp,coef_n1,coef_n2,h1,z_mirror,n_division,directory_r,prefix,directory_save)
   
    n_sub = length(list_z_im); % Number of subvolumes
    dz_im = round(abs(z_im(2)-z_im(1)),1); % Pixel size in yz
    n_zone = 2^(n_division*2); % Number of zones
    dx_im = round(abs(x_im(2)-x_im(1)),1); % Pixel size in x,y
    Nz_im = length(z_im);
    % Input and output k:
    kout_x = k(:,1); kout_y = k(:,2);
    kin_x = k(:,1).'; kin_y = k(:,2).';
    N = length(k); % Number of k
    % The Fourier components
    fx = kout_x-kin_x;
    fy = kout_y-kin_y;
    % Build the image of each subvolume
    Psi = single(zeros(length(x_im),length(y_im),length(z_im))); % Initialize the volumetric complex image
    n_freq = length(list_k0); % Number of frequency    
    freqc = 3e2*list_k0(round(n_freq/2))/2/pi; % The central frequency
    % Scan each frequency and build the volumetric image at each frequency
    for i_freq = [1:n_freq]
        tic
        k0 = list_k0(i_freq); freq = 3e2*k0/2/pi; dfreq = freq - freqc; 
        
        % The refractive indices of glass and the medium at the
        % corresponding frequency
        n1 = coef_n1(1) + coef_n1(2)*dfreq + coef_n1(3)*dfreq.^2 + coef_n1(4)*dfreq.^3 + coef_n1(5)*dfreq.^4;
        n_media = coef_n2(1) + coef_n2(2)*dfreq + coef_n2(3)*dfreq.^2 + coef_n2(4)*dfreq.^3 + coef_n2(5)*dfreq.^4;
        
        % Load the reflection matrix at this frequency
        r = load(""+directory_r+""+prefix+""+i_freq+".mat").r_pad;
        
        % Find kz and Fourier components fz 
        % In air
        kz_in0 = sqrt(k0^2-kin_x.^2-kin_y.^2);
        kz_out0 = -sqrt(k0^2-kout_x.^2-kout_y.^2);
        fz0 = kz_out0-kz_in0;
        % In glass
        kz_in1 = sqrt(k0^2*n1^2-kin_x.^2-kin_y.^2);
        kz_out1 = -sqrt(k0^2*n1^2-kout_x.^2-kout_y.^2);
        fz1 = kz_out1-kz_in1;   
        % In the sample
        kz_in = sqrt(k0^2*n_media^2-kin_x.^2-kin_y.^2);
        kz_out = -sqrt(k0^2*n_media^2-kout_x.^2-kout_y.^2);
        fz = kz_out-kz_in;
        
        % Correct index mismatch
        r = r(:).*exp(1i.*mod(-fz0(:)*z_mirror+fz1(:)*h1,2*pi));
        
        % Compensate dispersion and normalize with laser amplitude
        r = r.*exp(-1i*wrapToPi(phase_list(i_freq)))./list_amp(i_freq);
        r = reshape(r,N,N);
        
        % Update the reflection matrix with the phase correction of the
        % corresponding subvolume and zone
        psi_s = zeros(length(x_im),length(y_im),length(z_im)); % Initialize the single-frequency volumetric complex image
        for subvolume_id = 1:n_sub
            % z coordinate of the subvolume
            z_sub = list_z_im{subvolume_id,1};
            Nz_sub = length(z_sub); % Number of pixels in z of the zone
            % Build the single-frequency image of each zone in the subvolume
            for zone_id = 1:n_zone
                
                fprintf("Working with frequency"+i_freq+", subvolume "+subvolume_id+", zone "+zone_id+"\n")
                
                % Zone coordinate
                x_zone = list_x_im{zone_id,n_division+1};
                y_zone = list_y_im{zone_id,n_division+1};
                Nx_im_zone = length(x_zone); % Number of pixels in x and y of the zone
                Ny_im_zone = length(y_zone);
                [X_zone,Y_zone,Z_zone] = meshgrid(x_zone,y_zone,z_sub); % Make grid points
                nx = round((x_zone+dx_im/2)/dx_im); % Zone coordinate in pixels
                ny = round((y_zone+dx_im/2)/dx_im);
                nz = ceil((z_sub-list_z_im{1,1}(1)+dz_im/2)/dz_im);
                % Find the final aberration phase of the zone by summing the phase
                % at each division step
                phi_in = zeros(N,1); % Initialize the input and output aberration phase
                phi_out = zeros(N,1);
                for division_step = n_division:-1:0
                    if division_step == n_division
                        zone = zone_id; 
                    else
                        zone = ceil(zone/4);
                    end
                    c_in_step = load(""+directory_save+"c_in_"+division_step+"_zone_"+zone+"_subvolume_"+subvolume_id+"_out_of_"+n_sub+".mat").c_in;
                    c_out_step = load(""+directory_save+"c_out_"+division_step+"_zone_"+zone+"_subvolume_"+subvolume_id+"_out_of_"+n_sub+".mat").c_out;
                    Z_step = single(load(""+directory_save+"Z_"+division_step+"_re.mat").Z_re);
                    phi_in_step = Z_step*c_in_step;
                    phi_out_step = Z_step*c_out_step;
                    phi_in = phi_in+phi_in_step;
                    phi_out = phi_out+phi_out_step;
                end
                % Update the reflection matrix
                r_zone = exp(1i*phi_out).*r.*exp(1i*phi_in.');
                % Build the single-frequency zone image
                psi_zone = finufft3d3(fy(:),fx(:),fz(:),r_zone(:),1,1e-2,Y_zone(:),X_zone(:),Z_zone(:));
                psi_zone = reshape(psi_zone,Ny_im_zone,Nx_im_zone,Nz_sub);
                % Stitching window in x and y
                nbxy = 2*round(overlap2xy/dx_im); % Number of pixels overlapped
                windowx = ones(Ny_im_zone,Nx_im_zone,Nz_sub); % The window to stitch in x
                windowy = ones(Ny_im_zone,Nx_im_zone,Nz_sub); % The window to stitch in x
                if ~ismember(min(x_im,[],'all'),x_zone) % If the zone is not on the left of the image
                    windowx(:,1:nbxy,:) = (0:1/nbxy:1-1/nbxy).*ones(Ny_im_zone,nbxy,Nz_sub); % In the overlapping area, the weight of the window declines linearly from 1 to 0
                end
                if ~ismember(max(x_im,[],'all'),x_zone) % If the zone is not on the right of the image
                    windowx(:,end-nbxy+1:end,:) = fliplr((0:1/nbxy:1-1/nbxy).*ones(Ny_im_zone,nbxy,Nz_sub));
                end
                if ~ismember(min(y_im,[],'all'),y_zone) % If the zone is not at the top of the image
                    windowy(1:nbxy,:,:) = transpose(0:1/nbxy:1-1/nbxy).*ones(nbxy,Nx_im_zone,Nz_sub);
                end
                if ~ismember(max(y_im,[],'all'),y_zone) % If the zone is not at the bottom of the image
                    windowy(end-nbxy+1:end,:,:) = flipud(transpose(0:1/nbxy:1-1/nbxy).*ones(nbxy,Nx_im_zone,Nz_sub));
                end
                psi_zone = psi_zone.*windowx.*windowy; % Apply the stitching window
                % Stitching window in z
                nbz = 2*round(overlap2z/dz_im);
                windowz = single(ones(Ny_im_zone,Nx_im_zone,Nz_sub));
                for jj = 1:nbz % Scan each depth in the overlapping region and assign the corresponding weight to that depth
                    if subvolume_id ~= 1 % If it is not the top subvolumes (the volume with smallest z)
                        windowz(:,:,jj) = single((ones(Ny_im_zone,Nx_im_zone)*(jj-1)*1/nbz)).^0.75;
                    end
                    if subvolume_id ~= n_sub % If it is not the bottom subvolumes
                        windowz(:,:,end-jj+1) = single((ones(Ny_im_zone,Nx_im_zone)*(jj-1)*1/nbz)).^0.75;
                    end
                end
                psi_zone = psi_zone.*windowz; % Apply the stitching window
                psi_expand = zeros(length(x_im),length(y_im),length(z_im));
                psi_expand(ny,nx,nz) = psi_zone; % Put the zone image to its corresponding position in the big image
                psi_s = psi_s+psi_expand;
            end
        end
        Psi = Psi+psi_s;
        toc

    end
    % Finally, the volumetric image is
    I = abs(Psi).^2;
    
    for ii = 1:Nz_im
        I_2D = I(:,:,ii);
        M = sum(I_2D,'all');
        if M == 0
            if ii ~= 1 & ii ~= Nz_im
                I(:,:,ii) = 1/2*(I(:,:,ii+1)+I(:,:,ii-1));
            end
        end
    end
    
    save(""+directory_save+"./I_3D_with_overlap_"+n_sub+"_subvolumes.mat",'I')
end
