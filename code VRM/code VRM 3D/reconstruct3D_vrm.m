% Function to reconstruct the fully corrected 2D SMT image
function [I] = reconstruct3D_vrm(list_z_im,x_im,y_im,dz_im,z_im,overlap2,k,list_k0,list_amp,directory_r,prefix,directory_save,coef_n2,coef_n1,z_ref_air,h1)
    
    n_sub = length(list_z_im); % Number of subvolumes
    % Input and output k:
    kx_out = k(:,1); ky_out = k(:,2);
    kx_in = k(:,1).'; ky_in = k(:,2).';
    % The Fourier components
    fx = kx_out-kx_in;
    fy = ky_out-ky_in;
    % Build the image of each subvolume
    Psi = single(zeros(length(x_im),length(y_im),length(z_im))); % Initialize the volumetric complex image
    n_freq = length(list_k0); % Number of frequency    
    freqc = 3e2*list_k0(round(n_freq/2))/2/pi;
    Nz_im = length(z_im);
    % Scan each frequency and build the volumetric image at each frequency
    for i_freq = 1:n_freq
        k0 = list_k0(i_freq); freq = 3e2*k0/2/pi; dfreq = freq - freqc; 
        
        % The refractive indices of glass and the medium at the
        % corresponding frequency
        n1 = coef_n1(1) + coef_n1(2)*dfreq + coef_n1(3)*dfreq.^2 + coef_n1(4)*dfreq.^3 + coef_n1(5)*dfreq.^4;
        n_media = coef_n2(1) + coef_n2(2)*dfreq + coef_n2(3)*dfreq.^2 + coef_n2(4)*dfreq.^3 + coef_n2(5)*dfreq.^4;
        
        % Load the reflection matrix at this frequency 
        r = single(load(""+directory_r+""+prefix+""+i_freq+".mat").r_pad)/list_amp(i_freq);
        
        % Find kz and Fourier components fz 
        % In air
        kz_in0 = (sqrt(k0^2-kx_in.^2-ky_in.^2));
        kz_out0 = (-sqrt(k0^2-kx_out.^2-ky_out.^2));
        fz0 = (kz_out0-kz_in0);
        % In glass
        kz_in1 = (sqrt(k0^2*n1^2-kx_in.^2-ky_in.^2));
        kz_out1 = (-sqrt(k0^2*n1^2-kx_out.^2-ky_out.^2));
        fz1 = (kz_out1-kz_in1);   
        % In the sample
        kz_in = (sqrt(k0^2*n_media^2-kx_in.^2-ky_in.^2));
        kz_out = (-sqrt(k0^2*n_media^2-kx_out.^2-ky_out.^2));
        fz = (kz_out-kz_in);
        
        % Correct index mismatch
        r = r.*exp(1i.*mod(-fz0*(z_ref_air+h1)+fz1*h1,2*pi));
  
        % Update the reflection matrix with the phase correction of the
        % corresponding subvolume and zone
        psi_s = single(zeros(length(x_im),length(y_im),length(z_im))); % Initialize the single-frequency volumetric complex image
        for id_sub = 1:n_sub
            fprintf("Build 3D image, frequency"+i_freq+", subvolume "+id_sub+"\n")
                        
            % Correct dispersion
            phi_in_dis = load(""+directory_r+"r_vrm_"+i_freq+"_subvolume_"+id_sub+".mat").phi_in;
            phi_out_dis = load(""+directory_r+"r_vrm_"+i_freq+"_subvolume_"+id_sub+".mat").phi_out;
            r_dis = exp(-1i*phi_out_dis).*r.*exp(-1i*phi_in_dis.');
        
            % z coordinate of the subvolume
            z_sub = list_z_im{id_sub,1};
            Nz_sub = length(z_sub); % Number of pixels in z of the zone
            % Build the single-frequency image of the subvolume
            Nx = length(x_im); 
            Ny = length(y_im);
            [X,Y,Z_sub] = meshgrid(x_im,y_im,z_sub); % Make grid points
            nz = ceil((z_sub-list_z_im{1,1}(1)+dz_im/2)/dz_im);
            
            % Wavefront correction
            phi_in = load(""+directory_save+"phi_in_subvolume_"+id_sub+"_vrm.mat").phi_in;
            phi_out = load(""+directory_save+"phi_out_subvolume_"+id_sub+"_vrm.mat").phi_out;
            r_wavefront = exp(1i*phi_out).*r_dis.*exp(1i*phi_in.');
            
            % Build the single-frequency zone image
            psi_subvolume = finufft3d3(fy(:),fx(:),fz(:),r_wavefront(:),1,1e-2,Y(:),X(:),Z_sub(:));
            psi_subvolume = reshape(psi_subvolume,Ny,Nx,Nz_sub);
            
            % Stitching window in z
            nb = 2*round(overlap2/dz_im); % Number of pixels overlapped
            windowz = single(ones(Ny,Nx,Nz_sub));
            for jj = 1:nb % Scan each depth in the overlapping region and assign the corresponding weight to that depth
                if id_sub ~= 1 % If it is not the top subvolumes (the volume with smallest z)
                    windowz(:,:,jj) = single(ones(Ny,Nx)*(jj-1)*1/nb).^0.75;
                end
                if id_sub ~= n_sub % If it is not the bottom subvolumes
                    windowz(:,:,end-jj+1) = single(ones(Ny,Nx)*(jj-1)*1/nb).^0.75;
                end
            end
            psi_subvolume = psi_subvolume.*windowz; % Apply the stitching window
            psi_expand = single(zeros(length(x_im),length(y_im),length(z_im)));
            psi_expand(:,:,nz) = psi_subvolume; % Put the zone image to its corresponding position in the big image
            psi_s = psi_s+psi_expand;
        end
        Psi = Psi+psi_s;
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

    save(""+directory_save+"./I_3D_vrm.mat",'I')
end
