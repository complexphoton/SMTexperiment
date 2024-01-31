function dispersion_2nd_time(list_z,list_x,list_y,dx,overlap2xy,n_division,k,list_k0,coef_n1,coef_n2,directory_save,directory_r,prefix,z_mirror,h1)
    % Load the volume dispersion phase
    phase_list = load(""+directory_save+"./phase_list_"+im_case+"_single.mat").phase_list;
    % The function to compensate for dispersion
    n_sub = length(list_z); % Number of subvolumes
    n_zone = 2^(n_division*2); % Number of zones
    % Input and output k:
    kout_x = k(:,1); kout_y = k(:,2);
    kin_x = k(:,1).'; kin_y = k(:,2).';
    N = length(k); % Number of k
    % The Fourier components
    fx = kout_x-kin_x;
    fy = kout_y-kin_y;
    % Build the image of each subvolume
    n_freq = length(list_k0); % Number of frequency    
    freqc = 3e2*list_k0(round(n_freq/2))/2/pi; % The central frequency

    % 2. 2nd dispersion for each subvolume
    for subvolume_id = 1:n_sub
        % 2.1. Precompute
        Omega = single(zeros(n_freq,2));
        z_sub = list_z{subvolume_id,1}; Nz_sub = length(z_sub); % Number of pixels in z of the zone
        Psi = single(zeros(length(x)*length(y)*length(z_sub),n_freq));
        for i_freq = [1:n_freq]
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
            psi_s = zeros(length(x),length(y),length(z_sub)); % Initialize the single-frequency volumetric complex image

            % Build the single-frequency image of each zone in the subvolume
            for zone_id = 1:n_zone

                fprintf("Working with frequency"+i_freq+", subvolume "+subvolume_id+", zone "+zone_id+"\n")

                % Zone coordinate
                x_zone = list_x{zone_id,n_division+1};
                y_zone = list_y{zone_id,n_division+1};
                Nx_zone = length(x_zone); % Number of pixels in x and y of the zone
                Ny_zone = length(y_zone);
                [X_zone,Y_zone,Z_zone] = meshgrid(x_zone,y_zone,z_sub); % Make grid points
                nx = round((x_zone+dx/2)/dx); % Zone coordinate in pixels
                ny = round((y_zone+dx/2)/dx);
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
                psi_zone = reshape(psi_zone,Ny_zone,Nx_zone,Nz_sub);
                % Stitching window in x and y
                nbxy = 2*round(overlap2xy/dx); % Number of pixels overlapped
                windowx = ones(Ny_zone,Nx_zone,Nz_sub); % The window to stitch in x
                windowy = ones(Ny_zone,Nx_zone,Nz_sub); % The window to stitch in x
                if ~ismember(min(x_im,[],'all'),x_zone) % If the zone is not on the left of the image
                    windowx(:,1:nbxy,:) = (0:1/nbxy:1-1/nbxy).*ones(Ny_zone,nbxy,Nz_sub); % In the overlapping area, the weight of the window declines linearly from 1 to 0
                end
                if ~ismember(max(x_im,[],'all'),x_zone) % If the zone is not on the right of the image
                    windowx(:,end-nbxy+1:end,:) = fliplr((0:1/nbxy:1-1/nbxy).*ones(Ny_zone,nbxy,Nz_sub));
                end
                if ~ismember(min(y_im,[],'all'),y_zone) % If the zone is not at the top of the image
                    windowy(1:nbxy,:,:) = transpose(0:1/nbxy:1-1/nbxy).*ones(nbxy,Nx_zone,Nz_sub);
                end
                if ~ismember(max(y_im,[],'all'),y_zone) % If the zone is not at the bottom of the image
                    windowy(end-nbxy+1:end,:,:) = flipud(transpose(0:1/nbxy:1-1/nbxy).*ones(nbxy,Nx_zone,Nz_sub));
                end
                psi_zone = psi_zone.*windowx.*windowy; % Apply the stitching window
                psi_expand = zeros(length(x),length(y),length(z_sub));
                psi_expand(ny,nx,:) = psi_zone; % Put the zone image to its corresponding position in the big image
                psi_s = psi_s+psi_expand;
            end

            Psi(:,i_freq) = psi_s(:);
            Omega(i_freq,:) = [((freq-freqc)/freqc)^2 ((freq-freqc)/freqc)^3];
        end

        % 2.2. Do a serial scanning? 
        % Sweep a2
        list_M_a2 = [];
        for a2 = -1500:5:1500
            phase_a2 = Omega(:,1)*a2;
            I_a2 = single(abs(Psi*exp(-1i*phase_a2)).^2); % Image with the corresponding a2 and a1 fixed at a1_max
            M_a2 = sum(I_a2,'all');
            list_M_a2 = [list_M_a2 [M_a2;a2]]; 
            figure(11)
            plot(list_M_a2(2,:),list_M_a2(1,:))
        end
        a2_max_sub = list_M_a2(2,list_M_a2(1,:)==max(list_M_a2(1,:)));
        a2_max_sub = a2_max_sub(1);
        % Sweep a3
        list_M_a3 = [];
        for a3 = -1500:5:1500
            phase_a3 = Omega(:,1:2)*[a2_max_sub;a3];
            I_a3 = single(abs(Psi*exp(-1i*phase_a3)).^2); % Image with the corresponding a2 and a1 fixed at a1_max
            M_a3 = sum(I_a3,'all');
            list_M_a3 = [list_M_a3 [M_a3;a3]]; 
            figure(12)
            plot(list_M_a3(2,:),list_M_a3(1,:))
        end
        a3_max_sub = list_M_a3(2,list_M_a3(1,:)==max(list_M_a3(1,:)));
        a3_max_sub = a3_max_sub(1);

        delta_phase_list_2_sub = Omega*[a2_max_sub;a3_max_sub];
        phase_list_2 = phase_list+delta_phase_list_2_sub;

        % Save the phase list
        save(""+directory_save+"./phase_list_"+im_case+"_2nd_time_subvolume_"+subvolume_id+".mat",'phase_list_2')

        % Do a 2D scanning?
        a2_list = single(-100:5:100);
        a3_list = single(-1000:10:00);
        [A2,A3] = meshgrid(a2_list,a3_list);
        A2 = A2(:); A3 = A3(:);
        N_a = length(A2(:));
        list_M_A = [];
        for ii = 1:N_a
            phase_A = Omega*[A2(ii);A3(ii)];
            I_A = single(abs(Psi*exp(-1i*phase_A)).^2); % Image with the corresponding a2 and a1 fixed at a1_max
            M_A = sum((I_A/10e11).^2,'all');
            list_M_A = [list_M_A [M_A;ii]]; 
            figure(13)
            plot(list_M_A(2,:),list_M_A(1,:))
        end
        ii_max = list_M_A(2,list_M_A(1,:) == max(list_M_A(1,:)));
        ii_max = ii_max(1);

        phase_A = Omega*[A2(ii_max);A3(ii_max)];
        phase_list_2 = phase_list+phase_A;

        save(""+directory_save+"./phase_list_"+im_case+"_2nd_time_subvolume_"+subvolume_id+"_2d_scan.mat",'phase_list_2')

    end
 end

