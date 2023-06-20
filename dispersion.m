% The function to compensate for dispersion
function [phase_list,z_target] = dispersion(x,y,z,list_amp,directory_r,prefix,h1,z_mirror,coef_n1,coef_n2,list_k0,k,im_case,directory_save)
    %% 1. Setting
    % The number of frequency
    n_freq = length(list_k0);
    % The central frequency
    freqc = 3e2*list_k0(round(n_freq/2))/2/pi;
    % Input and output k
    kx_out = k(:,1);
    ky_out = k(:,2);
    kx_in = k(:,1).';
    ky_in = k(:,2).';
    % The Fourier components
    fx = kx_out-kx_in;
    fy = ky_out-ky_in;
    % The coordinate is built into grid points
    [X,Y,Z] = meshgrid(x,y,z);
    Nx = length(x);
    Ny = length(y);
    Nz = length(z);
    
    %% 2. Build matrices
    % Build the matrix Psi whose each column contains a single-frequency
    % image and matrix Omega that contains the frequency coefficients
    % (omega_i-omega_c)^j/omega_c
    Psi = [];
    Omega = [];
    for i_freq = 1:n_freq
        fprintf("Working with frequency"+i_freq+"\n")
        k0 = list_k0(i_freq); freq = 3e2*k0/2/pi; dfreq = freq - freqc; 
        
        % The refractive indices of glass and the medium at the
        % corresponding frequency
        n1 = coef_n1(1) + coef_n1(2)*dfreq + coef_n1(3)*dfreq.^2 + coef_n1(4)*dfreq.^3 + coef_n1(5)*dfreq.^4;
        n_media = coef_n2(1) + coef_n2(2)*dfreq + coef_n2(3)*dfreq.^2 + coef_n2(4)*dfreq.^3 + coef_n2(5)*dfreq.^4;
        
        % Load the reflection matrix at this frequency
        r = load(""+directory_r+""+prefix+""+i_freq+".mat").r_pad;
        
        % Find kz and Fourier components fz 
        % In air
        kz_in0 = sqrt(k0^2-kx_in.^2-ky_in.^2);
        kz_out0 = -sqrt(k0^2-kx_out.^2-ky_out.^2);
        fz0 = kz_out0-kz_in0;
        % In glass
        kz_in1 = sqrt(k0^2*n1^2-kx_in.^2-ky_in.^2);
        kz_out1 = -sqrt(k0^2*n1^2-kx_out.^2-ky_out.^2);
        fz1 = kz_out1-kz_in1;   
        % In the sample
        kz_in = sqrt(k0^2*n_media^2-kx_in.^2-ky_in.^2);
        kz_out = -sqrt(k0^2*n_media^2-kx_out.^2-ky_out.^2);
        fz = kz_out-kz_in;
        
        % Correct index mismatch
        r = r(:).*exp(1i.*(-fz0(:)*z_mirror+fz1(:)*h1));
        
        % Build the image at this frequency and scale it with the laser
        % amplitude
        psi_s = finufft3d3(fy(:),fx(:),fz(:),r(:),1,1e-2,Y(:),X(:),Z(:))/list_amp(i_freq);
        % Put it to the corresponding column of Psi
        Psi = [Psi psi_s];
        % Build the corresponding row of Omega
        Omega = vertcat(Omega,[(freq-freqc)/freqc ((freq-freqc)/freqc)^2 ((freq-freqc)/freqc)^3]);
    end
    
    %% 3. Scan the dispersion coefficients a = [a_1; a_2; a_3]
    % First, find the normalization coefficients
    I_a = abs(sum(Psi,2)).^2; % Sum all the rows of Psi together, we obtain the image before dispersion compensation
    I0 = max(I_a,[],'all')/1000; % Normalization coefficient
    
    % Sweep a1 (longitudinal shift coefficient)
    list_M_a1 = []; % List of FOM with different a1

    for a1 = -1000:10:1000
        phase_a1 = Omega(:,1)*a1;
        psi_a1 = Psi*exp(-1i*phase_a1); % Image with the corresponding a1
        I_a1 = abs(psi_a1).^2;
        M_a1 = sum(I_a1.*log(I_a1/I0),'all');
        list_M_a1 = [list_M_a1 [M_a1;a1]]; 
        clear psi_a1
    end
    % Find the coefficient that gives the highest sharpness
    a1_max = list_M_a1(2,list_M_a1(1,:)==max(list_M_a1(1,:)));

    % Sweep a2 (controlling width of the pulse)
    list_M_a2 = []; % List of FOM with different a2

    for a2 = -1000:10:1000
        phase_a2 = Omega(:,1:2)*[a1_max; a2];
        psi_a2 = Psi*exp(-1i*phase_a2); % Image with the corresponding a2 and a1 fixed at a1_max
        I_a2 = abs(psi_a2).^2;
        M_a2 = sum(I_a2.*log(I_a2/I0),'all');
        list_M_a2 = [list_M_a2 [M_a2;a2]]; 
        clear psi_a2
    end
    % Find the coefficient that gives the highest sharpness
    a2_max = list_M_a2(2,list_M_a2(1,:)==max(list_M_a2(1,:)));

    % Sweep a3 (asymmetric pulse deformation)
    list_M_a3 = []; % List of FOM with different a3

    for a3 = -1000:10:1000
        phase_a3 = Omega*[a1_max; a2_max; a3];
        psi_a3 = (Psi./list_amp)*exp(-1i*phase_a3); % Image with the corresponding a3 and a1, a2 fixed at a1_max, a2_max
        I_a3 = abs(psi_a3).^2;
        M_a3 = sum(I_a3.*log(I_a3/I0),'all');
        list_M_a3 = [list_M_a3 [M_a3;a3]]; 
        clear psi_a3
    end
    % Find the coefficient that gives the highest sharpness
    a3_max = list_M_a3(2,list_M_a3(1,:)==max(list_M_a3(1,:)));
    
    %% 4. Local optimization 
    clear opt
    a_init = [a1_max;a2_max;a3_max]; % Initial guess from the scanned result
    opt.algorithm = NLOPT_LD_LBFGS;
    my_func = @(a) sharpness_figure_of_merit(a, Psi, -Omega);
    opt.min_objective = @(a) my_func(a);
    % Convergence criteria
    opt.ftol_rel = 1e-4;
    opt.xtol_rel = 1e-4;
    % Bounds, restricting the algorithm from driving the coefficients to
    % far away
    opt.lower_bounds = [a1_max-100 a2_max-100 a3_max-100].';    
    opt.upper_bounds = [a1_max+100 a2_max+100 a3_max+100].';  
    opt.maxeval = 200;
    opt.verbose = 1;
    % Run optimization with the initial guess
    [a, ~] = nlopt_optimize(opt, a_init);
    phase_list = Omega*a;
    
    % Save the phase list
    save(""+directory_save+"./phase_list_"+im_case+".mat",'phase_list')
    
    %% 5. Find the depth of the image if the imaging case is 2D
    if im_case == "2D"
        % Corrected image
        psi_cor = Psi*exp(-1i*phase_list);
        psi_cor = reshape(psi_cor,Nx,Ny,Nz);
        list_M = []; % list of FOM of 2D image at each depth
        for ii = 1:Nz % Scan each depth and find the FOM of the 2D images
            I_2D = abs(psi_cor(:,:,ii).^2);
            M = sum(I_2D.*log(I_2D/I0),'all');
            list_M = [list_M M];    
        end
        % Search for target depth
        z_target = z(list_M == max(list_M,[],'all'));
    elseif im_case == "3D"
        z_target = [];
    end
    clear Psi Omega psi_s psi_cor
end
