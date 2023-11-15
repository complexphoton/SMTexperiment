% The function to compensate for dispersion
function [phase_list,z_target] = dispersion_class(x,y,z_class,list_amp,directory_r,prefix,list_k0,k,im_case,directory_save)
    % 1. Setting
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
    [X,Y,Z_class] = meshgrid(single(x),single(y),single(z_class));
    Nx = length(x);
    Ny = length(y);
    Nz = length(z_class);
    
    % 2. Build matrices
    % Build the matrix Psi whose each column contains a single-frequency
    % image and matrix Omega that contains the frequency coefficients
    % (omega_i-omega_c)^j/omega_c
    Psi = [];
    Omega = [];
    for i_freq = 1:n_freq
        fprintf("Working with frequency"+i_freq+"\n")
        k0 = list_k0(i_freq); freq = 3e2*k0/2/pi;
        
        % Load the reflection matrix at this frequency
        r = single(load(""+directory_r+""+prefix+""+i_freq+".mat").r_pad);
        
        % In the sample
        kz_in = sqrt(k0^2-kx_in.^2-ky_in.^2);
        kz_out = -sqrt(k0^2-kx_out.^2-ky_out.^2);
        fz0 = kz_out-kz_in;
        
        % Build the image at this frequency and scale it with the laser
        % amplitude
        psi_s = finufft3d3(fy(:),fx(:),fz0(:),r(:),1,1e-2,Y(:),X(:),Z_class(:))/list_amp(i_freq);
        
        % Put it to the corresponding column of Psi
        Psi = [Psi psi_s];
        
        % Build the corresponding row of Omega
        Omega = vertcat(Omega,(freq-freqc)/freqc);
    end
    
    % 3. Scan the dispersion coefficients a1
    % First, find the normalization coefficients
    I_a = abs(sum(Psi,2)).^2; % Sum all the rows of Psi together, we obtain the image before dispersion compensation
    I0 = max(I_a,[],'all')/1000; % Normalization coefficient
    
    % Sweep a1 (longitudinal shift coefficient)
    list_M_a1 = []; % List of FOM with different a1

    for a1 = 0:10:4000
        phase_a1 = Omega(:,1)*a1;
        psi_a1 = Psi*exp(-1i*phase_a1); % Image with the corresponding a1
        I_a1 = abs(psi_a1).^2;
        M_a1 = sum(I_a1.*log(I_a1/I0),'all');
        list_M_a1 = [list_M_a1 [M_a1;a1]]; 
        clear psi_a1
    end
    
    figure(3)
    plot(list_M_a1(2,:),list_M_a1(1,:))
    % Find the coefficient that gives the highest sharpness
    a1_max = list_M_a1(2,list_M_a1(1,:)==max(list_M_a1(1,:)));
    
    % 4. Local optimization 
    clear opt
    a_init = double([a1_max]); % Initial guess from the scanned result
    opt.algorithm = NLOPT_LD_LBFGS;
    my_func = @(a) sharpness_figure_of_merit(double(a), double(Psi), -double(Omega));
    opt.min_objective = @(a) my_func(a);
    % Convergence criteria
    opt.ftol_rel = 1e-4;
    opt.xtol_rel = 1e-4;
    % Bounds, restricting the algorithm from driving the coefficients to
    % far away
    opt.lower_bounds = double(a1_max)-100;    
    opt.upper_bounds = double(a1_max)+100;  
    opt.maxeval = 200;
    opt.verbose = 1;
    % Run optimization with the initial guess
    [a, ~] = nlopt_optimize(opt, a_init);
    phase_list = Omega*a;
    
    % Save the phase list
    save(""+directory_save+"./phase_list_"+im_case+"_class.mat",'phase_list')
    
    % 5. Find the depth of the image if the imaging case is 2D
    if im_case == "2D"
        % Corrected image
        psi_cor = Psi*exp(-1i*Omega*a);
        psi_cor = reshape(psi_cor,Nx,Ny,Nz);
        I_cor = abs(psi_cor).^2;
        list_M = []; % list of FOM of 2D image at each depth
        for ii = 1:Nz % Scan each depth and find the FOM of the 2D images
            I_2D = I_cor(:,:,ii);
            M = sum(I_2D.*log(I_2D/I0),'all');
            list_M = [list_M M];    
        end
        % Search for target depth
        z_target = z_class(list_M == max(list_M,[],'all'));
        
        I_target = I_cor(:,:,list_M == max(list_M,[],'all'));
        
        figure(1)
        imagesc(fliplr(interp2(I_target,3)))
        axis image
        colormap('hot')
        
        figure(2)
        plot(z_class,list_M)
        
    elseif im_case == "3D"
        z_target = [];
    end

end
