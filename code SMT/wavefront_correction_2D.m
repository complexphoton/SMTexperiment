% Function to correct wavefront
function wavefront_correction_2D(x,y,r_z,division_step,Z,zone_id,subvolume_id,directory_save)
        
    fprintf("Correcting zone "+zone_id+" division step "+division_step+" subvolume "+subvolume_id+".\n");
    % Inputs: 
    % x,y spatial coordinates of the zone
    % r_z: updated time-gated reflection matrix of the bigger zone containing this zone
    % Z: Zernike matrix of this zone
    Nx = length(x);
    Ny = length(y);
    % Build grid points 
    [X,Y] = (meshgrid(x,y));
    X = single(X(:));
    Y = single(Y(:));
    
    % 1. Load k, do basis conversion
    % 1.1. Load k
    % If we are correcting a smaller zone, we have to load the k of the
    % bigger zone
    if division_step > 0
        prev_step = division_step-1;
        k_big = single(load(""+directory_save+"k_"+prev_step+".mat").k);
        N_big = length(k_big);
        % Input and output k are
        kout_x_big = k_big(:,1); kout_y_big = k_big(:,2);
        kin_x_big = k_big(:,1).'; kin_y_big = k_big(:,2).';
    end
    % Load the k of the current zone
    k =  single(load(""+directory_save+"k_"+division_step+".mat").k);
    N = length(k);
    % Input and output k are:
    kout_x = k(:,1); kout_y = k(:,2);
    kin_x = k(:,1).'; kin_y = k(:,2).';
    % The Fourier components
    fx = single(kout_x-kin_x);
    fy = single(kout_y-kin_y);
    % 1.2. Do basis conversion
    if division_step > 0
        % 1.2.1. Convert to truncated spatial basis
        % Convert the input to spatial basis
        r_kout_r = zeros(N_big,Nx*Ny); % The reflection matrix whose input is spatial basis, output is angular basis
        parfor ii = 1:N_big
            r_kout_r(ii,:) = finufft2d3(kin_y_big.',kin_x_big.',r_z(ii,:).',-1,1e-2,Y,X).';
        end
        r_kout_r = single(r_kout_r/N_big);
        % Convert the output to spatial basis
        r_r = zeros(Nx*Ny,Nx*Ny); % The matrix whose both input and output are in spatial basis
        parfor ii = 1:Nx*Ny
            r_r(:,ii) = finufft2d3(-kout_y_big,-kout_x_big,r_kout_r(:,ii),-1,1e-2,Y,X);
        end
        r_r = single(r_r);
        clear r_kout_r
        % 1.2.2. Convert back to angular basis
        % Convert the input to angular basis
        r_r_kin = zeros(Nx*Ny,N); % The matrix whose input is in angular basis, output is spatial basis
        parfor ii = 1:Nx*Ny
            r_r_kin(ii,:) = finufft2d3(Y,X,r_r(ii,:).',1,1e-2,kin_y.',kin_x.').';
        end
        r_r_kin = single(r_r_kin/(Nx*Ny));
        clear r_r
        % Convert the output to angular basis
        r_k = zeros(N,N); % The matrix in angular basis
        parfor ii = 1:N
            r_k(:,ii) = finufft2d3(Y,X,r_r_kin(:,ii),1,1e-2,-kout_y,-kout_x);
        end
        r_k = single(r_k);
        clear r_r_kin
    elseif division_step == 0
        r_k = r_z;
        clear r_z
    end
    
    % 2. The initial image
    psi_in =  finufft2d3(fy(:),fx(:),r_k(:),1,1e-2,Y,X);
    I_in = abs(psi_in).^2;
    % Normalize to avoid the FOM getting too large
    I0 = max(I_in,[],'all')/1000; % Normalization constant
    I_in = I_in/I0; 
    
    % FOM of the initial image
    global M 
    M = -sum(I_in.*log(I_in),'all');
    
    % 3. Initialize the Zernike coefficients and other necessary matrices
    n_order = width(Z); % The number of Zernike polynomials
    c_in =  (zeros(n_order,1)); % Input and output Zernike coefficients
    c_out =  (zeros(n_order,1)); % Only this thing need to be double, the rest can be single.
    % The input and output-updated matrices
    r_in_update = r_k;
    r_out_update = r_k.';
    % To be more convenient when working with the input and output-updated
    % matrices, we define
    r_in = r_k; r_out = r_k.';
    
    % 4. Optimization
    delta_M = 1; % The change of FOM after each optimization loop
    while delta_M > 0.05
        M_prev = M; % Update the FOM of the previous loop
        
        % 4.1. Optimize input
        % 4.1.1. Build matrix Psi_in containing the image from each input
        Psi_in = zeros(Nx*Ny,N);
        parfor ii = 1:N
            Psi_in(:,ii) = finufft2d3(fy(:,ii),fx(:,ii),r_in_update(:,ii),1,1e-2,Y,X);
        end
        Psi_in =  (Psi_in/sqrt(I0)); % Normalize

        % 4.1.2. Optimize
        opt.algorithm =  NLOPT_LD_LBFGS; % Choose LBFGS algorithm
        my_func = @(c_in) sharpness_figure_of_merit(c_in,Psi_in,Z);
        opt.min_objective = @(c_in) my_func(c_in);
        % Convergence criteria of the optimizaton step
        opt.ftol_rel =  1e-4;        
        opt.xtol_rel =  1e-4;   
        opt.maxeval =  500; % The optimization will stop after a number of evaluation
        opt.verbose =  1;
        % The lower and upper bounds of the Zernike coefficients. This
        % bound is imposed to prevent the algorithm froom drive the Zernike
        % coefficients to very high values, which is not true in reality.
        % For us, the bound somewhere between +-1 and +-10 works fine.
        opt.lower_bounds =  -5*ones(n_order,1);
        opt.upper_bounds =  5*ones(n_order,1);
        % Initial guess
        c_in_init = double(c_in);
        % Run optimization
        [c_in,~] = nlopt_optimize(opt,c_in_init);
       
        % 4.1.3. Update the output-updated matrix
        phi_in = Z*c_in;
        r_out_update = exp(1i*phi_in).*r_out;

        % 4.2. Optimize output
        % 4.2.1. Build matrix Psi_out containing the image from each output
        Psi_out = zeros(Nx*Ny,N);
        parfor ii = 1:N
            Psi_out(:,ii) = (finufft2d3(fy(ii,:).',fx(ii,:).',r_out_update(:,ii),1,1e-2,Y,X));
        end
        Psi_out =  (Psi_out/sqrt(I0)); % Normalize
      
        % 4.2.2. Optimize
        my_func = @(c_out) sharpness_figure_of_merit(c_out,Psi_out,Z);
        opt.min_objective = @(c_out) my_func(c_out);
        % Initial guess
        c_out_init = c_out;
        % Run optimization
        [c_out,~] = nlopt_optimize(opt,c_out_init);
        
        phi_out = Z*c_out;
        r_in_update = exp(1i*phi_out).*r_in;
        
        % 4.3. Update the FOM change
        delta_M = (M_prev-M)/abs(M_prev); 
    end
    
    % 5. Update the reflection matrix of the zone, save the matrix and the Zernike coefficients
    phi_in = Z*c_in; phi_out = Z*c_out;
    r_update =  (exp(1i*phi_out).*r_k.*exp(1i*phi_in.'));   
   
    psi =  (finufft2d3(fy(:),fx(:),r_update(:),1,1e-2,Y,X));
    I = abs(psi).^2;
    I_zone = reshape(I,Ny,Nx);

    figure(1)
    imagesc(fliplr(interp2(I_zone,2.5)))
    axis image
    colormap('hot')
    
    save(""+directory_save+"c_in_"+division_step+"_zone_"+zone_id+"_subvolume_"+subvolume_id+".mat",'c_in');
    save(""+directory_save+"c_out_"+division_step+"_zone_"+zone_id+"_subvolume_"+subvolume_id+".mat",'c_out');
    save(""+directory_save+"r_update_"+division_step+"_zone_"+zone_id+"_subvolume_"+subvolume_id+".mat",'r_update');

end
