% Function to correct wavefront
function wavefront_correction_vrm_2D(x,y,k,r_z,directory_save)
        
    % Inputs: 
    % x,y spatial coordinates of the zone
    % r_z: updated time-gated reflection matrix of the bigger zone containing this zone
    Nx = length(x); Ny = length(y);
    % Build grid points 
    [X,Y] = meshgrid(single(x),single(y));
    X = X(:); Y = Y(:);
    
    % 1. Do basis conversion
    % Input and output k are:
    N = length(k); 
    kout_x = k(:,1); kout_y = k(:,2);
    kin_x = k(:,1).'; kin_y = k(:,2).';
    % The Fourier components
    fx = kout_x-kin_x;
    fy = kout_y-kin_y;
    
    % 2. The initial image
    psi_in = finufft2d3(fy(:),fx(:),r_z(:),1,1e-2,Y,X);
    I_in = abs(psi_in).^2;
    
    % Normalize to avoid the FOM getting too large
    I0 = 1; % Normalization constant
    I_in = I_in/I0; 
    
    % FOM of the initial image
    global M
    M = sum(I_in,'all');
    
    I_in = reshape(I_in,Ny,Nx);
    
    figure(1)
    imagesc(fliplr(interp2(I_in,3)))
    axis image
    colormap('hot')
    
    % 3. Initialize
    phi_in_tot = zeros(N,1); phi_out_tot = zeros(N,1); 
    % The input and output-updated matrices
    r_in_update = r_z;
    r_out_update = r_z.';
    % To be more convenient when working with the input and output-updated
    % matrices, we define
    fx_in = fx; fy_in = fy;
    fx_out = fx.'; fy_out = fy.';
    
    % 4. Optimization
    max_phi = pi; % The maximum value of |phi_in| and |phi_out| in each iteration
    while max_phi > pi/18
        M_prev = M; % Update the FOM of the previous loop
        
        % 4.1. Optimize input
        % 4.1.1. Build matrix Psi_in containing the image from each input
        Psi_in = single(zeros(Nx*Ny,N));
        for ii = 1:N
            Psi_in(:,ii) = finufft2d3(fy_in(:,ii),fx_in(:,ii),r_in_update(:,ii),1,1e-2,Y,X);
        end
        Psi_in = Psi_in/sqrt(I0); % Normalize
        
        % 4.1.3. Update the output-updated matrix
        phi_in = angle(Psi_in'*sum(Psi_in,2));
        r_in_update = r_in_update.*exp(1i*phi_in.');
        r_out_update = exp(1i*phi_in).*r_out_update;

        phi_in_tot = phi_in+phi_in_tot;
        
        % 4.2. Optimize output
        % 4.2.1. Build matrix Psi_out containing the image from each output
        Psi_out = single(zeros(Nx*Ny,N));
        for ii = 1:N
            Psi_out(:,ii) = finufft2d3(fy_out(:,ii),fx_out(:,ii),r_out_update(:,ii),1,1e-2,Y,X);
        end
        Psi_out = Psi_out/sqrt(I0); % Normalize
       
        % 4.2.3. Update the input-updated matrix
        phi_out = angle(Psi_out'*sum(Psi_out,2));
        r_in_update = exp(1i*phi_out).*r_in_update;
        r_out_update = r_out_update.*exp(1i*phi_out.');
        
        phi_out_tot = phi_out+phi_out_tot;
                
        % 4.3. Update the maximum phase change
        max_phi = max(vertcat(abs(phi_in),abs(phi_out)),[],'all'); 
    end
    
    % 5. Save the phases    
    phi_in = phi_in_tot; phi_out = phi_out_tot;

    save(""+directory_save+"phi_in_vrm.mat",'phi_in');
    save(""+directory_save+"phi_out_vrm.mat",'phi_out');
end
