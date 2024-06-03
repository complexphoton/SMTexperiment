% Function to correct wavefront
function wavefront_correction_vrm_3D(x,y,k,r_z,n_slice,subvolume_id,directory_save)
        
    fprintf("Correcting subvolume "+subvolume_id+".\n");
    % Inputs: 
    % x,y spatial coordinates of the zone
    % r_z: updated time-gated reflection matrix of the bigger zone containing this zone
    Nx = length(x); Ny = length(y);
    % Build grid points 
    [X,Y] = meshgrid(single(x),single(y));
    X = X(:); Y = Y(:);
    % Number of angles
    N = length(k); 
    
    % 1. Do basis conversion
    % Input and output k are:
    kout_x = k(:,1); kout_y = k(:,2);
    kin_x = k(:,1).'; kin_y = k(:,2).';
    % The Fourier components
    fx = kout_x-kin_x;
    fy = kout_y-kin_y;
    
    % 2. The initial image
    I_in = [];
    for slice_id = 1:n_slice
        r_slice = r_z((slice_id-1)*N+1:slice_id*N,:);
        psi_in_slice = finufft2d3(fy(:),fx(:),r_slice(:),1,1e-2,Y,X);
        I_in_slice = abs(psi_in_slice).^2;
        I_in = vertcat(I_in,I_in_slice);
    end
    
    % Normalize to avoid the FOM getting too large
    I0 = 1; % Normalization constant
    I_in = I_in/I0; 
    
    % FOM of the initial image
    global M
    M = sum(I_in,'all');
    
    I_in = reshape(I_in,Ny,Nx,n_slice);
    I_in_2D = I_in(:,:,ceil(n_slice/2));
    
    figure(1)
    imagesc(fliplr(interp2(I_in_2D,3)))
    axis image
    colormap('hot')
    
    % 3. Initialize
    phi_in_tot = zeros(N,1); phi_out_tot = zeros(N,1); 
    % The input and output-updated matrices
    r_in_update = r_z;
    r_out_update = [];
    for slice_id = 1:n_slice
        r_slice = r_z((slice_id-1)*N+1:slice_id*N,:);
        r_out_update = vertcat(r_out_update,r_slice.');
    end
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
        Psi_in = [];
        for slice_id = 1:n_slice
            r_in_update_slice = r_in_update((slice_id-1)*N+1:slice_id*N,:);
            parfor ii = 1:N
                Psi_in_slice(:,ii) = finufft2d3(fy_in(:,ii),fx_in(:,ii),r_in_update_slice(:,ii),1,1e-2,Y,X);
            end
            Psi_in = vertcat(Psi_in,Psi_in_slice);
        end
        Psi_in = Psi_in/sqrt(I0); % Normalize
        
        % 4.1.3. Update the output-updated matrix
        phi_in = angle(Psi_in'*sum(Psi_in,2));
        r_in_update = r_in_update.*exp(1i*phi_in.');
        r_out_update = repmat(exp(1i*phi_in),n_slice,1).*r_out_update;

        phi_in_tot = phi_in+phi_in_tot;
        
        % 4.2. Optimize output
        % 4.2.1. Build matrix Psi_out containing the image from each output
        Psi_out = [];
        for slice_id = 1:n_slice
            r_out_update_slice = r_out_update((slice_id-1)*N+1:slice_id*N,:);
            parfor ii = 1:N
                Psi_out_slice(:,ii) = finufft2d3(fy_out(:,ii),fx_out(:,ii),r_out_update_slice(:,ii),1,1e-2,Y,X);
            end
            Psi_out = vertcat(Psi_out,Psi_out_slice);
        end
        Psi_out = Psi_out/sqrt(I0); % Normalize
       
        % 4.2.3. Update the input-updated matrix
        phi_out = angle(Psi_out'*sum(Psi_out,2));
        r_in_update = repmat(exp(1i*phi_out),n_slice,1).*r_in_update;
        r_out_update = r_out_update.*exp(1i*phi_out.');
        
        phi_out_tot = phi_out+phi_out_tot;
        
        psi = Psi_out*exp(1i*phi_out);
        I = abs(psi).^2;
                
        % 4.3. Update the maximum phase change
        max_phi = max(vertcat(abs(phi_in),abs(phi_out)),[],'all'); 
    end
    
    % 5. Update the reflection matrix of the zone, save the matrix and the
    % phases
    phi_in = phi_in_tot; phi_out = phi_out_tot;

    I = reshape(I,Ny,Nx,n_slice);
    I_2D = I(:,:,ceil(n_slice/2));

    figure
    imagesc(fliplr(interp2(I_2D,2.5)));
    axis image
    colormap('hot')
   
    save(""+directory_save+"phi_in_subvolume_"+subvolume_id+"_vrm.mat",'phi_in');
    save(""+directory_save+"phi_out_subvolume_"+subvolume_id+"_vrm.mat",'phi_out');
end
