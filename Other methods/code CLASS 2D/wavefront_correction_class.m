% Function to correct wavefront
function wavefront_correction_class(x_zone,y_zone,k,r_z,n_division,zone_id,subvolume_id,directory_save)
        
    fprintf("Correcting zone "+zone_id+" subvolume "+subvolume_id+".\n");
    % Inputs: 
    % x,y spatial coordinates of the zone
    % r_z: updated time-gated reflection matrix of the bigger zone containing this zone
    Nx = length(x_zone); Ny = length(y_zone);
    % Build grid points 
    [X,Y] = meshgrid(single(x_zone),single(y_zone));
    X = X(:); Y = Y(:);
    
    % 1. Do basis conversion
    % Input and output k are:
    N = length(k); 
    kout_x = k(:,1); kout_y = k(:,2);
    kin_x = k(:,1).'; kin_y = k(:,2).';
    % The Fourier components
    fx = kout_x-kin_x;
    fy = kout_y-kin_y;
    
    dx = abs(x_zone(1)-x_zone(2));
    dk = abs(k(1,2)-k(2,2));

    % 1.2. Do basis conversion
    if n_division > 0        
        % 1.2.1. Convert to truncated spatial basis
        % Convert the input to spatial basis
        r_kout_r = single(zeros(N,Nx*Ny)); % The reflection matrix whose input is spatial basis, output is angular basis
        parfor ii = 1:N
            r_kout_r(ii,:) = finufft2d3(kin_y.',kin_x.',r_z(ii,:).',-1,1e-2,Y,X).';
        end
        r_kout_r = 1/(2*pi)*r_kout_r*dk^2;
        % Convert the output to spatial basis
        r_r = single(zeros(Nx*Ny,Nx*Ny)); % The matrix whose both input and output are in spatial basis
        parfor ii = 1:Nx*Ny
            r_r(:,ii) = finufft2d3(-kout_y,-kout_x,r_kout_r(:,ii),-1,1e-2,Y,X);
        end
        r_r = 1/(2*pi)*r_r*dk^2;
        % 1.2.2. Convert back to angular basis
        % Convert the input to angular basis
        r_r_kin = single(zeros(Nx*Ny,N)); % The matrix whose input is in angular basis, output is spatial basis
        parfor ii = 1:Nx*Ny
            r_r_kin(ii,:) = finufft2d3(Y,X,r_r(ii,:).',1,1e-2,kin_y.',kin_x.').';
        end
        r_r_kin = 1/(2*pi)*r_r_kin*dx^2;
        % Convert the output to angular basis
        r_k = single(zeros(N,N)); % The matrix in angular basis
        parfor ii = 1:N
            r_k(:,ii) = finufft2d3(Y,X,r_r_kin(:,ii),1,1e-2,-kout_y,-kout_x);
        end
        r_k = 1/(2*pi)*r_k*dx^2;
    elseif n_division == 0
        r_k = r_z;
        clear r_z
    end
    
    % 2. The initial image
    psi_in = finufft2d3(fy(:),fx(:),r_k(:),1,1e-2,Y,X);
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
    r_in_update = r_k;
    r_out_update = r_k.';
    % To be more convenient when working with the input and output-updated
    % matrices, we define
    fx_in = fx; fy_in = fy;
    fx_out = fx.'; fy_out = fy.';
    
    % 4. Optimization
    delta_M = 1; % The change of FOM after each optimization loop
    while delta_M > 0.05
        M_prev = M; % Update the FOM of the previous loop
        
        % 4.1. Optimize input
        % 4.1.1. Build matrix Psi_in containing the image from each input
        Psi_in = single(zeros(Nx*Ny,N));
        parfor ii = 1:N
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
        parfor ii = 1:N
            Psi_out(:,ii) = finufft2d3(fy_out(:,ii),fx_out(:,ii),r_out_update(:,ii),1,1e-2,Y,X);
        end
        Psi_out = Psi_out/sqrt(I0); % Normalize
       
        % 4.2.3. Update the input-updated matrix
        phi_out = angle(Psi_out'*sum(Psi_out,2));
        r_in_update = exp(1i*phi_out).*r_in_update;
        r_out_update = r_out_update.*exp(1i*phi_out.');
        
        phi_out_tot = phi_out+phi_out_tot;
        
        psi = Psi_out*exp(1i*phi_out);
        I = abs(psi).^2;
        M = sum(I,'all');
                
        % 4.3. Update the FOM change
        delta_M = (M-M_prev)/abs(M_prev); 
    end
    
    % 5. Update the reflection matrix of the zone, save the matrix and the phases
    r_update = r_in_update;   
    
    phi_in = phi_in_tot; phi_out = phi_out_tot;
   
    psi = finufft2d3(fy(:),fx(:),r_update(:),1,1e-2,Y,X);
    I = abs(psi).^2;
    I_zone = reshape(I,Ny,Nx);

    figure(2)
    imagesc(fliplr(interp2(I_zone,2.5)));
    axis image
    colormap('hot')
   
    save(""+directory_save+"phi_in_"+n_division+"_zone_"+zone_id+"_subvolume_"+subvolume_id+"_class.mat",'phi_in');
    save(""+directory_save+"phi_out_"+n_division+"_zone_"+zone_id+"_subvolume_"+subvolume_id+"_class.mat",'phi_out');
end
