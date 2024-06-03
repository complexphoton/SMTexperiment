function dispersion_VRM_no_index_mismatch_correction(directory_r,directory_save,z_max,z_prime_max,z_ref_air,list_k0,k,coef_n1,coef_n2,h1)
    
    % 1. Setting
    % k space
    kx = k(:,1); ky = k(:,2); N = length(k);
    % The number of frequency
    n_freq = length(list_k0);
    n1 = coef_n1(1);
    n2 = coef_n2(1);
    
    % The central frequency
    i_freq_center = round(n_freq/2); k0 = list_k0(i_freq_center); 
    omega = 3e2*k0;
    
    % Load reflection matrix of the center frequency
    r_center = load(""+directory_r+"r_truncated_"+i_freq_center+".mat").r;
    
    % Propagate to z_air_max
    kz = (sqrt(k0^2-kx.^2-ky.^2));
    fz0 = (-kz-kz.');
    time = (n2*z_max-(z_prime_max+z_ref_air+h1)+n1*h1)/300;
    r_center = exp(-2*1i*omega*time)*r_center.*exp(1i*fz0*z_prime_max);
    
    % 2. Build matrices
    for i_freq = 1:n_freq
        fprintf("Working with frequency"+i_freq+"\n")
        
        k0 = list_k0(i_freq); omega = 300*k0;
        
        % Load the reflection matrix at this frequency
        r = load(""+directory_r+"r_truncated_"+i_freq+".mat").r;
        % Propagate to z_sub_ref
        kz = (sqrt(k0^2-kx.^2-ky.^2));
        fz0 = (-kz-kz.');
        time = (n2*z_max-(z_prime_max+z_ref_air+h1)+n1*h1)/300;
        r = exp(-2*1i*omega*time)*r.*exp(1i*fz0*z_prime_max);
        
        % Do the iterative dispersion correction
        max_phi = pi; N = width(r);
        phi_out = zeros(N,1); phi_in = zeros(N,1);
        while max_phi > pi/18
            
            % Correct output dispersion
            phi_out_list = angle(diag(r*r_center'));
            
            r = diag(exp(-1i*phi_out_list))*r;
            
            % Correct input dispersion
            phi_in_list = angle(diag(transpose(r)*conj(r_center)));
            
            r = r*diag(exp(-1i*phi_in_list));
            
            phi_list = vertcat(phi_in_list(phi_in_list ~= 0), phi_out_list(phi_out_list ~= 0));
            max_phi = max(abs(phi_list),[],'all')
            
            phi_out = phi_out+phi_out_list; phi_in = phi_in+phi_in_list;
        end
        
        % Save the corrected matrix
        save(""+directory_save+"r_vrm_"+i_freq+".mat",'r','phi_in','phi_out')
    end
    
end
