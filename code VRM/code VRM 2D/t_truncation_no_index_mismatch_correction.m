function t_truncation_no_index_mismatch_correction(directory_r,prefix,list_k0,z_truncate_min,z_truncate_max,z_ref_air,z_prime_max,k,coef_n2,coef_n1,h1)

    % 1. Info
    % Central frequency
    n_freq = length(list_k0);    

    % k list
    kx_in = k(:,1).'; ky_in = k(:,2).'; kx_out = k(:,1); ky_out = k(:,2);
    N = length(k);
    n2 = coef_n2(1);
    n1 = coef_n1(1);
    
    % 2. Pick the range of z that doesn't contain the objective lens and build the time-gated reflection matrix
    k0_min = list_k0(end); k0_max = list_k0(1);
    Delta_k0 = k0_max-k0_min;
    Delta_omega = Delta_k0*300;
    dt_truncation = 2*pi/Delta_omega; 
    dz_truncation = dt_truncation*(300/n2)/2; 
    z_truncate = single(z_truncate_min:dz_truncation:z_truncate_max); Nz_truncate = length(z_truncate); % z truncate is z in medium
    
    % 3. Do conversion from frequency to z
    R_z = cell(Nz_truncate,1);
    for ii = 1:Nz_truncate
        R_z{ii,1} = single(zeros(N,N));
    end
    for i_freq = 1:n_freq            
        k0 = list_k0(i_freq); omega = k0*300;
        
        r = single(load(""+directory_r+""+prefix+""+i_freq+".mat").r_pad);
        
        % Find kz and Fourier components fz 
        % In air
        kz_in0 = (sqrt(k0^2-kx_in.^2-ky_in.^2));
        kz_out0 = (-sqrt(k0^2-kx_out.^2-ky_out.^2));
        fz0 = (kz_out0-kz_in0);
       
        for ii = 1:Nz_truncate
            z_ii = z_truncate(ii);
            time = (n2*z_ii-(z_prime_max+z_ref_air+h1)+n1*h1)/300;
            fprintf("Building rz at frequency "+i_freq+", z = "+z_ii+".\n")
            R_z{ii,1} = R_z{ii,1}+1/n_freq*exp(-2*1i*omega*time)*r.*exp(1i*mod(fz0*z_prime_max,2*pi));
        end
    end
    
    % 4. Rebuild the wavelength reflection matrix
    for i_freq = 1:n_freq
        
        k0 = list_k0(i_freq); omega = k0*300; 
        
        % Fourier components
        kz_in0 = (sqrt(k0^2-kx_in.^2-ky_in.^2));
        kz_out0 = (-sqrt(k0^2-kx_out.^2-ky_out.^2));
        fz0 = (kz_out0-kz_in0);

        r = single(zeros(N,N));
        for ii = 1:Nz_truncate
            z_ii = z_truncate(ii); 
            time = (n2*z_ii-(z_prime_max+z_ref_air+h1)+n1*h1)/300;
            fprintf("Building R at frequency "+i_freq+", z = "+z_ii+".\n")
            
            r_z = R_z{ii};

            r = r+exp(2*1i*omega*time)*r_z.*exp(-1i*mod(fz0*z_prime_max,2*pi));
        end

        save(""+directory_r+"r_truncated_"+i_freq+".mat",'r')
        
    end
    
    clear R_z
   
end