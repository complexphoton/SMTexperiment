function t_truncation(directory_r,prefix,list_k0,z_truncate_min,z_truncate_max,z_center,k,coef_n2,coef_n1,z_ref_air,h1,id_sub)

    % 1. Info
    % Central frequency
    n_freq = length(list_k0);    
    freqc = 3e2*list_k0(round(n_freq/2))/2/pi;

    % k list
    kx_in = k(:,1).'; ky_in = k(:,2).'; kx_out = k(:,1); ky_out = k(:,2);
    N = length(k);

    % 2. Pick the range of z that doesn't contain the objective lens and build the time-gated reflection matrix
    k0_min = list_k0(end); k0_max = list_k0(1);
    Delta_k0 = k0_max-k0_min;
    Delta_omega = Delta_k0*300;
    dt_truncation = 2*pi/Delta_omega; 
    dz_truncation = dt_truncation*(300/coef_n2(1))/2; 
    z_truncate = single(z_truncate_min:dz_truncation:z_truncate_max); Nz_truncate = length(z_truncate);
    
    % 3. Do conversion from frequency to z
    R_z = cell(Nz_truncate,1);
    for ii = 1:Nz_truncate
        R_z{ii,1} = single(zeros(N,N));
    end
    for i_freq = 1:n_freq            
        k0 = list_k0(i_freq); omega = k0*300; freq = 3e2*k0/2/pi; dfreq = freq - freqc; 
        r = single(load(""+directory_r+""+prefix+""+i_freq+".mat").r_pad);
        
        % Refractive index
        n1 = coef_n1(1) + coef_n1(2)*dfreq + coef_n1(3)*dfreq.^2 + coef_n1(4)*dfreq.^3 + coef_n1(5)*dfreq.^4;
        n_media = coef_n2(1) + coef_n2(2)*dfreq + coef_n2(3)*dfreq.^2 + coef_n2(4)*dfreq.^3 + coef_n2(5)*dfreq.^4;
        
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
        
        for ii = 1:Nz_truncate
            z_ii = z_truncate(ii);
            time = (z_ii-z_center)/(300/coef_n2(1));
            fprintf("Building rz at frequency "+i_freq+", z = "+z_ii+".\n")
            R_z{ii,1} = R_z{ii,1}+1/n_freq*exp(-2*1i*omega*time)*r.*exp(1i*mod(fz*z_center,2*pi));
        end
    end
    
    % 4. Rebuild the wavelength reflection matrix
    for i_freq = 1:n_freq
        
        k0 = list_k0(i_freq); omega = k0*300; freq = 3e2*k0/2/pi; dfreq = freq - freqc; 
        
        % Refractive index
        n_media = coef_n2(1) + coef_n2(2)*dfreq + coef_n2(3)*dfreq.^2 + coef_n2(4)*dfreq.^3 + coef_n2(5)*dfreq.^4;
               
        % Fourier components
        kz_in = (sqrt(k0^2*n_media^2-kx_in.^2-ky_in.^2));
        kz_out = (-sqrt(k0^2*n_media^2-kx_out.^2-ky_out.^2));
        fz = (kz_out-kz_in);

        r = single(zeros(N,N));
        for ii = 1:Nz_truncate
            z_ii = z_truncate(ii);
            time = (z_ii-z_center)/(300/coef_n2(1));
            fprintf("Building R at frequency "+i_freq+", z = "+z_ii+", subvolume "+id_sub+".\n")
            
            r_z = R_z{ii};

            r = r+exp(2*1i*omega*time)*r_z.*exp(-1i*mod(fz*z_center,2*pi));
        end

        save(""+directory_r+"r_truncated_"+i_freq+"_subvolume_"+id_sub+".mat",'r')
        
    end
    
    clear R_z
   
end