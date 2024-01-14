% The function that builds the deep-resolved reflection
% matrix for 2D imaging
function build_r_z_2D(z_target,list_amp,phase_list,directory_r,prefix,h1,z_mirror,coef_n1,coef_n2,list_k0,directory_save,k)
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
    N = length(k); % Number of input/output channels
    
    % 2. Build matrices
    r_z = zeros(N,N);
    for i_freq = 1:n_freq
        fprintf("Working with frequency"+i_freq+"\n")
        k0 = list_k0(i_freq); freq = 3e2*k0/2/pi; dfreq = freq - freqc; 
        
        % The refractive indices of glass and the medium at the
        % corresponding frequency
        n1 = coef_n1(1) + coef_n1(2)*dfreq + coef_n1(3)*dfreq.^2 + coef_n1(4)*dfreq.^3 + coef_n1(5)*dfreq.^4;
        n_media = coef_n2(1) + coef_n2(2)*dfreq + coef_n2(3)*dfreq.^2 + coef_n2(4)*dfreq.^3 + coef_n2(5)*dfreq.^4;
        
        % Load the reflection matrix at this frequency
        load_data = tic;
        r = load(""+directory_r+""+prefix+""+i_freq+"_single_precision.mat").r_pad;
        toc(load_data)
        
        compute = tic;
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
        
        % Correct index mismatch then propagate to z_target
        r = r(:).*exp(1i.*mod(-fz0(:)*z_mirror+fz1(:)*h1+fz(:)*z_target,2*pi));
        
        % Reshape it into a square matrix
        r = reshape(r,N,N);
        
        % Compensate dispersion, normalize with laser amplitude then sum r
        % together 
        r_z = r_z+r*exp(-1i*phase_list(i_freq))/list_amp(i_freq);
        toc(compute)
    end
    
    % Save and clear
    save(""+directory_save+"./r_z_2D.mat",'r_z')
    clear r_z
end
