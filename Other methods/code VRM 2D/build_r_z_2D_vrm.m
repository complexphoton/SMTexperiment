% The function that builds the deep-resolved reflection
% matrix for 2D imaging
function build_r_z_2D_vrm(z_target,list_amp,list_k0,directory_save,k)
    % 1. Setting
    % The number of frequency
    n_freq = length(list_k0);
    
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
        k0 = list_k0(i_freq); 
        
        % Load the dispersion corrected image
        r = load(""+directory_save+"r_2D_"+i_freq+"_vrm.mat").r;
        
        % Find kz and Fourier components fz 
        % In the sample
        kz_in = sqrt(k0^2-kx_in.^2-ky_in.^2);
        kz_out = -sqrt(k0^2-kx_out.^2-ky_out.^2);
        fz0 = kz_out-kz_in;
        
        % Propagate to z_target
        r = r.*exp(1i.*mod(fz0*z_target,2*pi))./list_amp(i_freq);
        
        r_z = r_z+r;
    end
    
    % Save and clear
    save(""+directory_save+"./r_z_2D_vrm.mat",'r_z')
end
