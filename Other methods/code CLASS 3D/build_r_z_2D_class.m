% The function that builds the deep-resolved reflection
% matrix for 2D imaging
function build_r_z_2D_class(z_target,list_amp,phase_list,directory_r,prefix,list_k0,directory_save,k)
    % 1. Setting
    % The number of frequency
    n_freq = length(list_k0);
    
    % Input and output k
    kx_out = k(:,1);
    ky_out = k(:,2);
    kx_in = k(:,1).';
    ky_in = k(:,2).';
    N = length(k); % Number of input/output channels
    % The Fourier components
    fx = kx_out-kx_in;
    fy = ky_out-ky_in;
    
    % 2. Build matrices
    r_z = zeros(N,N);
    for i_freq = 1:n_freq
        fprintf("Working with frequency"+i_freq+"\n")
        k0 = list_k0(i_freq);
        
        r = load(""+directory_r+""+prefix+""+i_freq+".mat").r_pad;
        
        % Find kz and Fourier components fz 
        kz_in = sqrt(k0^2-kx_in.^2-ky_in.^2);
        kz_out = -sqrt(k0^2-kx_out.^2-ky_out.^2);
        fz0 = kz_out-kz_in;
        
        % Correct index mismatch then propagate to z_target
        r = r(:).*exp(1i.*(fz0(:)*z_target));
        
        % Reshape it into a square matrix
        r = reshape(r,N,N);
        
        % Compensate dispersion, normalize with laser amplitude then sum r
        % together 
        r_z = r_z+r*exp(-1i*phase_list(i_freq))/list_amp(i_freq);
    end
    
    % Save and clear
    save(""+directory_save+"./r_z_2D_class.mat",'r_z')
    clear r_z
end
