% The function that builds the deep-resolved reflection
% matrix for 3D imaging
function build_r_z_3D_class(list_z_class,list_amp,phase_list,directory_r,prefix,list_k0,directory_save,k)
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
    % List of central depth of the subvolumes
    z_mid_list = [];
    for id_sub = 1:length(list_z_class) % Find the middle z of each subvolume
        z_sub = list_z_class{id_sub};
        z_mid = z_sub(round(length(z_sub)/2));
        z_mid_list = vertcat(z_mid_list,z_mid);
    end
    
    % 2. Build matrices
    % This matrix is a concatenation of the r_z at different depths.
    % We call this matrix rz, while the depth-resolved matrix of a
    % single depth is r_z. Don't be confused.
    rz = zeros(N*length(list_z_class),N);
    for i_freq = 1:n_freq
        fprintf("Working with frequency"+i_freq+"\n")
        k0 = list_k0(i_freq); 
        
        % Load the reflection matrix at this frequency
        r = load(""+directory_r+""+prefix+""+i_freq+".mat").r_pad;
        
        % Find kz and Fourier components fz 
        kz_in = sqrt(k0^2-kx_in.^2-ky_in.^2);
        kz_out = -sqrt(k0^2-kx_out.^2-ky_out.^2);
        fz0 = kz_out-kz_in;
        
        % Propagate to the chosen depths using a Kronecker product
        r = repmat(r,length(list_z_class),1).*exp(1i*mod(kron(z_mid_list,fz0),2*pi));
        
        % Compensate dispersion, normalize with laser amplitude then sum r
        % together 
        rz = rz+r*exp(-1i*wrapToPi(phase_list(i_freq)))/list_amp(i_freq);
    end
    
    % Save the r_z matrix of each subvolume and clear
    for id_sub = 1:length(list_z_class)
        r_z = rz((id_sub-1)*N+1:id_sub*N,:);
        save(""+directory_save+"./r_z_3D_class_"+id_sub+".mat",'r_z')
    end
    clear rz r_z
end
