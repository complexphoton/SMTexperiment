% The function that builds the deep-resolved reflection
% matrix for 3D imaging
function build_r_z_3D_vrm(z_sub,z_sub_ref,list_amp,list_k0,directory_save,directory_r,n_slice,k,id_sub,coef_n2)
    % 1. Setting
    % The number of frequency
    n_freq = length(list_k0);
    % Center frequency
    freqc = 3e2*list_k0(round(n_freq/2))/2/pi;
    % Input and output k
    kx_out = k(:,1);
    ky_out = k(:,2);
    kx_in = k(:,1).';
    ky_in = k(:,2).';
    N = length(k); % Number of input/output channels
    
    % 2. Build matrices for each subvolume
    % List of z that is corrected in this subvolume
    list_z = linspace(z_sub(1),z_sub(end),n_slice).';
    % This matrix is a concatenation of the r_z at different depths.
    % We call this matrix rz, while the depth-resolved matrix of a
    % single depth is r_z. Don't be confused.
    rz = zeros(N*length(list_z),N);
    for i_freq = 1:n_freq
        fprintf("Working with frequency"+i_freq+", subvolume "+id_sub+"\n")
        k0 = list_k0(i_freq); freq = 3e2*k0/2/pi; dfreq = freq - freqc; 
        
        % The refractive indices of the medium at the
        % corresponding frequency
        n_media = coef_n2(1) + coef_n2(2)*dfreq + coef_n2(3)*dfreq.^2 + coef_n2(4)*dfreq.^3 + coef_n2(5)*dfreq.^4;

        % Load the reflection matrix at this frequency. The reflection
        % matrix has the focal plane at z_sub_ref
        r = load(""+directory_r+"r_vrm_"+i_freq+"_subvolume_"+id_sub+".mat").r;

        kz_in = (sqrt(k0^2*n_media^2-kx_in.^2-ky_in.^2));
        kz_out = (-sqrt(k0^2*n_media^2-kx_out.^2-ky_out.^2));
        fz = (kz_out-kz_in);

        % Propagate to the chosen depths using a Kronecker product
        r = repmat(r,length(list_z),1).*exp(1i*mod(kron(list_z-z_sub_ref,fz),2*pi));

        % Compensate dispersion, normalize with laser amplitude then sum r
        % together 
        rz = rz+r./list_amp(i_freq);
    end

    % Save the rz matrix of each subvolume and clear
    save(""+directory_save+"./r_z_3D_vrm_"+id_sub+".mat",'rz')
    clear rz
    
end
