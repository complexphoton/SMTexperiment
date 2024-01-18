% The function that builds the deep-resolved reflection
% matrix for 3D imaging
function build_r_z_3D(list_z,list_amp,phase_list,directory_r,prefix,h1,z_mirror,coef_n1,coef_n2,list_k0,n_slice,n_sub,directory_save,k)
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
    R_z = cell(n_sub,1);
    for id_sub = 1:n_sub
        % List of z that is corrected in this subvolume
        z_sub = list_z{id_sub}; Nz_sub = length(z_sub); z_center = z_sub(round(Nz_sub/2));
        list_z_subvolume = linspace(z_sub(1),z_sub(end),n_slice).'; 
        r_z = zeros(N*n_slice,N);
        for i_freq = 1:n_freq
            tic
            fprintf("Working with frequency "+i_freq+", subvolume "+id_sub+"\n")
            k0 = list_k0(i_freq); freq = 3e2*k0/2/pi; dfreq = freq - freqc; 

            % The refractive indices of glass and the medium at the
            % corresponding frequency
            n1 = coef_n1(1) + coef_n1(2)*dfreq + coef_n1(3)*dfreq.^2 + coef_n1(4)*dfreq.^3 + coef_n1(5)*dfreq.^4;
            n_media = coef_n2(1) + coef_n2(2)*dfreq + coef_n2(3)*dfreq.^2 + coef_n2(4)*dfreq.^3 + coef_n2(5)*dfreq.^4;
        
            % Load the reflection matrix at this frequency
            r = load(""+directory_r+""+prefix+""+i_freq+".mat").r_pad;

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
            
            % Correct index mismatch and propagate to the chosen depths using a Kronecker product
            if n_slice ~= 1
                r = repmat(r.*exp(1i.*mod(-fz0*z_mirror+fz1*h1,2*pi)),n_slice,1).*exp(1i*mod(kron(list_z_subvolume,fz),2*pi));
            else
                r = repmat(r.*exp(1i.*mod(-fz0*z_mirror+fz1*h1,2*pi)),n_slice,1).*exp(1i*mod(fz*z_center,2*pi));
            end
            % Compensate dispersion, normalize with laser amplitude then sum r
            % together 
            r_z = r_z+r*exp(-1i*phase_list(i_freq))/list_amp(i_freq);
            toc
        end
        R_z{id_sub,1} = r_z;
    end
    for id_sub = 1:n_sub
        r_z = R_z{id_sub,1};
        save(""+directory_save+"./r_z_3D_"+id_sub+"_"+n_slice+"_slices_"+n_sub+"_subvolumes.mat",'r_z')
    end
end
