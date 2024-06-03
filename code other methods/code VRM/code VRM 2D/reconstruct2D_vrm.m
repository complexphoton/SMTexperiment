% Function to reconstruct the fully corrected 2D SMT image
function [I] = reconstruct2D_vrm(x_im,y_im,k,directory_r,z_max,z_prime_max,z_ref_air,h1,coef_n1,coef_n2,list_k0,directory_save,prefix)

    n_freq = length(list_k0);
    N = length(k);
    n1 = coef_n1(1); n2 = coef_n2(1);
    % Input and output k:
    kx_out = k(:,1); ky_out = k(:,2);
    kx_in = k(:,1).'; ky_in = k(:,2).';
    % The Fourier components
    fx = kx_out-kx_in;
    fy = ky_out-ky_in;
    % Space
    Nx = length(x_im); % Number of pixels in x and y of the zone
    Ny = length(y_im);
    [X,Y] = meshgrid(single(x_im),single(y_im)); % Make grid points

    % Build r_z
    r_z = single(zeros(N,N));
    
    for i_freq = 1:n_freq
        fprintf("Build image at frequency "+i_freq+".\n")
        k0 = list_k0(i_freq); omega = 300*k0;
        kz_in = (sqrt(k0^2-kx_in.^2-ky_in.^2));
        kz_out = (-sqrt(k0^2-kx_out.^2-ky_out.^2));
        fz0 = (kz_out-kz_in);
        
        r = single(load(""+directory_r+""+prefix+""+i_freq+".mat").r_pad);

        time = (n2*z_max-(z_prime_max+z_ref_air+h1)+n1*h1)/300;

        % Propagate to z_prime_max
        r = exp(-2*1i*mod(omega*time,2*pi))*r.*exp(1i*mod(fz0*z_prime_max,2*pi));
        
        % Apply VRM dispersion phase
        phi_in_dis = load(""+directory_save+"r_vrm_"+i_freq+".mat").phi_in;
        phi_out_dis = load(""+directory_save+"r_vrm_"+i_freq+".mat").phi_out;
        r = exp(-1i*phi_out_dis).*r.*exp(-1i*phi_in_dis.');
        
        % Sum over frequencies
        r_z = r_z+r;
    end
    
    % Find the aberration phase of the zone 
    phi_in = load(""+directory_save+"phi_in_vrm.mat").phi_in;
    phi_out = load(""+directory_save+"phi_out_vrm.mat").phi_out;
   
    % Update the reflection matrix
    r_vrm = exp(1i*phi_out).*r_z.*exp(1i*phi_in.');
    % Build the image
    % With wavefront correction
    psi = finufft2d3(fy(:),fx(:),r_vrm(:),1,1e-2,Y(:),X(:));
    I = abs(psi).^2;
    I = reshape(I,Ny,Nx);
    % Without wavefront correction but with dispersion compensation
    psi_disp = finufft2d3(fy(:),fx(:),r_z(:),1,1e-2,Y(:),X(:));
    I_disp = abs(psi_disp).^2;
    I_disp = reshape(I_disp,Ny,Nx);
    
    % Save the image
    save(""+directory_save+"./I_2D_vrm.mat",'I','I_disp')
    save(""+directory_save+"./r_vrm.mat",'r_vrm')

end

