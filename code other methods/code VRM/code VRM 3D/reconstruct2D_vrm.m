% Function to reconstruct the fully corrected 2D SMT image
function [I] = reconstruct2D_vrm(x_im,y_im,k,r_z,directory_save)

    I = zeros(length(y_im),length(x_im));
    subvolume_id = 1;
    % Input and output k:
    kout_x = k(:,1); kout_y = k(:,2);
    kin_x = k(:,1).'; kin_y = k(:,2).';
    % The Fourier components
    fx = kout_x-kin_x;
    fy = kout_y-kin_y;
    % Space
    Nx = length(x_im); % Number of pixels in x and y of the zone
    Ny = length(y_im);
    [X,Y] = meshgrid(single(x_im),single(y_im)); % Make grid points

    % Find the aberration phase of the zone 
    phi_in = load(""+directory_save+"phi_in_"+subvolume_id+"_vrm.mat").phi_in;
    phi_out = load(""+directory_save+"phi_out_"+subvolume_id+"_vrm.mat").phi_out;
    % Update the reflection matrix
    r_update = exp(1i*phi_out).*r_z.*exp(1i*phi_in.');
    % Build the image
    psi = finufft2d3(fy(:),fx(:),r_update(:),1,1e-2,Y(:),X(:));
    I = abs(psi).^2;
    I = reshape(I,Ny,Nx);
    
    % Save the image
    save(""+directory_save+"./I_2D_vrm.mat",'I')
end

