function [z_target] = find_z_target(x,y,z,list_k0,k,im_case,directory_save)
    % 1. Setting
    % The number of frequency
    n_freq = length(list_k0);

    % Input and output k
    kx_out = k(:,1);
    ky_out = k(:,2);
    kx_in = k(:,1).';
    ky_in = k(:,2).';
    % The Fourier components
    fx = single(kx_out-kx_in);
    fy = single(ky_out-ky_in);
    % The coordinate is built into grid points
    [X,Y,Z] = meshgrid(x,y,z);
    Nx = length(x);
    Ny = length(y);
    Nz = length(z);
    % 2. Build 3D image
    psi = zeros(Nx,Ny,Nz);
    for i_freq = 1:n_freq        
        fprintf("Working with frequency"+i_freq+"\n")
        k0 = list_k0(i_freq); 
        
        % Load the reflection matrix at this frequency
        r = load(""+directory_save+"r_"+im_case+"_"+i_freq+"_vrm.mat").r;
        
        % Find kz and Fourier components fz 
        kz_in = (sqrt(k0^2-kx_in.^2-ky_in.^2));
        kz_out = (-sqrt(k0^2-kx_out.^2-ky_out.^2));
        fz0 = (kz_out-kz_in);
        
        % Build the image at this frequency and scale it with the laser
        % amplitude
        fft_time = tic;
        psi_s = finufft3d3(fy(:),fx(:),fz0(:),r(:),1,1e-2,Y(:),X(:),Z(:));
        toc(fft_time)
        
        psi = psi+reshape(psi_s,Nx,Ny,Nz);
        
    end
    I = abs(psi).^2;
    
    % 3. Search for target depth
    list_M = []; % list of FOM of 2D image at each depth
    for ii = 1:Nz % Scan each depth and find the FOM of the 2D images
        I_2D = I(:,:,ii);
        M = sum(I_2D,'all');
        list_M = [list_M M];    
    end
    % Search for target depth
    z_target = z(list_M == max(list_M,[],'all'));
end