clear 
clc

begin_all = tic;

%% 0. Provide input info

% 1. System coordinate: the system coordinate (x, y, z) will be expressed as
% lx_offset +- lx_im, where lx_im is the length of the image in x and
% lx_offset is the central coordinate is x. The same goes for y and z. Here
% the z is a range of depths surrounding the estimated depth of the image
% z_offset. Unit: micron

lx_offset = single(25); lx_im = single(50);
ly_offset = single(25); ly_im = single(50);
lz_offset = single(1525); lz_im = single(110);

% 2. Provide the path to input data. We recommend saving all required input 
% data at the same directory, including the reflection matrices, k list, k0 
% list, dispersion coefficients of glass and sample, etc
directory_r = ["/media/minh/WD_BLACK/SMT 3D Tb3 1350 and 1150/Tb3 1150/"];

% 3. List of k0 at all frequencies.
list_k0 = single(load(""+directory_r+"list_k0.mat").list_k0); 
n_freq = length(list_k0);
k0_max = max(list_k0,[],'all');

% 4. Dispersion coefficients of glass is coef_n1 and sample is coef_n2
coef_n1 = single(load(""+directory_r+"coef_n1.mat").coef_n1);
coef_n2 = single(load(""+directory_r+"coef_n2_3D.mat").coef_n2);

% 4. Thickness of the coverslip h1 and mirror depth z_mirror (micron)
h1 = 154; z_ref_air = 1150-h1;

% 5. Specify the imaging case, either 2D or 3D. If 2D, after correcting the
% dispersion, we need to find the target's depth, and when reconstructing
% the image, we only have to reconstruct a 2D image
im_case = "3D";

% 6. In our case, ehe reflection matrices are stored in "/directory_r/" as
% ("/prefix"-frequency_id-".mat").r_pad. So we also need to specify 
% the prefix. Also, for the best convenience, in the future, before using
% this code, change the name of the file containing the reflection matrices
% into that format
prefix = "r_pad_";

% 7. Specify the pixel size. % Pixel size dx = dy, dz when running the
% optimization. When reconstructing the image, we can use a different pixel
% size dx_im = dy_im, dz_im (micron)
dx = 0.7; dx_im = 0.2;
dz = 1; dz_im = 0.2;

% 8. The image is divided into smaller zones during the wavefront
% correction procedure. Those zones are overlapped. Specify the
% half-of-overlapping-area here (micron)
overlap2 = 0; 

% 9. If the case is "3D", we need to specify the number of subvolumes in the z
% axis. If "2D" then the number of subvolumes is 1
n_sub = 16;

% 10. List of laser amplitude at all frequencies
list_amp = load(""+directory_r+"/list_amp.mat").list_amp_ave;

% 11. Specify a saving directory to store the intermediate variables which
% can be very heavy. We will store those heavy intermediate variables after
% obtaining them, then clear them.
directory_save = ["/media/minh/WD_BLACK/SMT 3D Tb3 1350 and 1150/Tb3 1150/Intermediate result/"];

% 12. Specify the list of k_parallel_in and k_parallel_out (kin/out_x/y).
% When saving these k, please save them as a nx2 matrix where the first
% column is kx, the 2nd column is ky. For each constant value of kx, the
% values of ky increases from negative to positive. Arranging like that
% helps us easily find the spacing in k space. The reflection matrices are
% also arranged in that order.
k = single(load(""+directory_r+"k.mat").k);

% 13. NA of the system
NA = 0.5;

% 14. % The outer part of the image is always darker than the central part
% due to vignetting. To compensate, we multiply the image with a
% Gaussian mask with the following parameters
x0 = 30; y0 = 22; % Central of the mask
wx0 = 40; wy0 = 22; % Width
rho = 0;
% The mask will be:
% [exp(-1/(2*(1-rho^2))*(((gridX-x0)/wx0).^2+((gridY-y0)/wy0).^2)-...
% 2*rho.*(gridX-x0)/wx0.*(gridY-y0)/wy0)'.^2

% 15. If the target is 3D, how many slices (around the center slice) do you
% want to correct for each subvolume?
n_slice = 5;
if mod(n_slice,2) == 0 
% If the number of slices is an even number, it's a bit inconvenient. It's 
% best to be an odd number so that the numbers of slices above and below the
% center depth of the subvolume are the same
    n_slice = n_slice+1;
end

% 16. Truncation range in z for 2D imaging case
if im_case == "2D"
    range_z = 9;
end
% That's all one need to put in in order to reconstruct the fully corrected
% image. 


%% I. System coordinates, zoning, converting z axis

fprintf("System coordinates, zoning, converting z axis.\n")

% 1. Coordinate
% 1.1. Below is the lists of coordinate x, y, z when running the computation
% and x_im, y_im, z_im when reconstructing the image
x = single(lx_offset-lx_im/2+dx/2:dx:lx_offset+lx_im/2); Nx = length(x);
y = single(ly_offset-ly_im/2+dx/2:dx:ly_offset+ly_im/2); Ny = length(y);
z = single(lz_offset-lz_im/2+dz/2:dz:lz_offset+lz_im/2); Nz = length(z);
x_im = single(lx_offset-lx_im/2+dx_im/2:dx_im:lx_offset+lx_im/2); Nx_im = length(x_im);
y_im = single(ly_offset-ly_im/2+dx_im/2:dx_im:ly_offset+ly_im/2); Ny_im = length(y_im);
z_im = single(lz_offset-lz_im/2+dz_im/2:dz_im:lz_offset+lz_im/2); Nz_im = length(z_im);

% 1.2. Zoning in z for 3D case
[list_z, list_z_im] = zoning_z(z,z_im,dz,dz_im,n_sub,overlap2); % list_z is a cell containing the z of each

% 1.3. If we are doing 2D imaging, scan the z axis to find the depth with
% the maximum intensity to set as the reference plane of VRM
kx_out = k(:,1); ky_out = k(:,2); kx_in = kx_out.'; ky_in = ky_out.';
fx = kx_out-kx_in; fy = ky_out-ky_in;
[X,Y,Z] = meshgrid(x,y,z);
freqc = 3e2*list_k0(round(n_freq/2))/2/pi;
if im_case == "2D"
    psi = single(zeros(Nx,Ny,Nz));
    for i_freq = 1:n_freq
        r = single(load(""+directory_r+""+prefix+""+i_freq+".mat").r_pad);
        
        fprintf("Building uncorrected image at frequency"+i_freq+"\n")
        k0 = list_k0(i_freq); freq = 3e2*k0/2/pi; dfreq = freq - freqc; 
        
        % The refractive indices of glass and the medium at the
        % corresponding frequency
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
        
        % Correct index mismatch. We use mod function here to make the
        % phase term inside the exp function small, otherwise, it will be
        % more time consuming.
        r = r.*exp(1i.*mod(-fz0*(z_ref_air+h1)+fz1*h1,2*pi));
        
        % Build the image at this frequency and scale it with the laser
        % amplitude
        psi_s = finufft3d3(fy(:),fx(:),fz(:),r(:),1,1e-2,Y(:),X(:),Z(:))/list_amp(i_freq);
        
        psi = psi+reshape(psi_s,Nx,Ny,Nz);
    end
    
    I = abs(psi).^2;
    clear psi
    
    list_M_z = [];
    for ii = 1:Nz
        I_2D = I(:,:,ii);
        list_M_z = [list_M_z M];
        
        figure(1)
        imagesc(interp2(I_2D,3))
        axis image
        colormap('hot')
        title(""+z(ii)+"")
        pause(0.1)
    end
   
    clear I
    
    z_max = z(list_M_z == max(list_M_z,[],'all'));

end

%% Do VRM for each subvolume
for id_sub = 1:n_sub
    
    z_sub = list_z{id_sub+1,1}; Nz_sub = length(z_sub);
    z_sub_min = min(z_sub,[],'all');
    z_sub_max = max(z_sub,[],'all');
    z_center = z_sub(round(Nz_sub/2));
    
    % II. Do z truncation
    if im_case == "2D"
        z_truncate_min = z_max-range_z/2;
        z_truncate_max = z_max+range_z/2;
    elseif im_case == "3D"
        % Due to dispersion, the particles will be shifted and/or spread to
        % a larger volume than the subvolume we are imaging
        z_truncate_min = z_sub_min-5;
        z_truncate_max = z_sub_max+5;
    end
    
    fprintf("Do time domain truncation for subvolume "+id_sub+".\n")
    t_truncation(directory_r,prefix,list_k0,z_truncate_min,z_truncate_max,z_center,k,coef_n2,coef_n1,z_ref_air,h1,id_sub);
        
    % After this stage, the reflection matrices already has index
    % mismatch correction. So, in the next step, no need to do that
    % again. 
    
    % III. Compensate for dispersion and build time-gated reflection matrices
    [X,Y,Z_sub] = meshgrid(x,y,z_sub);
    freqc = 3e2*list_k0(round(n_freq/2))/2/pi;

    % Build 3D image
    psi_sub = single(zeros(Nx,Ny,Nz_sub));
    for i_freq = 1:n_freq
        fprintf("Build pre-vrm image at frequency "+i_freq+".\n")
        
        k0 = list_k0(i_freq); freq = 3e2*k0/2/pi; dfreq = freq - freqc; 
        n_media = coef_n2(1) + coef_n2(2)*dfreq + coef_n2(3)*dfreq.^2 + coef_n2(4)*dfreq.^3 + coef_n2(5)*dfreq.^4;
        
        kz_in = (sqrt(k0^2*n_media^2-kx_in.^2-ky_in.^2));
        kz_out = (-sqrt(k0^2*n_media^2-kx_out.^2-ky_out.^2));
        fz = (kz_out-kz_in);
        
        r = load(""+directory_r+"r_truncated_"+i_freq+"_subvolume_"+id_sub+".mat").r;
        
        psi_sub = psi_sub+reshape(finufft3d3(fx(:),fy(:),fz(:),r(:),1,1e-2,X(:),Y(:),Z_sub(:)),Nx,Ny,Nz_sub);
        
    end
    
    I_sub = abs(psi_sub).^2;
    clear psi_sub
    
    list_M = [];
    for ii = 1:Nz_sub
        I_2D = I_sub(:,:,ii);
        M = sum(I_2D,'all');
        list_M = [list_M M];
    end
    
    % The depth with the maximum intensity will be comes the new reference
    % depth when doing VRM for this subvolume
    z_sub_ref = z_sub(list_M == max(list_M,[],'all'));
    
    I_max = I_sub(:,:,z_sub == z_sub_ref);
    
    figure(2)
    imagesc(fliplr(interp2(I_max,3)))
    axis image
    colormap('hot')
    
    clear I_sub
    
    % Do VRM dispersion
    if im_case == "2D"
        % Correct dispersion with VRM
        fprintf("Correcting dispersion. \n")
        dispersion_VRM(directory_r,z_sub_ref,list_k0,k,id_sub,coef_n2);
        % After this step, the focal plane is at z_sub_ref and the
        % depth with the highest signal is also z_sub_ref as well.
        % Build time-gated r matrix
        fprintf("Building time-gated matrices. \n")
        build_r_z_2D_vrm(list_amp,list_k0,directory_save,directory_r,k,id_sub);
    elseif im_case == "3D"
        fprintf("Correcting dispersion. \n")
        dispersion_VRM(directory_r,z_sub_ref,list_k0,k,id_sub,coef_n2);
        fprintf("Building time-gated matrices. \n")    
        build_r_z_3D_vrm(z_sub,z_sub_ref,list_amp,list_k0,directory_save,directory_r,n_slice,k,id_sub,coef_n2);
    end

    % III. Wavefront correction
    if im_case == "2D"
        r_z = single(load(""+directory_save+"r_z_2D_vrm_"+id_sub+".mat").r_z);
        wavefront_correction_vrm_2D(x,y,k,r0,id_sub,directory_save);
    elseif im_case == "3D"
        r_z = load(""+directory_save+"r_z_3D_vrm_"+id_sub+".mat").rz;
        wavefront_correction_vrm_3D(x,y,k,r_z,n_slice,id_sub,directory_save);
    end
    
end

%% IV. Reconstruct image
% Build mask to compensate vignetting
[mask] = build_mask(x_im,y_im,wx0,wy0,x0,y0,rho);

fprintf("Reconstructing the image. \n")
if im_case == "2D"
    [I] = reconstruct2D_vrm(x_im,y_im,k,r_z,directory_save);
    I_show = I.*mask/max(I.*mask,[],'all');
    I_show = fliplr(I_show);
    % Display image
    figure
    imagesc(I_show);
    colormap('hot')
    axis image
    caxis([0 0.8])
    set(gca,'Visible','off')
elseif im_case == "3D"
    [I] = reconstruct3D_vrm(list_z_im,x_im,y_im,dz_im,z_im,overlap2,k,list_k0,list_amp,directory_r,prefix,directory_save,coef_n2,coef_n1,z_ref_air,h1);
    
    % Choose the z (depth, in pixel coordinate) to display 
    % Show enface image: choose the z (depth, in pixel coordinate) to display 
    n_z_show = n_sub;
    z_display = round(Nz_im/(n_z_show+1)*[1:1:n_z_show]);
    for ii = z_display
        I_z = I(:,:,ii);
        I_z = I_z/max(I,[],'all');
        figure
        imagesc(fliplr(I_z));
        colormap('hot')
        axis image
        caxis([0 0.1])
        set(gca,'Visible','off')
    end
    
    % Show cross section image, choose the constant x to display
    n_x_show = 4;
    x_display = round(length(x_im)/(n_x_show+1)*[1:1:n_x_show]);
    for ii = x_display
        I_x = I(ii,:,:);
        I_x = reshape(I_x(:),Ny_im,Nz_im); 
        I_x = I_x/max(I,[],'all');
        figure
        imagesc((I_x));
        colormap('hot')
        axis image
        caxis([0 0.1])
        set(gca,'Visible','off')
    end
end

