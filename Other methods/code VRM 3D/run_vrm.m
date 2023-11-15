clear 
clc

begin_all = tic;

%% 0. Provide input info

% 1. System coordinate: the system coordinat (x, y, z) will be expressed as
% lx_offset +- lx_im, where lx_im is the length of the image in x and
% lx_offset is the central coordinate is x. The same goes for y and z. Here
% the z is a range of depths surrounding the estimated depth of the image
% z_offset. Unit: micron

lx_offset = 25; lx_im = 50;
ly_offset = 25; ly_im = 50;
lz_offset = 1220; lz_im = 110;

% 2. Provide the path to input data. We recommend saving all required input 
% data at the same directory, including the reflection matrices, k list, k0 
% list, dispersion coefficients of glass and sample, etc
directory_r = ["/media/minh/WD_BLACK/smt exp data 3d/"];

% 3. List of k0 at all frequencies.
list_k0 = load(""+directory_r+"list_k0_3D.mat").list_k0;
k0_max = max(list_k0,[],'all');

% 3. Dispersion coefficients of glass is coef_n1 and sample is coef_n2
coef_n1 = load(""+directory_r+"coef_n1.mat").coef_n1;
coef_n2 = load(""+directory_r+"coef_n2_3D.mat").coef_n2;

% 4. Thickness of the coverslip h1 and mirror depth z_mirror (micron)
h1 = 148; z_ref_air = 934-h1;

% 5. In our case, ehe reflection matrices are stored in "/directory_r/" as
% ("/prefix"-frequency_id-".mat").r_pad. So we also need to specify 
% the prefix. Also, for the best convenience, in the future, before using
% this code, change the name of the file containing the reflection matrices
% into that format
prefix = "r_pad_3D_freq_";

% 6. Specify the pixel size. % Pixel size dx = dy, dz when running the
% optimization. When reconstructing the image, we can use a different pixel
% size dx_im = dy_im, dz_im (micron)
dx = 0.5; dx_im = 0.2;
dz = 1; dz_im = 0.2;

% 7. The image is divided into smaller zones during the wavefront
% correction procedure. Those zones are overlapped. Specify the
% half-of-overlapping-area here (micron)
overlap2 = 1.2; 

% 8. Specify the imaging case, either 2D or 3D. If 2D, after correcting the
% dispersion, we need to find the target's depth, and when reconstructing
% the image, we only have to reconstruct a 2D image
im_case = "3D";

% If the case is "3D", we need to specify the number of subvolumes in the z
% axis. If "2D" then the number of subvolumes is 1
n_sub = 7;

% 9. The image will go through one correction of the full image, then
% n_divison division steps
n_division = 0;

% 11. List of laser amplitude at all frequencies
list_amp = load(""+directory_r+"/list_amp_3D.mat").list_amp;

% 12. Specify a saving directory to store the intermediate variables which
% can be very heavy. We will store those heavy intermediate variables after
% obtaining them, then clear them.
directory_save = ["/media/minh/WD_BLACK/smt exp data 3d/Intermediate data/vrm/"];

% 13. Specify the list of k_parallel_in and k_parallel_out (kin/out_x/y).
% When saving these k, please save them as a nx2 matrix where the first
% column is kx, the 2nd column is ky. For each constant value of kx, the
% values of ky increases from negative to positive. Arranging like that
% helps us easily find the spacing in k space. The reflection matrices are
% also arranged in that order.
k = single(load(""+directory_r+"kout_max_3D.mat").kout_max);

% 14. NA of the system
NA = 0.5;

% 15. % The outer part of the image is always darker than the central part
% due to vignetting. To compensate, we multiply the image with a
% Gaussian mask with the following parameters
x0 = 30; y0 = 22; % Central of the mask
wx0 = 40; wy0 = 22; % Width
rho = 0;
% The mask will be:
% [exp(-1/(2*(1-rho^2))*(((gridX-x0)/wx0).^2+((gridY-y0)/wy0).^2)-...
% 2*rho.*(gridX-x0)/wx0.*(gridY-y0)/wy0)'.^2

% 16. If the target is 3D, how many slices (around the center slice) do you
% want to correct for each subvolume?
n_slice = 5;
if mod(n_slice,2) == 0 
% If the number of slices is an even number, it's a bit inconvenient. It's 
% best to be an odd number so that the numbers of slices above and below the
% center depth of the subvolume are the same
    n_slice = n_slice+1;
end

% That's all one need to put in in order to reconstruct the fully corrected
% image. 

%% I. System coordinates, zoning, Zernike matrices

fprintf("System coordinates, zoning, Zernike matrices.\n")

% 1. Dealing with the full image
% 1.1. Below is the lists of coordinate x, y, z when running the computation
% and x_im, y_im, z_im when reconstructing the image
x = single(lx_offset-lx_im/2+dx/2:dx:lx_offset+lx_im/2); Nx = length(x);
y = single(ly_offset-ly_im/2+dx/2:dx:ly_offset+ly_im/2); Ny = length(y);
z = single(lz_offset-lz_im/2+dz/2:dz:lz_offset+lz_im/2); Nz = length(z);
x_im = single(lx_offset-lx_im/2+dx_im/2:dx_im:lx_offset+lx_im/2); Nx_im = length(x_im);
y_im = single(ly_offset-ly_im/2+dx_im/2:dx_im:ly_offset+ly_im/2); Ny_im = length(y_im);
z_im = single(lz_offset-lz_im/2+dz_im/2:dz_im:lz_offset+lz_im/2); Nz_im = length(z_im);

% 1.2. Convert to uncorrected coordinate
n0 = 1; n1 = coef_n1(1); n2 = coef_n2(1);
% Average
term1_ave_freq = 0; term2_ave_freq = 0; n_freq = length(list_k0); N = length(k);
for i_freq = 1:n_freq
    k0 = list_k0(i_freq);
    term1 = 0; term2 = 0; n_angle = 0;
    for ii = 1:N
        kx_ii = k(ii,1); ky_ii = k(ii,2);
        kt = sqrt(kx_ii^2+ky_ii^2);
        
        if kt/k0 <= NA
            sin_theta_air = kt/k0;
            theta_air = asin(sin_theta_air);
            theta_1 = asin(sin(theta_air)/n1);
            theta_2 = asin(sin(theta_air)/n2);
            if sin_theta_air ~= 0
                term1 = term1+tan(theta_air)/tan(theta_2);
                term2 = term2+(tan(theta_air)-tan(theta_1))/tan(theta_2);
            else
                term1 = term1+1;
                term2 = term2+1;
            end
            n_angle = n_angle+1;
        end
    end
    
    term1_ave = term1/n_angle; 
    term2_ave = term2/n_angle;
    term1_ave_freq = term1_ave_freq+term1_ave;
    term2_ave_freq = term2_ave_freq+term2_ave;
end

term1_ave_freq = term1_ave_freq/n_freq;
term2_ave_freq = term2_ave_freq/n_freq;

z_air_min = (min(z,[],'all')-h1*term2_ave_freq)/term1_ave_freq;
z_air_max = (max(z,[],'all')-h1*term2_ave_freq)/term1_ave_freq;
z_air = linspace(z_air_min,z_air_max,Nz);
z_air_im = linspace(z_air_min,z_air_max,Nz_im);
z_apparent = z_air-z_ref_air;
z_apparent_im = z_air_im-z_ref_air;

% Now, build the 3D image and scan for the depth with the maximum
% reflectance
[X,Y,Z_a] = meshgrid(single(x),single(y),single(z_apparent)); 
psi_3D = zeros(Nx,Ny,Nz);
for i_freq = 1:n_freq
    % Load the reflection matrix at this frequency
    fprintf("Working with frequency"+i_freq+"\n")
    k0 = list_k0(i_freq); freq = 3e2*k0/2/pi;

    % Load the reflection matrix at this frequency
    r = load(""+directory_r+""+prefix+""+i_freq+".mat").r_pad;

    % In the sample
    kz_in = sqrt(k0^2-kx_in.^2-ky_in.^2);
    kz_out = -sqrt(k0^2-kx_out.^2-ky_out.^2);
    fz0 = kz_out-kz_in;

    % Build the image at this frequency and scale it with the laser
    % amplitude
    psi_s = finufft3d3(fy(:),fx(:),fz0(:),r(:),1,1e-2,Y(:),X(:),Z_a(:))/list_amp(i_freq);

    psi_3D = psi_3D+reshape(psi_s,Ny,Nx,Nz);
end
I_3D = abs(psi_3D).^2;

% Scan for the depth with the highest reflectance
list_M = []; % list of FOM of 2D image at each depth
for ii = 1:Nz % Scan each depth and find the FOM of the 2D images
    I_2D = I_3D(:,:,ii);
    M = sum(I_2D,'all');
    list_M = [list_M M];    
end
% Search for target depth
z_max = z_a(list_M == max(list_M,[],'all'));

% The z coordinate for vrm image is
z_vrm = z_apparent-z_max;
z_vrm_im = z_apparent_im-z_max;
dz_vrm = abs(z_vrm(2)-z_vrm(1));
dz_vrm_im = abs(z_vrm_im(2)-z_vrm_im(1));

save(""+directory_save+"./z_vrm_"+im_case+".mat",'z_vrm')

% 1.3. Zoning in z for 3D case
if im_case == "3D"
    [list_z_vrm, list_z_vrm_im] = zoning_z(z_vrm,z_vrm_im,dz_vrm,dz_vrm_im,n_sub,overlap2); % list_z is a cell containing the z of each
end

%% II. Compensate for dispersion and build time-gated reflection matrices
if im_case == "2D"
    % Correct dispersion with VRM
    fprintf("Correcting dispersion. \n")
    dispersion_VRM(directory_r,prefix,list_k0,im_case,directory_save);
    % Find target depth (it will probabaly be z_vrm = 0 anyway but recheck
    % just to make sure)
    [z_target] = find_z_target(x,y,z_vrm,list_k0,k,im_case,directory_save);
    % Build time-gated r matrix
    fprintf("Building time-gated matrices. \n")
    build_r_z_2D_vrm(z_target,list_amp,list_k0,directory_save,k);
elseif im_case == "3D"
    fprintf("Correcting dispersion. \n")
    dispersion_VRM(directory_r,prefix,list_k0,im_case,directory_save);
    fprintf("Building time-gated matrices. \n")    
    build_r_z_3D_vrm(list_z_vrm,dz_vrm,list_amp,list_k0,directory_save,n_sub,n_slice,k);
end

%% III. Wavefront correction
for subvolume_id = 1:n_sub
    if im_case == "2D"
        r_z = single(load(""+directory_save+"r_z_2D_vrm.mat").r_z);
        wavefront_correction_vrm_2D(x,y,k,r_z,subvolume_id,directory_save);
    elseif im_case == "3D"
        r_z = load(""+directory_save+"r_z_3D_vrm_"+subvolume_id+".mat").rz;
        wavefront_correction_vrm_3D(x,y,k,rz,n_slice,subvolume_id,directory_save);
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
    [I] = reconstruct3D_vrm(list_z_vrm_im,x_im,y_im,z_vrm_im,overlap2,k,list_k0,phase_list,list_amp,directory_r,prefix,directory_save); 
    % Choose the z (depth, in pixel coordinate) to display 
    % Show enface image: choose the z (depth, in pixel coordinate) to display 
    n_z_show = 8;
    z_display = round(length(z_im)/(n_z_show+1)*[1:1:n_z_show]);
    for ii = z_display(4)
        I_z = I(:,:,ii);
        I_z = I_z/max(I_z,[],'all');
        figure
        imagesc(fliplr(I_z));
        colormap('hot')
        axis image
        caxis([0 0.5])
        set(gca,'Visible','off')
    end
    
    % Show cross section image, choose the constant x to display
    n_x_show = 3;
    x_display = round(length(x_im)/(n_x_show+1)*[1:1:n_x_show]);
    for ii = x_display
        I_x = I(ii,:,:);
        I_x = reshape(I_x(:),Ny_im,Nz_im); 
        I_x = I_x/max(I_x,[],'all');
        figure
        imagesc(fliplr(I_x));
        colormap('hot')
        axis image
        caxis([0 0.5])
        set(gca,'Visible','off')
    end
end

