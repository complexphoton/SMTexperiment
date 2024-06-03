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
lz_offset = single(995); lz_im = single(50);

% 2. Provide the path to input data. We recommend saving all required input 
% data at the same directory, including the reflection matrices, k list, k0 
% list, dispersion coefficients of glass and sample, etc
directory_r = ["/media/minh/WD_BLACK/smt exp data 2d/"];

% 3. List of k0 at all frequencies.
list_k0 = load(""+directory_r+"list_k0_2D.mat").list_k0;
n_freq = length(list_k0);
k0_max = max(list_k0,[],'all');

% 4. Dispersion coefficients of glass is coef_n1 and sample is coef_n2
coef_n1 = single(load(""+directory_r+"coef_n1.mat").coef_n1);
coef_n2 = load(""+directory_r+"coef_n2_2D.mat").coef_n2;

% 4. Thickness of the coverslip h1 and mirror depth z_mirror (micron)
h1 = 141.5; z_ref_air = 820-h1;

% 5. In our case, ehe reflection matrices are stored in "/directory_r/" as
% ("/prefix"-frequency_id-".mat").r_pad. So we also need to specify 
% the prefix. Also, for the best convenience, in the future, before using
% this code, change the name of the file containing the reflection matrices
% into that format
prefix = "r_pad_2D_freq_";

% 6. Specify the pixel size. % Pixel size dx = dy, dz when running the
% optimization. When reconstructing the image, we can use a different pixel
% size dx_im = dy_im, 
dx = 0.7; dx_im = 0.2;
dz = 1; 

% 7. List of laser amplitude at all frequencies
list_amp = load(""+directory_r+"/list_amp_2D.mat").list_amp;

% 8. Specify a saving directory to store the intermediate variables which
% can be very heavy. We will store those heavy intermediate variables after
% obtaining them, then clear them.
directory_save = ["/media/minh/WD_BLACK/smt exp data 2d/Intermediate data VRM 2D/"];

% 9. Specify the list of k_parallel_in and k_parallel_out (kin/out_x/y).
% When saving these k, please save them as a nx2 matrix where the first
% column is kx, the 2nd column is ky. For each constant value of kx, the
% values of ky increases from negative to positive. Arranging like that
% helps us easily find the spacing in k space. The reflection matrices are
% also arranged in that order.
k = single(load(""+directory_r+"k_2D.mat").k); N = length(k); kx = k(:,1); ky = k(:,2);

% 10. NA of the system
NA = 0.5;

% 11. % The outer part of the image is always darker than the central part
% due to vignetting. To compensate, we multiply the image with a
% Gaussian mask with the following parameters
x0 = 30; y0 = 22; % Central of the mask
wx0 = 40; wy0 = 22; % Width
rho = 0;
% The mask will be:
% [exp(-1/(2*(1-rho^2))*(((gridX-x0)/wx0).^2+((gridY-y0)/wy0).^2)-...
% 2*rho.*(gridX-x0)/wx0.*(gridY-y0)/wy0)'.^2

% 12. Truncation range in z for 2D imaging case
range_z = 9;
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

% Convert z coordinate from medium to air
n0 = 1; n1 = coef_n1(1); n2 = coef_n2(1);
% Average
term1_ave_freq = 0; term2_ave_freq = 0;
for i_freq = 1:n_freq
    k0 = list_k0(i_freq);
    term1 = 0; term2 = 0; n_angle = 0;
    for ii = 1:N
        kx_ii = kx(ii); ky_ii = ky(ii);
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

z_prime = (z-h1*term2_ave_freq)/term1_ave_freq-z_ref_air;

%% 1.3. Scan the z axis to find the depth with
% the maximum intensity to set as the reference plane of VRM
kx_out = k(:,1); ky_out = k(:,2); kx_in = kx_out.'; ky_in = ky_out.';
fx = kx_out-kx_in; fy = ky_out-ky_in;

[X,Y,Z_prime] = meshgrid(x,y,z_prime);
[X,Y,Z] = meshgrid(x,y,z);

psi = single(zeros(Nx,Ny,Nz));
for i_freq = 1:n_freq
    r = single(load(""+directory_r+""+prefix+""+i_freq+".mat").r_pad);

    fprintf("Building uncorrected image at frequency "+i_freq+"\n")
    k0 = list_k0(i_freq); omega = k0*300;

    % Find kz and Fourier components fz 
    % In air
    kz_in0 = (sqrt(k0^2-kx_in.^2-ky_in.^2));
    kz_out0 = (-sqrt(k0^2-kx_out.^2-ky_out.^2));
    fz0 = (kz_out0-kz_in0);

    % Build the image at this frequency and scale it with the laser
    % amplitude
    psi_s = finufft3d3(fy(:),fx(:),fz0(:),r(:),1,1e-2,Y(:),X(:),Z_prime(:))/list_amp(i_freq);

    % Compensate for time gate mismatch 
    time = (n2*Z-(Z_prime+z_ref_air+h1)+n1*h1)/300;
    
    psi = psi+exp(-2*1i*omega*time).*reshape(psi_s,Nx,Ny,Nz);
end

I = abs(psi).^2;
clear psi

list_M_z = [];

for ii = 1:Nz
    I_2D = I(:,:,ii);
    M = sum(I_2D,'all');
    list_M_z = [list_M_z M];
    figure(1)
    imagesc(fliplr(interp2(I_2D)))
    axis image
    colormap('hot')
end

z_prime_max = z_prime(list_M_z == max(list_M_z,[],'all'));
z_max = z(list_M_z == max(list_M_z,[],'all'));

% Image before doing anything
I_pre = I(:,:,list_M_z == max(list_M_z,[],'all'));

%% II. Do time truncation

z_truncate = z(z >= z_max-range_z/2 & z <= z_max+range_z/2);
z_truncate_min = z_truncate(1); z_truncate_max = z_truncate(end);

fprintf("Do time domain truncation.\n")
t_truncation_no_index_mismatch_correction(directory_r,prefix,list_k0,z_truncate_min,z_truncate_max,z_ref_air,z_prime_max,k,coef_n2,coef_n1,h1);

%% III. Compensate for dispersion and build time-gated reflection matrices
% Build 3D image
psi_after_truncation = single(zeros(Nx,Ny,Nz));
for i_freq = 1:n_freq
    fprintf("Build pre-vrm image at frequency "+i_freq+".\n")

    k0 = list_k0(i_freq); omega = k0*300;
    kz_in = (sqrt(k0^2-kx_in.^2-ky_in.^2));
    kz_out = (-sqrt(k0^2-kx_out.^2-ky_out.^2));
    fz0 = (kz_out-kz_in);

    r = load(""+directory_r+"r_truncated_"+i_freq+".mat").r;

    time = (n2*Z-(Z_prime+z_ref_air+h1)+n1*h1)/300;

    psi_after_truncation = psi_after_truncation+exp(-2*1i*omega*time).*reshape(finufft3d3(fx(:),fy(:),fz0(:),r(:),1,1e-2,X(:),Y(:),Z_prime(:)),Nx,Ny,Nz);

end

I_after_truncation = abs(psi_after_truncation).^2;
clear psi_after_truncation

list_M = [];
for ii = 1:Nz
    I_2D = I_after_truncation(:,:,ii);
    M = sum(I_2D,'all');
    list_M = [list_M M];
end

% The depth with the maximum intensity will be comes the new reference
% depth when doing VRM
z_max = z(list_M == max(list_M,[],'all'));
z_prime_max = z_prime(list_M == max(list_M,[],'all'));

I_max = I_after_truncation(:,:,list_M == max(list_M,[],'all'));

figure(2)
imagesc(fliplr(interp2(I_max,3)))
axis image
colormap('hot')

% Do VRM dispersion
fprintf("Correcting dispersion. \n")
dispersion_VRM_no_index_mismatch_correction(directory_r,directory_save,z_max,z_prime_max,z_ref_air,list_k0,k,coef_n1,coef_n2,h1);

% Build time-gated matrix at z_air_max
fprintf("Building time-gated matrices. \n")
r_z = single(zeros(N,N));
for i_freq = 1:n_freq
    r_vrm = load(""+directory_save+"r_vrm_"+i_freq+".mat").r;
    r_z = r_z+r_vrm;
end

%% IV. Wavefront correction
wavefront_correction_vrm_2D(x,y,k,r_z,directory_save);
    
%% V. Reconstruct image
% Build mask to compensate vignetting
[mask] = build_mask(x_im,y_im,wx0,wy0,x0,y0,rho);

fprintf("Reconstructing the image. \n")
[I] = reconstruct2D_vrm(x_im,y_im,k,directory_r,z_max,z_prime_max,z_ref_air,h1,coef_n1,coef_n2,list_k0,directory_save,prefix);

I_show = I.*mask;
I_show = fliplr(I_show);
% Display image
figure
imagesc(I_show);
cmap_heat
caxis([0 3e19])
axis image
colorbar
cbh = colorbar ; %Create Colorbar
cbh.Ticks = linspace(0, 0, 0) ; %Create 8 ticks from zero to 1
set(gca,'Visible','off')

function [] = cmap_heat
cmap = hot(256);
scale_factor = 1.3;
if scale_factor ~= 1
    np = size(cmap, 1);
    x1_cmap = linspace(-1,1,np);
    x2_cmap = linspace(-1,1,np);
    x2_cmap = sign(x2_cmap).*(abs(x2_cmap).^(scale_factor));
    cmap_original = cmap;
    cmap = zeros(np, 3);
    for i_cmap=1:3
        cmap(:,i_cmap) = interp1(x1_cmap,cmap_original(:,i_cmap),x2_cmap);
    end
end
colormap(cmap);
end
