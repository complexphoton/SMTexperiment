%% 0. Settings
% Pixel size
dx = single(0.5); dx_im = single(0.2);
dz = single(0.5); dz_im = single(0.2);

% Size of the image and depths of the target in SMT image
lx_offset = single(25); lx_im = single(50);
ly_offset = single(25); ly_im = single(50);
lz_offset = single(1525); lz_im = single(110);

% The x,y,z axes
z_im = single(lz_offset-lz_im/2+dz_im/2:dz_im:lz_offset+lz_im/2);  Nz_im = length(z_im);
x_im = single(lx_offset-lx_im/2+dx_im/2:dx_im:lx_offset+lx_im/2); Nx_im = length(x_im);
y_im = single(ly_offset-ly_im/2+dx_im/2:dx_im:ly_offset+ly_im/2); Ny_im = length(y_im);

% Directory to load the reflection matrices from
directory_r = ["/media/minh/WD_BLACK/SMT 3D Tb3 1350 and 1150/Tb3 1150/"];

% List of frequencies
list_k0 = single(load(""+directory_r+"list_k0.mat").list_k0);
k0_max = max(list_k0,[],'all'); n_freq = length(list_k0);

% Refractive indices of glass and the medium
coef_n1 = single(load(""+directory_r+"coef_n1.mat").coef_n1);
coef_n2 = single(load(""+directory_r+"coef_n2_3D.mat").coef_n2);

% Glass thickness and reference depth in air
h1 = single(154); z_ref_air = single(1150)-h1;

prefix = "r_pad_";

% Imaging scenario, either 2D or 3D
im_case = "3D";

% List of laser amplitudes at different frequencies. Normalize each reflection matrix with the corresponding amplitude
list_amp = single(load(""+directory_r+"/list_amp.mat").list_amp_ave);

% Save the results to
directory_save = ["/media/minh/WD_BLACK/SMT 3D Tb3 1350 and 1150/Tb3 1150/Intermediate result/other methods/"];

% List of transverse momentum
k = single(load(""+directory_r+"k.mat").k);
kx = k(:,1); ky = k(:,2); 
fx = kx-kx.'; fy = ky-ky.';
N = length(k);

% Numerical aperture
NA = single(0.5);
list_amp = single(load(""+directory_r+"/list_amp.mat").list_amp_ave);

%% Build image
i_freq = round(n_freq/2);
[X_im,Y_im,Z_im] = meshgrid(x_im,y_im,z_im);
r = single(load(""+directory_r+""+prefix+""+i_freq+".mat").r_pad)/list_amp(i_freq);
n_medium = coef_n1(1); n_glass = coef_n2(1);
kz = single(-sqrt(k0^2*n_medium^2-kx.^2-ky.^2);
kz_glass = single(-sqrt(k0^2*n_glass^2-kx.^2-ky.^2);
kz0 = single(-sqrt(k0^2-kx.^2-ky.^2));
fz = -kz-kz.';
fz_glass = -kz_glass-kz_glass.';
fz0 = -kz0-kz0.';
fx = kx-kx.'; fy = ky-ky.';

% Correct index-mismatch
r = r.*exp(-1i*mod(fz0*(z_ref_air+h1)+fz_glass*h1,2*pi));

I_rcm = abs(finufft3d3(fy(:),fx(:),fz(:),r(:),1,1e-2,Y_im(:),X_im(:),Z_im(:))).^2;
I_rcm = reshape(I_rcm,Ny_im,Nx_im,Nz_im);

%% Show image
for ii = 1:7
    I_xy = I_rcm(:,:,round(ii*Nz_im/8));
    I_yz = I_rcm(:,round(ii*Nx_im/8),:); I_yz = I_yz(:); I_yz = reshape(I_yz,Ny_im,Nz_im);

    figure
    imagesc(I_xy)
    axis image
    colormap('hot')

    figure
    imagesc(I_yz)
    axis image
    colormap('hot')
end

save(""+directory_save+"RCM.mat",'I_rcm')
