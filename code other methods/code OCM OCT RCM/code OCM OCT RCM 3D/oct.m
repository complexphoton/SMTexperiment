clear
clc
close all

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
NA_reduce = single(0.1);

% Angle, kx_center = 0; ky_center = k0*sin(alpha)
alpha = pi/18; 

% Dispersion compensation phase is the same as SMT
phase = load("/media/minh/WD_BLACK/SMT 3D Tb3 1350 and 1150/Tb3 1150/Intermediate result/phase_list_3D_single.mat").phase_list;

% Focal plane
zf = z_im(round(Nz_im/2));
list_amp = single(load(""+directory_r+"/list_amp.mat").list_amp_ave);

%% Build high resolution image
[X_im,Y_im] = meshgrid(x_im,y_im); 
psi_oct = zeros(Ny_im,Nx_im,Nz_im);
freqc = 3e2*list_k0(round(n_freq/2))/2/pi;
for i_freq = 1:n_freq
    i_freq
    k0 = list_k0(i_freq); freq = 3e2*k0/2/pi; dfreq = freq - freqc; 
    omega = 2*pi*freq;    
    
    % The refractive indices of glass and the medium at the
    % corresponding frequency
    n1 = coef_n1(1) + coef_n1(2)*dfreq + coef_n1(3)*dfreq.^2 + coef_n1(4)*dfreq.^3 + coef_n1(5)*dfreq.^4;
    n_media = coef_n2(1) + coef_n2(2)*dfreq + coef_n2(3)*dfreq.^2 + coef_n2(4)*dfreq.^3 + coef_n2(5)*dfreq.^4;

    % Find kz and Fourier components fz 
    % In air
    kz_in0 = (sqrt(k0^2-kx.^2-ky.^2)).';
    kz_out0 = (-sqrt(k0^2-kx.^2-ky.^2));
    fz0 = (kz_out0-kz_in0);
    % In glass
    kz_in1 = (sqrt(k0^2*n1^2-kx.^2-ky.^2)).';
    kz_out1 = (-sqrt(k0^2*n1^2-kx.^2-ky.^2));
    fz1 = (kz_out1-kz_in1);   
    % In the sample
    kz_in = (sqrt(k0^2*n_media^2-kx.^2-ky.^2)).';
    kz_out = (-sqrt(k0^2*n_media^2-kx.^2-ky.^2));
    fz = (kz_out-kz_in);
    
    r = single(load(""+directory_r+"r_pad_"+i_freq+".mat").r_pad)*exp(-1i*phase(i_freq));
    
    % Index mismatch correction
    r = r.*exp(1i.*mod(-fz0*(z_ref_air+h1)+fz1*h1,2*pi));
    
    % Remove the out-of-NA k components
    k_mask = zeros(N,1);
    k_mask(kx.^2+(ky-k0*sin(alpha)).^2 <= k0^2*NA_reduce^2) = 1;
    K_mask = k_mask.*k_mask.';
    
    r = r.*K_mask;
    
    % Time domain
    time = (z_im-zf)/(300/n_media);
    
    psi_oct_single_freq = exp(-1i*omega*time*2/cos(alpha)).*finufft2d3(fy(:),fx(:),r(:).*exp(1i*fz(:)*zf),1,1e-2,Y_im(:),X_im(:))/list_amp(i_freq);
    psi_oct = psi_oct+reshape(psi_oct_single_freq,Ny_im,Nx_im,Nz_im);
end

I_oct = abs(psi_oct).^2;
save(""+directory_save+"OCT.mat",'I_oct')

%% Show image
n = 8;
for ii = 1:n
    I_xy = I_oct(:,:,round(ii*Nz_im/(n+1)));
    I_yz = I_oct(:,round(ii*Nx_im/(n+1)),:); I_yz = I_yz(:); I_yz = reshape(I_yz,Ny_im,Nz_im);

    figure(ii)
    imagesc(I_xy)
    axis image
    colormap('hot')
   % caxis([0 0.5]*max(I_xy,[],'all'))
    
    figure(ii+n)
    imagesc(I_yz)
    axis image
    colormap('hot')
    caxis([0 1]*max(I_yz,[],'all'))
end

