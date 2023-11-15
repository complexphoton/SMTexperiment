clear
clc
close all
%% 0. Settings
% Size of the image and depth of the target in SMT image
lx_offset = single(25); lx_im = single(50);
ly_offset = single(25); ly_im = single(50);
z_target = 977.5;

% Directory to load the reflection matrices from
directory_r = ["/media/minh/WD_BLACK/smt exp data 2d/"];

% List of frequencies
list_k0 = single(load(""+directory_r+"list_k0_2D_single.mat").list_k0);
k0_max = max(list_k0,[],'all'); n_freq = length(list_k0);

% Refractive index of the medium and glass
coef_n1 = single(load(""+directory_r+"coef_n1_single.mat").coef_n1);
coef_n2 = single(load(""+directory_r+"coef_n2_2D_single.mat").coef_n2);

% Glass thickness and reference depth in air
h1 = single(141.5); z_ref_air = single(820)-h1;

prefix = "r_pad_2D_freq_";

% Pixel size
dx = single(0.5); dx_im = single(0.2);
dz = single(0.5); dz_im = single(0.2);

% Imaging scenario, either 2D or 3D
im_case = "2D";

% List of laser amplitude, normalize the reflection matrix of each frequency with the corresponding laser amplitude
list_amp = single(load(""+directory_r+"/list_amp_2D_single.mat").list_amp);

% Save the results to this directory
directory_save = ["/media/minh/WD_BLACK/smt exp data 2d/Intermediate data/other methods/"];

% List of transverse momentum
k = single(load(""+directory_r+"k_2D_single.mat").k);
kx = k(:,1); ky = k(:,2); 
N = length(k);

% Numerical aperture
NA = single(0.5);

% Parameters to generate mask that compensates for vignetting. May be used or not
x0 = single(30); y0 = single(22); % Central of the mask
wx0 = single(40); wy0 = single(22); % Width
rho = single(0);

% x, y axes
x = single(lx_offset-lx_im/2+dx/2:dx:lx_offset+lx_im/2); Nx = length(x);
y = single(ly_offset-ly_im/2+dx/2:dx:ly_offset+ly_im/2); Ny = length(y);
[X,Y] = meshgrid(x,y);

%% 1. Convert z_target to z'
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

z_air = (z_target-h1*term2_ave_freq)/term1_ave_freq;
z_ocm = z_air-z_ref_air;

%% Correct dispersion at this depth
Psi = zeros(Ny*Nx,n_freq);
Omega = zeros(n_freq,1);
omega_c = 3e2*list_k0(round(n_freq/2));
for i_freq = 1:n_freq
    i_freq
    k0 = list_k0(i_freq);
    omega = 3e2*list_k0(round(i_freq));
    Omega(i_freq,1) = (omega-omega_c)/omega_c;
    
    r = single(load(""+directory_r+"r_pad_"+im_case+"_freq_"+i_freq+".mat").r_pad);
    kz_in = single(sqrt(k0^2-kx.^2-ky.^2)).';
    kz_out = single(-sqrt(k0^2-kx.^2-ky.^2));
    fz0 = (kz_out-kz_in);
    fx = kx-kx.'; fy = ky-ky.';
    
    psi_ocm = exp(-2*1i*omega*z_ocm/300)*finufft2d3(fy(:),fx(:),r(:),1,1e-2,Y(:),X(:));
    Psi(:,i_freq) = psi_ocm;
end

%% Scan a1
M_list = [];
for a1 = -4000:10:4000
    a1
    phase = Omega*a1;
    psi = Psi*exp(-1i*phase);
    I = abs(psi).^2;
    M = sum(I.*log(I),'all');
    M_list = [M_list [M;a1]];
end

a1_list = M_list(2,:);
M_list = M_list(1,:);
    
a1_max = a1_list(:,M_list == max(M_list,[],'all'));
phase = Omega*a1_max;
psi = Psi*exp(-1i*phase);

I = abs(psi).^2;
I = reshape(I,Ny,Nx);

figure(1)
imagesc(interp2(I,3))
axis image
colormap('hot')


figure(2)
plot(a1_list,M_list)
%% Build high resolution image
x_im = single(lx_offset-lx_im/2+dx_im/2:dx_im:lx_offset+lx_im/2); Nx_im = length(x_im);
y_im = single(ly_offset-ly_im/2+dx_im/2:dx_im:ly_offset+ly_im/2); Ny_im = length(y_im);
[X_im,Y_im] = meshgrid(x_im,y_im);
psi_ocm = zeros(Ny_im,Nx_im);
for i_freq = 1:n_freq
    i_freq
    k0 = list_k0(i_freq);
    omega = 3e2*list_k0(round(i_freq));
    Omega(i_freq,1) = (omega-omega_c)/omega_c;
    
    r = single(load(""+directory_r+"r_pad_"+im_case+"_freq_"+i_freq+".mat").r_pad).*exp(-1i*phase(i_freq));
    kz_in = single(sqrt(k0^2-kx.^2-ky.^2)).';
    kz_out = single(-sqrt(k0^2-kx.^2-ky.^2));
    fz0 = (kz_out-kz_in);
    fx = kx-kx.'; fy = ky-ky.';
    
    psi_ocm = psi_ocm + reshape(exp(-2*1i*omega*z_ocm/300)*finufft2d3(fy(:),fx(:),r(:),1,1e-2,Y_im(:),X_im(:)),Ny_im,Nx_im);
    
end
I_ocm = abs(psi_ocm).^2;
figure(1)
imagesc(fliplr(I_ocm));
axis image
colormap('hot')

save(""+directory_save+"OCM.mat",'I_ocm','phase','z_ocm')
