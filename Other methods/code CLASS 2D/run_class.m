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
lz_offset = 990; lz_im = 50;

% 2. Provide the path to input data. We recommend saving all required input 
% data at the same directory, including the reflection matrices, k list, k0 
% list, dispersion coefficients of glass and sample, etc
directory_r = ["/media/minh/WD_BLACK/smt exp data 2d/"];

% 3. List of k0 at all frequencies.
list_k0 = load(""+directory_r+"list_k0_2D.mat").list_k0;
k0_max = max(list_k0,[],'all');

% 3. Dispersion coefficients of glass is coef_n1 and sample is coef_n2
coef_n1 = load(""+directory_r+"coef_n1.mat").coef_n1;
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
im_case = "2D";

% If the case is "3D", we need to specify the number of subvolumes in the z
% axis. If "2D" then the number of subvolumes is 1
n_sub = 1;

% 9. The image will go through one correction of the full image, then
% n_divison division steps
n_division = 3;

% 11. List of laser amplitude at all frequencies
list_amp = load(""+directory_r+"/list_amp_2D.mat").list_amp;

% 12. Specify a saving directory to store the intermediate variables which
% can be very heavy. We will store those heavy intermediate variables after
% obtaining them, then clear them.
directory_save = ["/media/minh/WD_BLACK/smt exp data 2d/Intermediate data/class/"];

% 13. Specify the list of k_parallel_in and k_parallel_out (kin/out_x/y).
% When saving these k, please save them as a nx2 matrix where the first
% column is kx, the 2nd column is ky. For each constant value of kx, the
% values of ky increases from negative to positive. Arranging like that
% helps us easily find the spacing in k space. The reflection matrices are
% also arranged in that order.
k = single(load(""+directory_r+"k_2D.mat").k);

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

% That's all one need to put in in order to reconstruct the fully corrected
% image. 

%% I. System coordinates, zoning, Zernike matrices

fprintf("System coordinates, zoning, Zernike matrices.\n")

% 1. Dealing with the full image
% 1.1. Below is the lists of coordinate x, y, z when running the computation
% and x_im, y_im, z_im when reconstructing the image
x = lx_offset-lx_im/2+dx/2:dx:lx_offset+lx_im/2; Nx = length(x);
y = ly_offset-ly_im/2+dx/2:dx:ly_offset+ly_im/2; Ny = length(y);
z = lz_offset-lz_im/2+dz/2:dz:lz_offset+lz_im/2; Nz = length(z);
x_im = lx_offset-lx_im/2+dx_im/2:dx_im:lx_offset+lx_im/2; Nx_im = length(x_im);
y_im = ly_offset-ly_im/2+dx_im/2:dx_im:ly_offset+ly_im/2; Ny_im = length(y_im);
z_im = lz_offset-lz_im/2+dz_im/2:dz_im:lz_offset+lz_im/2; Nz_im = length(z_im);

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
z_class = z_air-z_ref_air;
z_class_im = z_air_im-z_ref_air;
dz_class = abs(z_air(2)-z_air(1));
dz_class_im = abs(z_air_im(2)-z_air_im(1));
save(""+directory_save+"./z_class_"+im_case+".mat",'z_class')

%% 2. Divide the image into zones, obtaining the list of x,y,z for each
% division step. Each list will be stored in a 2D cell. Each column of the
% cell belongs to the same division step. 
% Initialize the cells
list_x = cell(1,1); list_y = cell(1,1);
list_x_im = cell(1,1); list_y_im = cell(1,1);
list_x{1,1} = x; list_y{1,1} = y;
list_x_im{1,1} = x_im; list_y_im{1,1} = y_im;
for division_step = 1:n_division
    [x_step,y_step,x_im_step,y_im_step] = zoning(dx,dx_im,list_x,list_y,list_x_im,list_y_im,division_step,overlap2); 
    % After this step, we obtain the list of x,y,z of the corresponding
    % division step as column cell. We have to put those cells in the
    % correct position in the list_x,y,z cells. The number of zones is
    % 2^(division_step*2)
    for ii = 1:2^(division_step*2)
        list_x{ii,division_step+1} = x_step{ii};
        list_y{ii,division_step+1} = y_step{ii};
        list_x_im{ii,division_step+1} = x_im_step{ii};
        list_y_im{ii,division_step+1} = y_im_step{ii};
    end
end
% Meanwhile, if we do 3D imaging the list of z in each subvolumes will be
% found by 
if im_case == "3D"
    [list_z_class, list_z_class_im] = zoning_z(z_class,z_class_im,dz_class,dz_class_im,n_sub,overlap2); % list_z is a cell containing the z of each
end

%% II. Correct dispersion and build time-gated reflection matrices
begin_compensation = tic;
if im_case == "2D"
    fprintf("Correcting dispersion. \n")
    [phase_list, z_target] = dispersion_class(x,y,z_class,list_amp,directory_r,prefix,list_k0,k,im_case,directory_save);
    fprintf("Building time-gated matrices. \n")
    build_r_z_2D_class(z_target,list_amp,phase_list,directory_r,prefix,list_k0,directory_save,k);
elseif im_case == "3D"
    fprintf("Correcting dispersion. \n")
    [phase_list, ~] = dispersion_class(x,y,z_class,list_amp,directory_r,prefix,list_k0,k,im_case,directory_save);
    fprintf("Building time-gated matrices. \n")    
    build_r_z_3D_class(list_z_class,list_amp,phase_list,directory_r,prefix,list_k0,directory_save,k);
end
toc(begin_compensation)

%% III. Wavefront correction
for subvolume_id = 1:n_sub
    if im_case == "2D"
        r_z = single(load(""+directory_save+"r_z_2D_class.mat").r_z);
    elseif im_case == "3D"
        r_z = load(""+directory_save+"r_z_3D_class_"+subvolume_id+".mat").r_z;
    end
    % We deal with each zone, there are 2^(division_step*2) zones
    n_zone = 2^(n_division*2);
    for zone_id = 1:n_zone
        % The corresponding x and y of the zone is
        x_zone = list_x{zone_id,n_division+1};
        y_zone = list_y{zone_id,n_division+1};
        % Correct the wavefront of the zone
        wavefront_correction_class(x_zone,y_zone,k,r_z,n_division,zone_id,subvolume_id,directory_save)
    end
end

%% IV. Reconstruct image
% Build mask to compensate vignetting
[mask] = build_mask(x_im,y_im,wx0,wy0,x0,y0,rho);

fprintf("Reconstructing the image. \n")
if im_case == "2D"
    
    [I] = reconstruct2D_class(list_x_im,list_y_im,x_im,y_im,overlap2,k,n_division,r_z,directory_save);
    
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
    [I] = reconstruct3D_class(list_x_im,list_y_im,list_z_class_im,x_im,y_im,z_class_im,overlap2,k,list_k0,phase_list,list_amp,n_division,directory_r,prefix,directory_save); 
    % Choose the z (depth, in pixel coordinate) to display 
    % Show enface image: choose the z (depth, in pixel coordinate) to display 
    n_z_show = 8;
    z_display = round(length(z_im)/(n_z_show+1)*[1:1:n_z_show]);
    for ii = z_display
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
        I_x = reshape(I_x(:),Ny_im,length(I(1,1,:))); 
        I_x = I_x/max(I_x,[],'all');
        figure
        imagesc(fliplr(I_x));
        colormap('hot')
        axis image
        caxis([0 0.5])
        set(gca,'Visible','off')
    end
end

