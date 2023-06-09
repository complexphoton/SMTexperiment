%% SMT image reconstruction code from reflection matrices in angular space 
% Required input info: System size (x, y, z), list of k0, kin, kout,
% dispersion coefficients of glass coverslip and sample, coverslip's
% thickness, depth of the reference mirror, list of laser amplitudes at all
% frequencies, NA of the system.

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
lz_offset = 975; lz_im = 10;

% 2. Provide the path to input data. We recommend saving all required input 
% data at the same directory, including the reflection matrices, k list, k0 
% list, dispersion coefficients of glass and sample, etc
directory_r = ["/media/minh/Data/Imaging scatter media/Aberration correction project/Yiwen's paper data/2D/"];

% 3. List of k0 at all frequencies.
list_k0 = load(""+directory_r+"list_k0_2D.mat").list_k0;
k0_max = max(list_k0,[],'all');

% 3. Dispersion coefficients of glass is coef_n1 and sample is coef_n2
coef_n1 = load(""+directory_r+"coef_n1.mat").coef_n1;
coef_n2 = load(""+directory_r+"coef_n2_2D.mat").coef_n2;

% 4. Thickness of the coverslip h1 and mirror depth z_mirror (micron)
h1 = 141.5; z_mirror = 820;

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
dz = 0.5; dz_im = 0.2;

% 7. The image is divided into smaller zones during the wavefront
% correction procedure. Those zones are overlapped. Specify the
% half-of-overlapping-area here (micron). The overlapping region of 
% two adjacent zones is 2*overlap2
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
n_division = 4;

% 10. At each optimization step, we optimize a certain number of Zernike
% polynomials. After each optimization step, we increase the number of
% Zernike modes. Specify the number of Zernike radial orders to start from
% and the number of incremental Zernike radial orders
rad_order_start = 11;
rad_order_inc = 5;

% 11. List of laser amplitude at all frequencies
list_amp = load(""+directory_r+"/list_amp_2D.mat").list_amp;

% 12. Specify a saving directory to store the intermediate variables which
% can be very heavy. We will store those heavy intermediate variables after
% obtaining them, then clear them.
directory_save = ["/media/minh/Data/Imaging scatter media/Aberration correction project/Yiwen's paper data/2D/Intermediate result/"];

% 13. Specify the list of k_parallel_in and k_parallel_out (kin/out_x/y).
% When saving these k, please save them as a nx2 matrix where the first
% column is kx, the 2nd column is ky. For each constant value of kx, the
% values of ky increases from negative to positive. Arranging like that
% helps us easily find the spacing in k space. The reflection matrices are
% also arranged in that order.
k = load(""+directory_r+"k_2D.mat").k;

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
x = lx_offset-lx_im/2+dx/2:dx:lx_offset+lx_im/2;
y = ly_offset-ly_im/2+dx/2:dx:ly_offset+ly_im/2;
z = lz_offset-lz_im/2+dz/2:dz:lz_offset+lz_im/2;
x_im = lx_offset-lx_im/2+dx_im/2:dx_im:lx_offset+lx_im/2;
y_im = ly_offset-ly_im/2+dx_im/2:dx_im:ly_offset+ly_im/2;
z_im = lz_offset-lz_im/2+dz_im/2:dz_im:lz_offset+lz_im/2;

% 1.2. Build Zernike matrix for the full image
division_step = 0;
rad_order = rad_order_start+rad_order_inc*division_step; % Number of radial orders
build_zernike(k0_max,k,NA,division_step,rad_order,directory_save,lx_im);

% 2. Divide the image into zones, obtaining the list of x,y,z for each
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
    % In addition to zoning, for each zone division, we use a different k
    % space and a different number of Zernike polynomials. So, we have to
    % build a matrix that stores the corresponding Zernike polynomials.
    % Each zone division will correspond to two Zernike matrices, one is
    % for optimization, stored as "/directory_save/Z_"+division_step+".mat"
    % and the other is for image reconstruction, stored as
    % "directory_save/Z_"+division_step+"_re.mat". The list of k in the new
    % k space will be saved as "directory_save/k_"+division_step+".mat".
    rad_order = rad_order_start+rad_order_inc*division_step; % Number of radial orders
    l_zone = max(list_x{1,division_step+1},[],'all')-min(list_x{1,division_step+1},[],'all');
    build_zernike(k0_max,k,NA,division_step,rad_order,directory_save,l_zone);
    % In this function, k0_max*NA is the maximum k_parallel, which is used
    % to normalize the k space, k is the list of original k_parallel, which
    % is used to optimize the full 2D image, as well as image
    % reconstruction. For division_step > 0 (dealing with zones instead of
    % full image), a new k space is built for optimization purpose. This
    % new k space, k_zone, has a bigger spacing, which we choose as
    % dk*(sqrt(2)^{division_step}) where dk is the spacing of the original
    % k space. One can access the build_zernike function to change this
    % number into whatever they like. But keep in mind that the k space is
    % a circle, and if dk is too large, a lot of k components near the edge
    % of the circle will be missed.
end
% Meanwhile, if we do 3D imaging the list of z in each subvolumes will be
% found by 
if im_case == "3D"
    [list_z, list_z_im] = zoning_z(z,dz,dz_im,n_sub,overlap2); % list_z is a cell containing the z of each
end

%% II. Compensate dispersion then build the time-gated matrices
begin_compensation = tic;
if im_case == "2D"
    fprintf("Correcting dispersion. \n")
    [phase_list, z_target] = dispersion(x,y,z,list_amp,directory_r,prefix,h1,z_mirror,coef_n1,coef_n2,list_k0,k,im_case,directory_save);
    fprintf("Building time-gated matrices. \n")
    build_r_z_2D(z_target,list_amp,phase_list,directory_r,prefix,h1,z_mirror,coef_n1,coef_n2,list_k0,directory_save,k);
elseif im_case == "3D"
    fprintf("Correcting dispersion. \n")
    [phase_list, ~] = dispersion(x,y,z,list_amp,directory_r,prefix,h1,z_mirror,coef_n1,coef_n2,list_k0,k,im_case,directory_save);
    frintf("Building time-gated matrices. \n")    
    build_r_z_3D(list_z,list_amp,phase_list,directory_r,prefix,h1,z_mirror,coef_n1,coef_n2,list_k0,directory_save,k);
end
% After this step, we obtain r_z, but it will be saved then cleared. The
% user can access the root script of these build_r_z functions to modify
% the names they want to save their r_z files as. By default, we let the
% 2D file be "r_z_2D.mat" and the 3D files at subvolume n as
% "r_z_3D_"+n+".mat".
toc(begin_compensation)

%% III. Correct the spatial wavefront
begin_opt = tic;
fprintf("Correcting wavefront. \n")
for subvolume_id = 1:n_sub
    % 3.1. Full 2D image correction
    % For each subvolumne, we first, correct for the full 2D image, the
    % result are the input and output Zernike coefficients, which will be
    % stored and cleared. 
    if im_case == "2D"
        r_z = load(""+directory_save+"r_z_2D.mat").r_z;
    elseif im_case == "3D"
        r_z = load(""+directory_save+"r_z_3D_"+subvolume_id+".mat").r_z;
    end
    % Variable division_step = 0 means we are not doing any zone division
    % yet. Variable zone_id = 1 since there is only one zone, which is the
    % full image itself.
    division_step = 0;
    zone_id = 1;
    % The corresponding Zernike matrix
    Z = load(""+directory_save+"Z_"+division_step+".mat").Z;
    % Then, we do wavefront correction. The function wavefront_correction
    % will do the following works: do basis conversion from k basis to
    % truncated spatial basis (if division_step > 0) then convert back,
    % optimize Zernike coefficients alternatively.
    wavefront_correction(x,y,r_z,division_step,n_division,Z,zone_id,subvolume_id,directory_save);
    % The input Zernike coefficients will be stored as the following format
    % "/directory_save/c_in_"+division_step+"_zone_"+zone_id+"_subvolume_"+subvolume_id+".mat".
    % The same applies to the output Zernike coefficients.
    % For the full image, by default, division_step = 0 and zone_id = 1.
    % For subsequent divisions, division_step = 1, zone_id = 1,2,3,4;
    % division_step = 2, zone_id = 1,2,...,16; division_step = 3, zone_id =
    % 1,2,...,64; etc. In addition to the Zernike coefficients, the updated
    % reflection matrices are also saved at the same directory as
    % "r_update_"+division_step+"_zone_"+zone_id+"_subvolume_"+subvolume_id+".mat"
    
    % 3.2. Correct the zones after zone division steps
    for division_step = 1:n_division
        % We deal with each zone, there are 2^(division_step*2) zones
        n_zone = 2^(division_step*2);
        % The Zernike matrix
        Z = load(""+directory_save+"Z_"+division_step+".mat").Z;
        for zone_id = 1:n_zone
            % Find the bigger zone. Since each zone will contain 4 smaller
            % zones, we have:
            zone_big = ceil(zone_id/4);
            % Load the updated r matrix of the corresponding bigger zone
            prev_step = division_step-1;
            r_big = load(""+directory_save+"r_update_"+prev_step+"_zone_"+zone_big+"_subvolume_"+subvolume_id+".mat").r_update;
            % The corresponding x and y of the zone is
            x_zone = list_x{zone_id,division_step+1};
            y_zone = list_y{zone_id,division_step+1};
            % Correct the wavefront of the zone
            wavefront_correction(x_zone,y_zone,r_big,division_step,n_division,Z,zone_id,subvolume_id,directory_save);
        end
    end
end
toc(begin_opt)

%% IV. Reconstruct the image
% Build mask to compensate vignetting
[mask] = build_mask(x_im,y_im,wx0,wy0,x0,y0,rho);

fprintf("Reconstructing the image. \n")
if im_case == "2D"
    [I] = reconstruct2D(list_x_im,list_y_im,x_im,y_im,overlap2,k,n_division,r_z,directory_save);
    I = fliplr(I.*mask);
    I = I/max(I,[],'all');
    % Display image
    figure
    imagesc(I);
    colormap('hot')
    axis image
    caxis([0 0.6])
    set(gca,'Visible','off')
elseif im_case == "3D"
    [I] = reconstruct3D(list_x_im,list_y_im,list_z_im,x_im,y_im,z_im,overlap2,k,list_k0,phase_list,list_amp,coef_n1,coef_n2,h1,z_mirror,n_division,directory_r,prefix,directory_save); 
    % Choose the z (depth, in pixel coordinate) to display 
    z_display = [round(length(z_im)/4) 2*round(length(z_im)/4) 3*round(length(z_im)/4)];
    for ii = z_display
        I_z = I(:,:,ii);
        I_z = I_z/max(I,[],'all');
        figure
        imagesc(fliplr(I_z.*mask));
        colormap('hot')
        axis image
        caxis([0 0.5])
        set(gca,'Visible','off')
    end
end

toc(begin_all)
