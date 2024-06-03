% This function builds the Zernike wavefront matrices for both optimization
% and reconstruction
function build_zernike(k0_max,k,NA,division_step,rad_order,directory_save,l_zone)
    %% 1. Build the Zernike matrix for reconstruction
    
    % 1.1. Build the normalized polar k space
    % The maximum k_parallel
    kt_max = k0_max*NA;
    n1 = length(k);
    % Convert k to normalized polar coordinate
    k_norm_pol = zeros(n1,2); % 1st column: radial coordinate, 2nd column: angular coordinate
    for ii = 1:n1
        % radial coordinate
        k_norm_pol(ii,1) = sqrt(k(ii,1).^2+k(ii,2).^2)/kt_max;
        sine_theta = k(ii,2)/sqrt(k(ii,1).^2+k(ii,2).^2); % Sine of the angular coordinate
        % Now, find the angular coordinate
        if k(ii,1) >= 0 && (k(ii,1) ~= 0 || k(ii,2) ~= 0)
            theta = asin(sine_theta);
        elseif k(ii,1) < 0 && (k(ii,1) ~= 0 || k(ii,2) ~= 0)
            theta = pi-asin(sine_theta);
        elseif k(ii,1)==0 && k(ii,2) == 0
            theta = 0;
        end
        % after this step, theta is from -pi/2 to 3*pi/2. wrap it to pi
        theta = wrapToPi(theta);
        k_norm_pol(ii,2) = theta;
    end
    
    % 1.2. Build the Zernike matrix
    % 1.2.1. Assign the radial and angular orders for the Zernike
    % polynomials
    order = []; % Initialize a matrix that contain the radial and angular orders 
    for ii = 1:rad_order
        for jj = -ii:2:ii
            order = vertcat(order,[jj ii]);
        end
    end
    n_order = length(order); % The number of Zernike polynomials
    % 1.2.2. Build matrix
    Z_re = zeros(n1,n_order);
    for jj = 1:n1
        rho = k_norm_pol(jj,1);
        theta = k_norm_pol(jj,2);
        m = order(:,1);
        n = order(:,2);
        Z_re(jj,:) = zernfun(n,m,rho,theta,'norm');
    end
    Z_re(:,1:2) = 0; % Remove tip and tilt
    
    % 1.3. Save the matrix
    save(""+directory_save+"Z_"+division_step+"_re.mat",'Z_re')

    %% 2. Build the Zernike matrix when running optimization
    if division_step == 0 % If we are correcting the full image, then Z is Z_re, k is the same
        Z = Z_re;
        save(""+directory_save+"Z_"+division_step+".mat",'Z')
        save(""+directory_save+"k_"+division_step+".mat",'k')
    else
        % 2.1. Build new k space
        if division_step > 1
            dk_zone = pi/l_zone; % Spacing in k space
        else 
            dk_zone = 2*pi/l_zone;
        end
        kx_zone = -kt_max+dk_zone:dk_zone:kt_max; % New list of kx and ky
        ky_zone = kx_zone;
        N = length(kx_zone); % Number of new kx or ky
        k_zone = zeros(N^2,2); % Initialize the list of new kx and ky. 1st column is kx, 2nd column is ky
        for ii = 1:N
            k_zone((ii-1)*N+1:ii*N,1) = kx_zone(ii)*ones(N,1);
            k_zone((ii-1)*N+1:ii*N,2) = ky_zone;
        end
        % Filter the out-of-NA values
        k_zone = [k_zone, sqrt(k_zone(:,1).^2+k_zone(:,2).^2)];
        k_zone(k_zone(:,3) > kt_max,:) = [];
        k_zone = k_zone(:,1:2);
        % Save the new k space
        k = k_zone;
        save(""+directory_save+"k_"+division_step+".mat",'k')
        
        % 2.2. Now, build the Zernike matrix in the same way as part 1 of
        % this function
        
        % 2.2.1. Build the normalized polar k space
        % The maximum k_parallel
        n2 = length(k_zone);
        % Convert k to normalized polar coordinate
        k_norm_pol_zone = zeros(n2,2); % 1st column: radial coordinate, 2nd column: angular coordinate
        for ii = 1:n2
            % radial coordinate
            k_norm_pol_zone(ii,1) = sqrt(k_zone(ii,1).^2+k_zone(ii,2).^2)/kt_max;
            sine_theta = k_zone(ii,2)/sqrt(k_zone(ii,1).^2+k_zone(ii,2).^2); % Sine of the angular coordinate
            % Now, find the angular coordinate
            if k_zone(ii,1) >= 0 && (k_zone(ii,1) ~= 0 || k_zone(ii,2) ~= 0)
                theta = asin(sine_theta);
            elseif k_zone(ii,1) < 0 && (k_zone(ii,1) ~= 0 || k_zone(ii,2) ~= 0)
                theta = pi-asin(sine_theta);
            elseif k_zone(ii,1)==0 && k_zone(ii,2) == 0
                theta = 0;
            end
            % after this step, theta is from -pi/2 to 3*pi/2. wrap it to pi
            theta = wrapToPi(theta);
            k_norm_pol_zone(ii,2) = theta;
        end

        % 1.2. Build the Zernike matrix
        Z = zeros(n2,n_order);
        for jj = 1:n2
            rho = k_norm_pol_zone(jj,1);
            theta = k_norm_pol_zone(jj,2);
            m = order(:,1);
            n = order(:,2);
            Z(jj,:) = zernfun(n,m,rho,theta,'norm');
        end
        Z(:,1:2) = 0; % Remove tip and tilt

        % 1.3. Save the matrix
        save(""+directory_save+"Z_"+division_step+".mat",'Z')
    end
    clear
end