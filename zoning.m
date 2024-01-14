% Function to divide the image into equally-sized overlapped zones. The
% number of zones and their sizes depend on the division step.
function [x_step,y_step,x_im_step,y_im_step] = zoning(dx,dx_im,list_x,list_y,list_x_im,list_y_im,division_step,overlap2)
    %% The inputs:
    % list_x,list_y is the coordinate of the entire image when running optimization
    % list_x_im, list_y_im is the coordinate of the entire image when reconstructing
    % division_step determines the number of zones
    % overlap2 is half the size of the overlapping region, in micron

    % The outputs:
    % x_step, y_step, x_im_step, y_im_step are all column cells that store the
    % coordinates of each zones

    % The number of zones in x and y
    n_zone_x = 2^division_step;
    n_zone_y = n_zone_x;
    n_zone = n_zone_x*n_zone_y;
    
    % Size of the overlapping region in micron
    d = 2*overlap2;
    
    % Initialize the cell that contains the zones' coordinates
    x_step = cell(n_zone,1);
    y_step = cell(n_zone,1);
    x_im_step = cell(n_zone,1);
    y_im_step = cell(n_zone,1);
    
    %% Assign the coordinates for each zone
    for zone_id = 1:n_zone
        % Find the bigger zone that contains this zone
        zone_big = ceil(zone_id/4);
        % Then, find the coordinate of that big zone
        x_big = list_x{zone_big,division_step};
        y_big = list_y{zone_big,division_step};
        x_im_big = list_x_im{zone_big,division_step};
        y_im_big = list_y_im{zone_big,division_step};
        % Size of the big zone in micron
        l_big = length(x_big)*dx;
        % So the size of the zone is
        l = (l_big+1*d)/2;
        % Given the size of the zone and the overlapping region, along with
        % the coordinates of the bigger zone, the coordinates of each zone
        % can be assigned as
        switch mod(zone_id,4)
            case 1 % If mod(zone_id,4) == 1, the zone is on the upper right of the big zone
                x_step{zone_id,1} = fliplr(x_big(end):-dx:x_big(end)-l);
                y_step{zone_id,1} = y_big(1):dx:y_big(1)+l;
                x_im_step{zone_id,1} = fliplr(x_im_big(end):-dx_im:x_im_big(end)-l);
                y_im_step{zone_id,1} = y_im_big(1):dx_im:y_im_big(1)+l;
            case 2 % If 2, the zone is on the lower right
                x_step{zone_id,1} = fliplr(x_big(end):-dx:x_big(end)-l);
                y_step{zone_id,1} = fliplr(y_big(end):-dx:y_big(end)-l);
                x_im_step{zone_id,1} = fliplr(x_im_big(end):-dx_im:x_im_big(end)-l);
                y_im_step{zone_id,1} = fliplr(y_im_big(end):-dx_im:y_im_big(end)-l);
            case 3 % If 3, the zone is on the upper left
                x_step{zone_id,1} = x_big(1):dx:x_big(1)+l;
                y_step{zone_id,1} = y_big(1):dx:y_big(1)+l;
                x_im_step{zone_id,1} = x_im_big(1):dx_im:x_im_big(1)+l;
                y_im_step{zone_id,1} = y_im_big(1):dx_im:y_im_big(1)+l;
            case 0 % If 0, the zone is on the lower left
                x_step{zone_id,1} = x_big(1):dx:x_big(1)+l;
                y_step{zone_id,1} = fliplr(y_big(end):-dx:y_big(end)-l);
                x_im_step{zone_id,1} = x_im_big(1):dx_im:x_im_big(1)+l;
                y_im_step{zone_id,1} = fliplr(y_im_big(end):-dx_im:y_im_big(end)-l);
        end
    end
end