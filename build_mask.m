function [mask] = build_mask(x_im,y_im,wx0,wy0,x0,y0,rho)
    [gridX, gridY] = meshgrid(x_im,y_im);
    fgauss = exp(-1/(2*(1-rho^2))*(...
            ((gridX-x0)/wx0).^2+((gridY-y0)/wy0).^2)-2*rho.*(gridX-x0)/wx0.*(gridY-y0)/wy0...
            );
    Igauss = fgauss.^2;
    mask = 1./Igauss;
end