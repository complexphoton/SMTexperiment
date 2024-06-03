function [M,grad] = sharpness_figure_of_merit(a, F, Z)
% The function to minimize M = -sum(I*log(I)), where I is the intensity image. For maximization of -M, one can just change the sign of the grad variable.
% a is the vector containing the zernike coefficients, Z is the zernike
% matrix
    % make sure a is a column vector
    a = a(:);
    % The aberration angle is
    phi = Z*a;
    % Calculate the FOM
    psi = F*exp(1i*phi);
    I = abs(psi).^2;
    
    global M 
    
    M = single(-sum(I.*log(I),'all'));
    
    a = double(a);
    
% Calculate the gradient
    if nargout>1
        grad = double(2*imag(((log(I)+1).*(psi))'*F.*(exp(1i*phi.')))*Z); 
    end
end
