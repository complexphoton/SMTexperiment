function dispersion_VRM(directory_r,prefix,list_k0,im_case,directory_save)
    
    % 1. Setting
    % The number of frequency
    n_freq = length(list_k0);
    % The central frequency
    i_freq_center = round(n_freq/2);
    r_center = load(""+directory_r+""+prefix+""+i_freq_center+".mat").r_pad;
    
    % 2. Build matrices
    for i_freq = 1:n_freq
        fprintf("Working with frequency"+i_freq+"\n")
        
        % Load the reflection matrix at this frequency
        r = load(""+directory_save+"r_"+im_case+"_"+i_freq+"_new_ref_trunc.mat").r_trunc;
 
        % Do the iterative dispersion correction
        mean_phi = pi; N = width(r);
        phi_out = zeros(N,1); phi_in = zeros(N,1);
        while mean_phi > pi/20
            
            % Correct output dispersion
            phi_out_list = angle(diag(r*r_center'));
            %
            r = diag(exp(-1i*phi_out_list))*r;
            
            % Correct input dispersion
            phi_in_list = angle(diag(transpose(r)*conj(r_center)));
            
            r = r*diag(exp(-1i*phi_in_list));
            
            phi_list = vertcat(phi_in_list(phi_in_list ~= 0), phi_out_list(phi_out_list ~= 0));
            mean_phi = mean(abs(phi_list),'all');
            
            phi_out = phi_out+phi_out_list; phi_in = phi_in+phi_in_list;
        end
        
        % Save the corrected matrix
        save(""+directory_save+"r_"+im_case+"_"+i_freq+"_vrm.mat",'r','phi_in','phi_out')
    end
    
end
