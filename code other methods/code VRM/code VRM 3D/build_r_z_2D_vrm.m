% The function that builds the deep-resolved reflection
% matrix for 2D imaging
function build_r_z_2D_vrm(list_amp,list_k0,directory_save,directory_r,k,id_sub)
    % 1. Setting
    % The number of frequency
    n_freq = length(list_k0);
   
    N = length(k); % Number of input/output channels
    
    % 2. Build matrices
    r_z = zeros(N,N);
    for i_freq = 1:n_freq
        fprintf("Building r_z frequency"+i_freq+"\n")
        
        % Load the dispersion corrected image
        r = load(""+directory_r+"r_vrm_"+i_freq+"_subvolume_"+id_sub+".mat").r;
        
        % Sum over frequency
        r_z = r_z+r/list_amp(i_freq);
    end
    
    % Save and clear
    save(""+directory_save+"./r_z_2D_vrm.mat",'r_z')
end
