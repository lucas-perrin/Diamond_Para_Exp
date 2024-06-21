function [nb_obs,C] = GetC(obs_map)

    C_flat = obs_map(:)';
    
    nb_obs = sum(C_flat ~= 0);
    
    C = zeros(nb_obs,length(C_flat));
    
    j = 1;
    
    for p = 1:length(C_flat)
        if C_flat(p) ~= 0
            C(j,p) = C_flat(p);
            j = j+1;
        end
    end
    
end