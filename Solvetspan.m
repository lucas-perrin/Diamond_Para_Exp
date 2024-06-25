function [X_sol_mat] = Solvetspan(X,tspan)
    sz = size(X(0));
    mx = sz(1); my = sz(2);
    nt = length(tspan);
    X_sol_mat = zeros(mx*my,nt);
    for i = 1:nt
        X_sol_mat(:,i) = reshape(X(tspan(i)),mx*my,1);
    end
    end