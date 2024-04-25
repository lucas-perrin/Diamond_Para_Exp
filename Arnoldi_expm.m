function [Q,h,exptAb,n_bk] = Arnoldi_expm(A,b,n,t,s,toler)
    % arnoldi_exp :
    % computes exp(t*A)*b where t is a scalar reprensenting time,
    % A is a square matrix m*m and b a vector of size m.
    % 
    % Input :
    % - matrix A;
    % - vector b;
    % - n : max size of the Krylov space;
    %
    % (optional input)
    % - time t : in order to compute exp(t*A)*b (default = 1);
    % - s : if entered, computes the rational Krylov approx with scaling s (default = +\infty);
    % - toler : stopping criteria (default = 1e-12);
    %
    % Output :
    % - matrix Q of size (n+1,n);
    % - matrix h of size (m,n+1);
    % - exptAb = exp(t*A)*b;
    % - n_bk : number of iteration it took before reaching tolereance and breaking loop
    %
    % This algorithm is much faster if you feed it with a sparse matrix.
    
    if nargin < 4
        t = 1;
    end
    
    if nargin < 5
        s_null = true;
        S = A;
    else
        s_null = false;
        S = (speye(length(A)) - (1/s).*A)\A;
    end
    
    if nargin < 6
        toler = 1e-12;
    end
    
    h = spalloc(n+1,n,(n+1)*(n+1)/2+n);
    Q = zeros(length(S),n+1);
    
    Q(:,1) = b./norm(b);
    
    bk = false;
    
    for k = 2:n+1
        v = S * Q(:,k-1);
        for j = 1:k
            h(j,k-1) = Q(:,j)' * v;
            v = v - h(j,k-1) * Q(:,j);
        end
        h(k,k-1) = norm(v);
        if h(k,k-1) > toler
            Q(:,k) = v / h(k,k-1);
        else
            bk = true;
            n_bk = k;
            break
        end
    end
    
    if ~ bk
        n_bk = n;
    end
    
    Q = Q(:,1:n_bk);
    h = h(1:n_bk,1:n_bk);
    
    if s_null
        exptAb = Q * expm(t.*h) * Q' * b;
    else
        exptAb = Q * expm(t.* inv(inv(h) + (1/s).*speye(n_bk))) * Q' * b;
    end
    
    end