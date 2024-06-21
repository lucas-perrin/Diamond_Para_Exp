% fonction x'(t) = F(x,t) observateur
%
    function res = F(x,t)

    dx    = 0.05;
    xspan = 0:dx:1;
    mx    = length(xspan);

    %% ----- State trajectory : continuous in time solution

    alpha = 1;

    x_ref     = @(x,t) sin(pi.*x').*(x'<1) * sin(alpha.*t/2+2);
    x_ref_dt   = @(x,t) (alpha/2).*sin(pi.*x').*(x'<1) * cos(t./2+2);

    %% ----- Semi-continuous discretization

    A = ((mx-1)^2).*tridiag(1,-2,1,mx); A(1:1) = 1; A(end,end) = 1;
    A(1,2) = 0; A(2,1) = 0; A(end,end-1) = 0; A(end-1,end) = 0;

    B = eye(mx);

    X_ref    = @(t) x_ref(xspan,t);
    X_ref_dt  = @(t) x_ref_dt(xspan,t);
    X_ref_dxx = @(t) A*X_ref(t);

    res = @(t) X_ref_dt(t) - alpha.*X_ref_dxx(t);
end
%
% fonction x'(t) = L(x,t,y) observateur
% 
function res=L(x,t,y)

    m_obs = %%%% RENTRER M_OBS
    liste = m_obs.*(1:mx);
    liste = liste(m_obs>0);
    
    Id = eye(mx);
    C = Id(:,liste);
    gain  = 1;
    
    res = F(x,t) - gain*C*(h(x) - y);
end
%
% fonction h(x) = C*x

function res=h(x)
    m_obs = %%%% RENTRER M_OBS
    liste = m_obs.*(1:mx);
    liste = liste(m_obs>0);
    Id = eye(mx);
    C = Id(:,liste);
    res = C'*x;
end