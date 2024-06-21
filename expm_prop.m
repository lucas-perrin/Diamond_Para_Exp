function [exptAv] = expm_prop(A,tspan,v)
    %EXPM_PROP Time exponential of a matrix using the built-in matlab function.
    %
    %   given a matrix A, a time span tspan, and a initial condition v, we
    %   compute the timewise exponential propagation of the matrix A using the
    %   built-in function of matlab "expm()" based on
    %expA = expm(A);
    exptAv = zeros(length(v),length(tspan));
    for i = 1:length(tspan)
        %exptAv(:,i) = expA^(tspan(i))*v;
        exptAv(:,i) = expm(tspan(i)*A)*v;
    end
end