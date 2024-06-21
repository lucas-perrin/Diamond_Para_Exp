function [X_para_end,time] = ParaExp_end(p,tspan,M,G,x0)

    times_T1 = zeros(1,p);
    times_T2 = zeros(1,p);
    time=0;

    tic;

    mx              = length(x0);
    nt              = length(tspan);
    size_int        = fix(nt/p);

    X_para      = zeros(mx,nt);
    X_para(:,1) = x0;

    ini_zero        = zeros(mx,1);

    x0p            = zeros(mx,p);
    x0p(:,1)       = x0;

    time=time+toc;

    for iter = 1:p
        
        tic;

    %----- Interval definition
        
        deb       = 1*(iter==1) + (iter-1)*size_int*(iter > 1);
        fin       = iter*size_int;
        if iter == p
            fin = length(tspan);
        end
        time_interval = tspan(deb:fin);
        
    %----- Type 1
                
        [~,X_T1_p] = odeRK4_inhom_ufunc(M,G,time_interval,ini_zero);
                
        X_para(:,deb:fin) = X_para(:,deb:fin) + X_T1_p;
        
        if iter < p
            x0p(:,iter+1) = X_T1_p(:,end);
        end
        
        times_T1(p) = toc;

    end

    for iter = 1:p
 
        tic;

    %----- Interval definition
        
        deb       = 1*(iter==1) + (iter-1)*size_int*(iter > 1);
        time_interval_till_end_from_zero = tspan(deb:end) - tspan(deb);
        
    %----- Type 2

        %[~,X_T2_p_end] = odeRK4_hom(M,time_interval_till_end_from_zero,x0p(:,iter));
        %[~,X_T2_p_end] = odeRK4_inhom_ufunc(M,@(t) zeros(mx,1),time_interval_till_end_from_zero,x0p(:,iter));
        [X_T2_p_end] = expm(time_interval_till_end_from_zero(end).*M)*x0p(:,iter);
                
        %X_para(:,end) = X_para(:,end) + X_T2_p_end;
        X_para(:,end) = X_para(:,end) + X_T2_p_end(:,end);

        times_T2(p) = toc;

    end

    X_para_end = X_para(:,end);

    time = time + max(times_T1) + max(times_T2);

end