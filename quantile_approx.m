function [invcdf_m, invcdf_c] = quantile_approx(invert, by, p_to, p_from, x_init, f, d_f, dd_f, ddd_f, max_overapprox, linearize_from)
    p = p_from:by:p_to;
    k = size(p, 2);
    x = zeros(k, 1);
    x(1) = x_init;
    
    beta_i = -log(p);
    
    d_x = @(x, beta) ...
        -exp(-beta) / f(x);
    dd_x = @(x, beta) ...
        (exp(-beta) - ...
         d_f(x)*(d_x(x, beta))^2 ) / f(x);
    ddd_x = @(x, beta) ...
        (-exp(-beta) - ...
         3*d_f(x)*d_x(x, beta)*dd_x(x, beta) - ...
         dd_f(x)*(d_x(x, beta))^3 ) / f(x);
    dddd_x = @(x, beta) ...
        (exp(-beta) - ...
         4*d_f(x)*d_x(x, beta)*ddd_x(x, beta) - ...
         3*d_f(x)*(dd_x(x, beta))^2 - ...
         6*dd_f(x)*(d_x(x, beta))^2* dd_x(x, beta) - ...
         ddd_f(x)*(d_x(x, beta))^4 ) / f(x);

    x_p_i = @(x, beta, p, p_0) ...
        x - ...
        d_x(x, beta) * log(p/p_0) + ...
        dd_x(x, beta) * (log(p/p_0))^2 / 2 - ...
        ddd_x(x, beta) * (log(p/p_0))^3 / 6 + ...
        dddd_x(x, beta) * (log(p/p_0))^4 / 24;

    for i = 1:(k-1)
        x(i+1) = x_p_i(x(i), beta_i(i), p(i+1), p(i));
    end
    
    [~, closestIndex] = min(abs(p - linearize_from));
    if linearize_from < p(closestIndex)
        closestIndex = closestIndex-1;
    end
       
    if invert
       x = flip(x);
       p = flip(1 - p);
    end
        
    invcdf_m = [];
    invcdf_c = [];
    
    
    if invert
        current_index = 1;
        index =  k - closestIndex;
    else
        current_index = closestIndex;
        index = k;
    end
    
    last_index = index;
   
    while current_index < last_index
        continue_trigger = 0;
        m = (x(index) - x(current_index)) / ((index - current_index)*by);
        c = x(current_index) - p(current_index) * m;
        
        if index == current_index+1
            invcdf_m = [invcdf_m; m];
            invcdf_c = [invcdf_c; c]; 
            if index ~= last_index
                current_index = index;
                index = last_index;
                continue
            else
                break
            end
        end
        
        for i = (current_index+1):(index-1)
            x_est = m * p(i) + c;
            x_error = abs(x(i) - x_est);
            if (x_error > max_overapprox)
                index = index - 1;
                continue_trigger = 1;
                break
            end
        end
        
        if continue_trigger == 1
            continue
        end
        
        invcdf_m = [invcdf_m; m];
        invcdf_c = [invcdf_c; c];
        current_index = index;
        index = last_index;
    end
end

