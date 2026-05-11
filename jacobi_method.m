function [x, iter] = jacobi_method(A, b, x0, tol, max_iter)
    n = length(b);
    x = x0;
    x_new = zeros(n, 1);
    
    for iter = 1:max_iter
        for i = 1:n
            sum_val = b(i);
            for j = 1:n
                if i ~= j
                    sum_val = sum_val - A(i, j) * x(j);
                end
            end
            x_new(i) = sum_val / A(i, i);
        end
        
        if norm(x_new - x, inf) < tol
            x = x_new;
            break;
        end
        
        x = x_new;
    end
end