function [x, iter] = gauss_seidel(A, b, x0, tol, max_iter)
    n = length(b);
    x = x0;
    
    for iter = 1:max_iter
        x_old = x;
        
        for i = 1:n
            sum_val = b(i);
            for j = 1:n
                if i ~= j
                    sum_val = sum_val - A(i, j) * x(j);
                end
            end
            x(i) = sum_val / A(i, i);
        end
        
        if norm(x - x_old, inf) < tol
            break;
        end
    end
end