function x = naive_gauss(A, b)
    [n, m] = size(A);
    if n ~= m
        error('계수 행렬 A는 정방행렬이어야 합니다.');
    end

    % Forward Elimination
    for k = 1:n-1          
        for i = k+1:n      
            
            if A(k,k) == 0
                error('피벗이 0입니다. 부분 피벗팅(Partial Pivoting)이 필요합니다.');
            end
            
            factor = A(i, k) / A(k, k);
            
            A(i, k:n) = A(i, k:n) - factor * A(k, k:n);
            
            b(i) = b(i) - factor * b(k);
        end
    end

    % Back Substitution
    x = zeros(n, 1);           
    x(n) = b(n) / A(n, n);      

    for i = n-1:-1:1             
        sum_val = 0;
        for j = i+1:n            
            sum_val = sum_val + A(i, j) * x(j);
        end
       
        x(i) = (b(i) - sum_val) / A(i, i);
    end
end