clear; clc; close all;

% 1. definition
A = [7 -2 1 0;
     1 -9 3 -1;
     2 0 10 1;
     1 -1 1 6];
b = [17; 13; 15; 10];

n = length(b);

% 1. Naive Gauss Elimination Method
fprintf('--- 1. Naive Gauss Elimination Method ---\n');
Aug = [A b]; % Augmented matrix
x_gauss = zeros(n, 1);

% Forward Elimination
for k = 1:n-1
    for i = k+1:n
        factor = Aug(i,k) / Aug(k,k);
        Aug(i,k:n+1) = Aug(i,k:n+1) - factor * Aug(k,k:n+1);
    end
end

% Back Substitution
x_gauss(n) = Aug(n,n+1) / Aug(n,n);
for i = n-1:-1:1
    x_gauss(i) = (Aug(i,n+1) - Aug(i,i+1:n) * x_gauss(i+1:n)) / Aug(i,i);
end

fprintf('Solution (x):\n');
disp(x_gauss);
fprintf('\n');

% Iterative Methods Setup 
tol = 1e-6;          % Tolerance 
max_iter = 100;      % Maximum iterations
x0 = zeros(n, 1);    % Initial guess 

% 2. Jacobi Iteration Method
fprintf('--- 2. Jacobi Iteration Method ---\n');
fprintf(' Iter |     x1     |     x2     |     x3     |     x4     |    Error    \n');
fprintf('--------------------------------------------------------------------------\n');

x_jacobi = x0;
x_new_J = zeros(n, 1);
err_J = 100;
iter_J = 0;
history_err_J = []; % For plotting

while err_J > tol && iter_J < max_iter
    iter_J = iter_J + 1;
    
    for i = 1:n
        sum_val = 0;
        for j = 1:n
            if j ~= i
                sum_val = sum_val + A(i,j) * x_jacobi(j); % Uses ONLY old values
            end
        end
        x_new_J(i) = (b(i) - sum_val) / A(i,i);
    end
    
    % Infinity Norm Error 
    err_J = norm(x_new_J - x_jacobi, inf); 
    history_err_J(iter_J) = err_J;
    
    fprintf('%4d  | %10.6f | %10.6f | %10.6f | %10.6f | %11.6e\n', ...
            iter_J, x_new_J(1), x_new_J(2), x_new_J(3), x_new_J(4), err_J);
            
    x_jacobi = x_new_J; % Update state 
end
fprintf('=> Jacobi converged in %d iterations.\n\n', iter_J);

% 3. Gauss-Seidel Method
fprintf('--- 3. Gauss-Seidel Method ---\n');
fprintf(' Iter |     x1     |     x2     |     x3     |     x4     |    Error    \n');
fprintf('--------------------------------------------------------------------------\n');

x_gs = x0;
x_new_GS = x0;
err_GS = 100;
iter_GS = 0;
history_err_GS = []; % Plotting

while err_GS > tol && iter_GS < max_iter
    iter_GS = iter_GS + 1;
    
    for i = 1:n
        sum_val = 0;
        for j = 1:n
            if j ~= i
                % Uses newly calculated values
                sum_val = sum_val + A(i,j) * x_new_GS(j); 
            end
        end
        x_new_GS(i) = (b(i) - sum_val) / A(i,i);
    end
    
    % Infinity Norm Error
    err_GS = norm(x_new_GS - x_gs, inf);
    history_err_GS(iter_GS) = err_GS;
    
    fprintf('%4d  | %10.6f | %10.6f | %10.6f | %10.6f | %11.6e\n', ...
            iter_GS, x_new_GS(1), x_new_GS(2), x_new_GS(3), x_new_GS(4), err_GS);
            
    x_gs = x_new_GS; % Update state
end
fprintf('=> Gauss-Seidel converged in %d iterations.\n\n', iter_GS);

% 4. Figure Plotting
figure('Name', 'Convergence Speed Comparison', 'Position', [100, 100, 700, 500]);
semilogy(1:iter_J, history_err_J, '-o', 'LineWidth', 2, 'DisplayName', 'Jacobi Method');
hold on;
semilogy(1:iter_GS, history_err_GS, '-s', 'LineWidth', 2, 'DisplayName', 'Gauss-Seidel Method');
grid on;
title('Convergence Speed: Jacobi vs. Gauss-Seidel');
xlabel('Number of Iterations');
ylabel('Error (Infinity Norm) - Log Scale');
legend('Location', 'northeast');
set(gca, 'FontSize', 12);
