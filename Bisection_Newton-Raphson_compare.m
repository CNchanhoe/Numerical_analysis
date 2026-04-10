clear; clc; close all;

% Function_definition 
f = @(x) exp(-2*x) - 7*x;         % f(x)
df = @(x) -2*exp(-2*x) - 7;       % f'(x)

tol = 1e-6; % Tolerance

% Bisection_method 
fprintf('Bisection Method\n');
fprintf('Iter |        a         |        b         |        c         |       f(c)       |      Error      \n');
fprintf('---------------------------------------------------------------------------------------------------\n');

% Initial_section_setting
a = 0; 
b = 1; 
iter_B = 0;
err_B = (b - a) / 2; % Initial_error
history_B = [];      % Error_history

while err_B > tol
    iter_B = iter_B + 1;
    c = (a + b) / 2; % Median_calculation
    
    % Error_calculation
    err_B = (b - a) / 2; 
    history_B(iter_B) = err_B;
    
    % Print_current_status
    fprintf('%4d | %16.10f | %16.10f | %16.10f | %16.10f | %16.10e \n', iter_B, a, b, c, f(c), err_B);
    
    % Termination_condition
    if f(c) == 0 || err_B < tol
        break;
    end
    
    % Reduce_section
    if f(a) * f(c) < 0
        b = c;
    else
        a = c;
    end
end
root_B = c;
fprintf('-> Bisection Method Root: %f (Iterations: %d)\n\n', root_B, iter_B);

% Newton-Raphson Method
fprintf('Newton-Raphson Method\n');
fprintf('Iter |        x_i       |      f(x_i)      |      Error      \n');
fprintf('-------------------------------------------------------------\n');

x0 = 0.5; % Initial_estimate_setting
iter_N = 0;
err_N = 1;    % Initial_error(bigger_than_tolerance)
history_N = []; % Error_history

while err_N > tol
    iter_N = iter_N + 1;
    
    % Recurrence_formula
    x1 = x0 - f(x0) / df(x0); 
    
    % Error_calculation(absolute_value)
    err_N = abs(x1 - x0);
    history_N(iter_N) = err_N;
    
    % Print_current_status
    fprintf('%4d | %16.10f | %16.10f | %16.10e \n', iter_N, x1, f(x1), err_N);
    
    % Update_variable
    x0 = x1;
end
root_N = x1;
fprintf('-> Newton-Raphson Root: %f (Iterations: %d)\n\n', root_N, iter_N);

% Figure Plotting 
figure('Name', 'Convergence Speed Comparison', 'Position', [100, 100, 700, 500]);
semilogy(1:iter_B, history_B, '-o', 'LineWidth', 2, 'DisplayName', 'Bisection Method');
hold on;
semilogy(1:iter_N, history_N, '-s', 'LineWidth', 2, 'DisplayName', 'Newton-Raphson Method');
grid on;
title('Convergence Speed: Bisection vs Newton-Raphson');
xlabel('Number of Iterations');
ylabel('Error (Log Scale)');
legend('Location', 'northeast');
set(gca, 'FontSize', 12);
