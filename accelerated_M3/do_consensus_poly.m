function [MSE, x_it] = do_consensus_poly(M, W, theta, x)

% Perform accelerated consensus iterations for the polynomial filter, 
% output the MSE curve and the history of consensus iterations
%% Inputs
% x - vector of initial sensor values
% theta - vestor of polynomial filter coefficients
% W - foundational consensus weight matrix
% M - number of polynomial filter consensus operations to be performed
%% Outputs
% MSE - Nit x 1 vector of mean squared error values at every iteration
% x_it - history of consensus iterations

p = length(theta);
N = length(x);
x_it = zeros(N, M);
x_it(:, 1) = x;
x_mean = mean(x);
MSE = zeros(M, 1);
one_vec = ones(1,N);
MSE(1) = one_vec * (x_it(:, 1) - x_mean).^2;

for i = 2:M
    
    if (mod(i, round(p)+1) == 0) 
        x_it(:, i) = x_it(:, i-p:i-1) * theta; 
        x_it(:, i) = W * x_it(:, i);
        MSE(i) = one_vec * (x_it(:, i) - x_mean).^2 ;
    else
        x_it(:, i) = W * x_it(:, i-1);
        MSE(i) = one_vec * (x_it(:, i) - x_mean).^2 ;
    end;
    
end;

return;