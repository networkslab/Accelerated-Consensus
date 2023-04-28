function [MSE, x_it] = do_consensus_poly_circ(M, W, theta, x)

% For the path graph, perform accelerated consensus iterations for the polynomial filter, 
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
x_it = zeros(N, p+1);
x_it(:, p) = x;
x_mean = mean(x);
MSE = zeros(M, 1);
one_vec = ones(1,N);  
MSE(1) = one_vec * (x - x_mean).^2;

for i = 2:M
    
    if (mod(i, round(p)+1) == 0) 
        x_it(:, p+1) = x_it(:, 1:p) * theta;  
        x_it(:, p+1) = W * x_it(:, p+1);
        MSE(i) = one_vec * (x_it(:, p+1) - x_mean).^2 ;
    else
        x_it(:, p+1) = W * x_it(:, p);
        MSE(i) = one_vec * (x_it(:, p+1) - x_mean).^2 ;
    end;
    x_it = circshift(x_it, [0, -1]);
end;

return;