function [MSE, x] = do_consensus(Nit, W, x)
% Perform consensus iterations, output the MSE curve and the final value
%% Inputs
% x - vector of initial sensor values
% W - consensus weight matrix
% Nit - number of iterations to be performed
%% Outputs
% MSE - Nit x 1 vector of mean squared error values at every iteration
% x - vector of sensor values after Nit consensus iterations

N = length(x);
x_mean = mean(x);
MSE = zeros(Nit, 1);
one_vec = ones(1,N); %  / N
MSE(1) = one_vec * (x - x_mean).^2;
for i = 2:Nit
    x = W * x;
    MSE(i) = one_vec * (x - x_mean).^2 ;
end;

return;