function [l2_est, l2_M3, MSE, x_it] = ...
    do_consensus_acc_circ(M, W, E, l2, theta, x, EstMethod, Nit, NitM)

% For the path graph, estimate the eigenvalue of the foundational weight 
% matrix (if instructed) calculate the optimal value of mixing parameter,
% perform accelerated consensus iterations, output the MSE curve 
% and the history of consensus iterations
%% Inputs
% NitM - number of iterations to calculate orthogonal iteration in DOI (outdated and no longer used)
% Nit - number of iterations in DOI to estimate the second largest eigenvalue
% EstMethod - eigenvalue estimation method. String, the following options supported:
%       'PWR' - old version of the power iteration method with l_2 norm
%       'PWRN' - new version of the power iteration method with l_\infty
%       norm and mean subtraction
%       'NONE' - no eigenvalue estimation, eigenvalue is known
% x - vector of initial sensor values
% theta - vestor of predictor coefficients
% l2 - second largest eigenvalue of the foundational weight matrix
% E - adjacency matrix
% W - foundational consensus weight matrix
% M - number of accelerated consensus operations to be performed
%% Outputs
% l2_est - the estimate of the second largest eigenvalue of the foundational weight matrix
% l2_M3 - the (absolute) value of the second largest eigenvalue of the
% accelerated operator
% MSE - Nit x 1 vector of mean squared error values at every iteration
% x_it - history of consensus iterations

N = length(x);
x_it = zeros(2*N, 1);
x_it(1:N, 1) = x;
x_mean = mean(x);
MSE = zeros(M, 1);
one_vec = ones(1,N); %  / N
MSE(1) = one_vec * (x_it(1:N, 1) - x_mean).^2;

IM3 = eye(N);
OM3 = zeros(N);

% Estimate \lambda_2
if (strcmp(EstMethod,'PWR'))
    % distributed eigenvalue estimation based on power iterations
    l2_est = estimate_l2_pwr_1(W, Nit, NitM);
elseif (strcmp(EstMethod,'PWRN'))
    % distributed eigenvalue estimation based on power iterations
    l2_est = estimate_l2_new(W, E, N, Nit);
%     l2_est1 = estimate_l2_pwr_4(W, Nit, NitM);
elseif (strcmp(EstMethod,'NONE'))
    l2_est = l2;
end;
[alp, l2_M3] = get_alpha(l2, l2_est, theta);
WM3 = [(1+alp*(theta(3)-1))*W + alp*theta(2)*IM3      alp*theta(1)*IM3
          IM3                                    OM3];
      

for i = 2:M
    if (i > 2)
        x_it = WM3 * x_it;
        MSE(i) = one_vec * (x_it(1:N) - x_mean).^2 ;
    else
        x_it(1:N) = W * x_it(1:N);
        x_it(N+1:2*N) = x_it(1:N);
        MSE(i) = one_vec * (x_it(1:N) - x_mean).^2 ;
    end;
end;

return;