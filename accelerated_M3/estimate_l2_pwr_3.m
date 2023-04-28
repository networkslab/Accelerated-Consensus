function l2_lanc = estimate_l2_pwr_3(Wmh, Nit, NitM)
% Estimate the second (by the modulus) largest eigenvalue of a given 
% stochastic matrix based on the baseline Distributed Orthogonal Iteration
% algortithm. See S. Boyd, A. Ghosh, B. Prabhakar, and D. Shah, 
% “Randomized gossip algorithms.” IEEE Trans. Inf. Theory, vol. 52, no. 6,
% pp. 2508–2530, Jun. 2006.
%% Inputs
% Wmh - stochastic matrix
% Nit - number of DOI iterations
% NitM - number of iterations to calculate the orthogonal projection (mean)
%% Outputs
% l2_lanc - the estimate of the second (by the modulus) largest eigenvalue
% of Wmh

E = double(abs(Wmh) > 0);
N = length(Wmh);
theta = [-1/3; 0; 4/3];

% initialize
lam_init = 1 - 1/sqrt(N);
x = 100*randn(N, 1);

v = x;
for i = 1:Nit
    u1 = v;
    v = Wmh * v;
    u2 = v;
    
    [dummy1, dummy2, dummy3, x_mean] = ...
        do_consensus_acc_circ(NitM, Wmh, E, lam_init, theta, v, 'NONE', 0);
    
    v = v - x_mean(1:N);
    v = v ./ max(abs(v));
end;

n = u2.^2;
d = u1.^2;

% Do second consensus to determine \sum_i (x_i(t)-x_i(t+K))^2
[dummy1, dummy2, dummy3, n] = ...
    do_consensus_acc_circ(NitM, Wmh, E, lam_init, theta, n, 'NONE', 0);
% Do third consensus to determine \sum_i (x_i(t-1)-x_i(t+K-1))^2
[dummy1, dummy2, dummy3, d] = ...
    do_consensus_acc_circ(NitM, Wmh, E, lam_init, theta, d, 'NONE', 0);
% Compute eigenvalues at individual nodes
l2_lanc = sqrt(n(1:N)./ d(1:N));
% Use max consensus to converge to a single value at all nodes
l2_lanc = min(l2_lanc);

l2_lanc = min(l2_lanc, 0.99999999);
l2_lanc = max(l2_lanc, 0.00000001);

return;