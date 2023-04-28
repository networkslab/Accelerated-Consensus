function l2_lanc = estimate_l2_pwr_4(Wmh, Nit, NitM)

% Estimate the second (by the modulus) largest eigenvalue of a given 
% stochastic matrix based on the modified and streamlined 
% Distributed Orthogonal Iteration algortithm (Algorithm 1 in the paper). 
%% Inputs
% Wmh - stochastic matrix
% Nit - number of DOI iterations
% NitM - number of iterations to calculate the orthogonal projection
% (mean). obsolete parameter
%% Outputs
% l2_lanc - the estimate of the second (by the modulus) largest eigenvalue
% of Wmh

N = length(Wmh);
x = randn(N, 1);
x = x - Wmh * x;

v = x;
for i = 1:Nit
    v = Wmh * v;
    v = v ./ max(abs(v));
end;

v = v - Wmh * v;
v = Wmh * v;
v = v ./ max(abs(v));
v = Wmh * v;

l2_lanc = max(abs(v));

l2_lanc = min(l2_lanc, 0.9999999999);
l2_lanc = max(l2_lanc, 0.0000000001);

return;