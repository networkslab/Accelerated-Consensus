function lambdaest = estimate_l2_new(W,E,N,L)

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

% Mark's code
T=10;

X = ceil(4*rand(N,1))./4;

X2 = X-W*X;

for ii = 1:L

    X2 = W*X2;
    if rem(ii,T)==1
     X2 = X2./max(X2);
    end;
end;

% X2 = X2 - W * X2;
% X2 = W * X2;
% X2 = X2 ./ max(abs(X2));
% X2 = W * X2;

lambdaest = max(W*X2)/max(X2);