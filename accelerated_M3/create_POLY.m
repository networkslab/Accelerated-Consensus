function [alpha] = create_POLY(W, p)
% Create the optimal weights for the polynomial filter, which are the 
% solution of the spectral radius minimization problem. The algorithm and
% the problem definition are discussed in
% E. Kokiopoulou and P. Frossard, “Polynomial filtering for fast convergence 
% in distributed consensus,” IEEE Trans. Signal Process., vol. 57, no. 1, 
% pp. 342–354, Jan. 2009.
% Uses the cvx package, available online under GPL license,
% http://www.stanford.edu/~boyd/cvx/
% See M. Grant and S. Boyd. Graph implementations for nonsmooth convex 
% programs, Recent Advances in Learning and Control (a tribute to M. Vidyasagar), 
% V. Blondel, S. Boyd, and H. Kimura, editors, pages 95-110, Lecture Notes 
% in Control and Information Sciences, Springer, 2008
%% Inputs
% W the initial weight matrix
% p the order of the polynomial filter
%% Outputs
% alpha the vector of the optimal weights of the polynomial filter


n = size(W, 1);

I = eye(n,n);
J = (1/n) * ones(n,n);

T = zeros(n,n,p+1);
T(:,:,1) = W;
for k = 2:p+1
    T(:,:,k) = W * T(:,:,k-1);%  W * T(:,:,k-1);
end;

cvx_begin sdp

    cvx_quiet(1)

    variable alpha(p+1,1)   % edge weights
    variable s        % epigraph variable
    minimize( s )
    expression L(n,n)
    L = zeros(n,n);
    for k = 1:p+1
            L = L + alpha(k)*T(:,:,k);
    end;
    subject to
        J - L <= +s * I;
        J - L >= -s * I;
        sum(alpha) == 1;
cvx_end
