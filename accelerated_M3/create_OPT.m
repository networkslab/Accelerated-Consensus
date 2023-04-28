function [Wopt, l2_opt] = create_OPT(E)
% Create the optimal symmetric weight matrix based on the adjacency matrix. The
% optimal weight matrix is the solution of the symmetric FDLA problem from
% L. Xiao and S. Boyd, “Fast linear iterations for distributed averaging.” 
% Sys. and Control Letters, vol. 53, no. 1, pp. 65–78, Sep. 2004.
% Uses the cvx package, available online under GPL license,
% http://www.stanford.edu/~boyd/cvx/
% See M. Grant and S. Boyd. Graph implementations for nonsmooth convex 
% programs, Recent Advances in Learning and Control (a tribute to M. Vidyasagar), 
% V. Blondel, S. Boyd, and H. Kimura, editors, pages 95-110, Lecture Notes 
% in Control and Information Sciences, Springer, 2008
%% Inputs
% E adjacency matrix
%% Outputs
% Wopt the optimal symmetric weight matrix, the solution of the symmetric
% FDLA problem
% l2_opt the value of the second (by the modulus) largest eigenvalue of Wopt

N = size(E, 1);
m = (sum(sum(E)) - N)/2;
A = zeros(N, m);

k = 1;
for i = 1:N
    for j = i+1:N
        if (E(i,j)==1)
            A(i,k) = 1;
            A(j,k) = -1;
            k = k+1;
        end;
    end;
end;

[n,m] = size(A);
I = eye(n,n);
J = I - (1/n) * ones(n,n);

cvx_quiet(1)
cvx_begin sdp
    variable w(m,1)   % edge weights
    variable s        % epigraph variable
    variable L(n,n) symmetric
    minimize( s )
    subject to
        L == A * diag(w) * A';
        J - L <= +s * I;
        J - L >= -s * I;
cvx_end

k = 1;
Wopt = zeros(N);
for i = 1:N
    for j = i+1:N
        if (E(i,j)==1)
            Wopt(i,j) = w(k);
            Wopt(j,i) = w(k);
            k = k+1;
        end;
    end;
end;
for ii = 1:N
    Wopt(ii,ii) = 1-sum(Wopt(ii,:));
end;

l2_opt = s;