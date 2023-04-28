function q = get_hermite(p, a)
% Create the sub-optimal weights for the polynomial filter based on the 
% Hermite interpolating polynomial. The algorithm and this sub-optimal
% initialization are discussed in 
% E. Kokiopoulou and P. Frossard, “Polynomial filtering for fast convergence 
% in distributed consensus,” IEEE Trans. Signal Process., vol. 57, no. 1, 
% pp. 342–354, Jan. 2009.
% Interfaces the difftable building code by Arpad Toth, arpi@elte.hu, available
% online at
% http://www.mathworks.com/matlabcentral/fileexchange/14353-hermite-interpolation
%% Inputs
% a the left border
% p the order of the polynomial filter
%% Outputs
% q the vector of the sub-optimal weights of the polynomial filter

A = zeros(2, p);
A(1, 1:2) = 1;
A(1, 3:p) = inf;
A(2, 1) = a;
A(2, 2) = 0;
 
q=difftable(A)';