function [Wmax, E] = create_E_chain(N, tau)
% Create the connectivity matrix for the path of N nodes
% For maximum degree matrix See L. Xiao, S. Boyd, and S. Lall, “A scheme for robust distributed sensor fusion based on average consensus,” in Proc.
% IEEE/ACM Int. Symp. on Information Processing in Sensor Networks, Los Angeles, CA, Apr. 2005.
%% Inputs
% N number of nodes
% tau connectivity radius (not used)
%% Outputs
% E adjacency matrix
% Wmax maximum degree matrix

E = zeros(N);
one_vec = ones(2,1);
E(1,1) = 1;
E(1,2) = 1;
for i=2:N-1
    E(i,i) = 1;
    E(i,i-1) = 1;
    E(i,i+1) = 1;
end
E(N,N) = 1;
E(N,N-1) = 1;

Wmax = E / N;
for i = 1:N
    Wmax(i,i) = 1 - sum(Wmax(i, :)) + 1/N;
end;

return;