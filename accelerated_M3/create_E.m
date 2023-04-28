function [Wmax, E, g] = create_E(N, tau)
% Create the connectivity matrix for the 2D random geometric graph
% For maximum degree matrix See L. Xiao, S. Boyd, and S. Lall, “A scheme for robust distributed sensor fusion based on average consensus,” in Proc.
% IEEE/ACM Int. Symp. on Information Processing in Sensor Networks, Los Angeles, CA, Apr. 2005.
%% Inputs
% N number of nodes
% tau connectivity radius
%% Outputs
% g matrix of mutual distances
% E adjacency matrix
% Wmax maximum degree matrix

%% 
u1 = rand(N,1);
u2 = rand(N,1);
g = [u1,u2];
E = zeros(N);
one_vec = ones(2,1);
for i=1:N
    temp1 = repmat(g(i,:), N, 1);
    dis = (temp1 - g).^2;
    dis = dis * one_vec;
    dis = sqrt(dis);
    E(i,:) = double(dis <= tau);
end
    

Wmax = E / N;
for i = 1:N
    Wmax(i,i) = 1 - sum(Wmax(i, :)) + 1/N;
end;

return;