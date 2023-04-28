function Wmh = create_MH(E)
% Create the Metropolis-Hastings weight matrix based on the adjacency matrix
% W_{ij} = 1 / (1 + max(d_i, d_j))
% See L. Xiao, S. Boyd, and S. Lall, “A scheme for robust distributed sensor fusion based on average consensus,” in Proc.
% IEEE/ACM Int. Symp. on Information Processing in Sensor Networks, Los Angeles, CA, Apr. 2005.
%% Inputs
% E adjacency matrix
%% Outputs
% Wmh the Metropolis-Hastings weight matrix

N = size(E, 1);
Wmh = zeros(N,N);
d = sum(E, 2) - 1;
for ii = 1:N
    for jj = 1:N
        if ii~=jj & (E(ii,jj)>0)
                Wmh(ii,jj) = 1/(1+max(d(ii),d(jj)));
        end;
    end;
end;
for ii = 1:N
    Wmh(ii,ii) = 1-sum(Wmh(ii,:));
end;