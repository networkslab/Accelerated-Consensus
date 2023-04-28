function T = get_sund_haj_time(E)
% Estimate the time that the algorithm proposed by Sundaram and Hajicostis 
% requires to achieve perfect convergence. The algorithm is described in
% S. Sundaram and C. Hadjicostis, “Distributed consensus and linear function 
% calculation in networks: An observability perspective,” in Proc. IEEE/ACM 
% Int. Symp. Information Proc. in Sensor Networks, Cambridge, MA, USA, Apr. 2007.
%% Inputs
% E the adjacency matrix
%% Outputs
% T the time required by the algorithm proposed by Sundaram and Hajicostis

N = size(E, 1);

d = sum(E, 2) - 1;

T = N - min(d);