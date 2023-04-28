% script for generating the path graph adjacency matrix, and various
% weights and weight matrices. Weight matrices generated: 
% Metropolis-Hastings, FDLA optimal. Weights generated: Polynomial weights
% for 3, 5, 7 order polynomial filter

clc;

l2_cutoff = 0.9999999999999;

% Network size
Nvec = [25; 50; 100; 150; 200; 250; 300; 400; 500; 600];
Nvec = 25;
% Number of trials for every network size
Mvec = 1;
% predictor parameters
theta = [-1/3; 0; 4/3];
% temporary eigenvalues

Mac = zeros(length(Nvec), 1);
% order of polynomial filtering
PH = 20;
PO20 = 6;
PO5 = 4;
PO3 = 2;

l2_mh_tens = zeros(max(Mvec), length(Nvec));
l2_opt_tens = zeros(max(Mvec), length(Nvec));

Wmh_tens = single(zeros(Nvec, Nvec, max(Mvec)));
Wopt_tens = zeros(Nvec, Nvec, max(Mvec));
E_tens = uint8(zeros(Nvec, Nvec, max(Mvec)));
g_tens = single(zeros(Nvec, 2, max(Mvec)));

alpha_polyO3_tens = zeros(PO3+1, max(Mvec));
alpha_polyO5_tens = zeros(PO5+1, max(Mvec));
alpha_polyO20_tens = zeros(PO20+1, max(Mvec));

opts.disp=0;
for j = 1:length(Nvec)
    j
    N = Nvec(j);
    tau = sqrt(2*log(N)/N);
    M = Mvec(j);
    for i = 1:M
        i
        % Create adjacency matrix and MD matrix 
        [Wmax, E] = create_E_chain(N, tau);
        % Compute MH matrix 
        Wmh = create_MH(E);
        m = (sum(sum(E)) - N)/2;
        A = zeros(N, m);
        if (N < 150)
            [V, Dmh] = eig(Wmh);
            temp_eig(2) = Dmh(N-1,N-1);
        else
            J = ones(N)/N;
            [V,Dmh,flag] = eigs(Wmh-J, 1,'LM', opts);
            temp_eig(2) = Dmh;
        end;
        
        if temp_eig(2) > l2_cutoff
            continue;
        end;
        Mac(j) = Mac(j) + 1;
        % Compute optimal matrix and store eigenvalues
        Wopt = Wmh;
        [Wopt, l2_opt(Mac(j), j)] = create_OPT(E);
        l2_mh(Mac(j), j) = temp_eig(2);
        
        [alpha_polyO3] = create_POLY(Wmh, PO3);
        [alpha_polyO5] = create_POLY(Wmh, PO5);
        [alpha_polyO20] = create_POLY(Wmh, PO20);
        
        
        l2_mh_tens(Mac(j)) = l2_mh(Mac(j), j);
        l2_opt_tens(Mac(j)) = l2_opt(Mac(j), j);

        Wmh_tens(:,:,Mac(j)) = Wmh;
        Wopt_tens(:,:,Mac(j)) = Wopt;
        E_tens(:,:,Mac(j)) = uint8(E);

        alpha_polyO3_tens(:,Mac(j)) = alpha_polyO3;
        alpha_polyO5_tens(:,Mac(j)) = alpha_polyO5;
        alpha_polyO20_tens(:,Mac(j)) = alpha_polyO20;
        
    end;
    
end;

filename = strcat('matrix_CIRC_N', num2str(Nvec));
save(filename, 'Mvec', 'l2_mh_tens', 'l2_opt_tens', ...
               'Wmh_tens', 'Wopt_tens', 'E_tens', ...
               'alpha_polyO3_tens', 'alpha_polyO5_tens', 'alpha_polyO20_tens');


