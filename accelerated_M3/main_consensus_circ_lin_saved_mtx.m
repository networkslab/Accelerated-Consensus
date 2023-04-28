% main script for the path graph. script for obtaining the MSE and averaging time
% characterization of the proposed algorithm. Run this script to obtain the
% simulation results used to create figures in the paper. The simulation
% results are saved in a separate file for a graph with the number of nodes
% N. The data saved by this script can be furtheer accessed by the plotting
% routines.


clc;
l2_cutoff = 0.9999999999999;

% Network size and number of iterations
Nvec=25; Mc=5000; 
% Nvec=50; Mc=10000; 
% Nvec=100; Mc=40000; 
% Nvec=150; Mc=100000; 
Nvec=200; Mc=200000;
% % Nvec = 200;
% % Mc = 5000; 
load(strcat('RGG_matrices/matrix_circ_N', num2str(Nvec),'.mat'));



% Number of trials for every network size
Mvec = 1*ones(length(Nvec),1); % 300
% predictor parameters
theta = [-1/3; 0; 4/3];

Mac = zeros(length(Nvec), 1);
% order of polynomial filtering
PH = 20;
PO20 = 6;
PO5 = 4;
PO3 = 2;
alpha_polyH = get_hermite(PH+1, 0.2);
alpha_polyH = flipud(alpha_polyH);

l2_mh = zeros(max(Mvec), length(Nvec));
l2_opt = zeros(max(Mvec), length(Nvec));

l2_mhM3 = zeros(max(Mvec), length(Nvec));
l2_optM3 = zeros(max(Mvec), length(Nvec));

l2_mh_est_max = zeros(max(Mvec), length(Nvec));
l2_mh_est_pwr = zeros(max(Mvec), length(Nvec));
l2_opt_est_max = zeros(max(Mvec), length(Nvec));
l2_opt_est_pwr = zeros(max(Mvec), length(Nvec));

l2_mh_est_maxM3 = zeros(max(Mvec), length(Nvec));
l2_mh_est_pwrM3 = zeros(max(Mvec), length(Nvec));
l2_opt_est_pwrM3 = zeros(max(Mvec), length(Nvec));
l2_opt_est_maxM3 = zeros(max(Mvec), length(Nvec));

MSE_mh = zeros(Mc, 1);
MSE_opt = zeros(Mc, 1);

MSE_mh_pwrM3 = zeros(Mc, 1);
MSE_mh_maxM3 = zeros(Mc, 1);
MSE_opt_pwrM3 = zeros(Mc, 1);
MSE_opt_maxM3 = zeros(Mc, 1);
MSE_mhM3 = zeros(Mc, 1);
MSE_optM3 = zeros(Mc, 1);

MSE_mh_polyH = zeros(Mc, 1);
MSE_mh_polyO3 = zeros(Mc, 1);
MSE_mh_polyO5 = zeros(Mc, 1);
MSE_mh_polyO20 = zeros(Mc, 1);

MSE_mh_it = zeros(Mc, 1);
MSE_opt_it = zeros(Mc, 1);
MSE_mh_pwrM3_it = single(zeros(Mc, max(Mvec)));
MSE_mh_maxM3_it = zeros(Mc, 1);
MSE_opt_pwrM3_it = single(zeros(Mc, max(Mvec)));
MSE_opt_maxM3_it = zeros(Mc, 1);
MSE_mhM3_it = zeros(Mc, 1);
MSE_optM3_it = zeros(Mc, 1);
MSE_mh_polyH_it = zeros(Mc, 1);
MSE_mh_polyO3_it = zeros(Mc, 1);
MSE_mh_polyO5_it = zeros(Mc, 1);
MSE_mh_polyO20_it = zeros(Mc, 1);

MSE_mh_pwrM3_SHT_it = zeros(max(Mvec), 1);
MSE_opt_pwrM3_SHT_it = zeros(max(Mvec), 1);

opts.disp=0;

mu = 1;
sig = 1;
for j = 1:length(Nvec)
    j
    N = Nvec(j);
    tau = sqrt(log(N)/N);
    M = Mvec(j);
    % Number of Consensus iterations for Lanczos algorithm
    Nit = floor((Nvec)^2); 
    NitM = floor(Nit);
    J = ones(N)/N;
    IM3 = eye(N);
    OM3 = zeros(N);
    
    % Compute MH matrix 
    E = single(E_tens(:,:,j));
    Wmh = Wmh_tens(:,:,j);
    Wopt = Wopt_tens(:,:,j);
    if l2_mh_tens(j) > l2_cutoff
        continue;
    end;
    Mac(j) = Mac(j) + 1;
    % Compute optimal matrix and store eigenvalues
    l2_opt(Mac(j), j) = l2_opt_tens(j);
    l2_mh(Mac(j), j) = l2_mh_tens(j);

    Tsund_haj = get_sund_haj_time(E);
    
    alpha_polyO3 = alpha_polyO3_tens(:, j);
    alpha_polyO5 = alpha_polyO5_tens(:, j);
    alpha_polyO20 = alpha_polyO20_tens(:, j);
    
    x = (1:N)';
    mean_x = mean(x);
    x = sqrt(sig) * (x-mean_x) / norm(x-mean_x) + mean_x;
    
    [MSE_mh_it(:,1), x_mh_it] = do_consensus(Mc, Wmh, x);
    [MSE_opt_it(:,1), x_opt_it] = do_consensus(Mc, Wopt, x);
    
    [dummy, l2_mhM3(1, j), MSE_mhM3_it(:,1), x_mhM3_it] = do_consensus_acc_circ(...
            Mc, Wmh, E, l2_mh(1, j), theta, x, 'NONE', Nit, NitM);
    [dummy, l2_mhM3(1, j), MSE_mh_maxM3_it(:,1), x_mhM3_it] = do_consensus_acc_circ(...
            Mc, Wmh, E, 1-1/sqrt(N), theta, x, 'NONE', Nit, NitM);
        
    [dummy, l2_optM3(1, j), MSE_optM3_it(:,1), x_opt_M3_it] = do_consensus_acc_circ(...
            Mc, Wopt, E, l2_opt(1, j), theta, x, 'NONE', Nit, NitM);
    [dummy, l2_optM3(1, j), MSE_opt_maxM3_it(:,1), x_opt_M3_it] = do_consensus_acc_circ(...
            Mc, Wopt, E, 1-1/sqrt(N), theta, x, 'NONE', Nit, NitM);
        
    [MSE_mh_polyH_it(:,1), x_mh_polyH_it] = ...
            do_consensus_poly_circ(Mc, Wmh, alpha_polyH, x);
    [MSE_mh_polyO3_it(:,1), x_mh_polyO3_it] = ...
        do_consensus_poly_circ(Mc, Wmh, alpha_polyO3, x);
    [MSE_mh_polyO5_it(:,1), x_mh_polyO5_it] = ...
        do_consensus_poly_circ(Mc, Wmh, alpha_polyO5, x);
    [MSE_mh_polyO20_it(:,1), x_mh_polyO20_it] = ...
        do_consensus_poly_circ(Mc, Wmh, alpha_polyO20, x);

    
    for i = 1:M
        i
        
        l2_mhM3(i, j) = l2_mhM3(1, j);
        l2_optM3(i, j) = l2_optM3(1, j);
        
        [l2_mh_est_pwr(i, j), l2_mh_est_pwrM3(i, j), MSE_mh_pwrM3_it(:,i), dummy] = do_consensus_acc_circ(...
            Mc, Wmh, E, l2_mh(i, j), theta, x, 'PWRN', Nit, NitM);
        
        [l2_opt_est_pwr(i, j), l2_opt_est_pwrM3(i, j), MSE_opt_pwrM3_it(:,i), dummy] = do_consensus_acc_circ(...
            Mc, Wopt, E, l2_opt(i, j), theta, x, 'PWRN', Nit, NitM);
        

        MSE_mh_pwrM3 = MSE_mh_pwrM3 + MSE_mh_pwrM3_it(:,i);
        MSE_opt_pwrM3 = MSE_opt_pwrM3 + MSE_opt_pwrM3_it(:,i);
    end;
    
end;
MSE_mh = MSE_mh_it;
MSE_opt = MSE_opt_it;

MSE_mh_pwrM3 = MSE_mh_pwrM3 ./ Mvec;
MSE_mh_maxM3 = MSE_mh_maxM3_it;
MSE_mhM3 = MSE_mhM3_it;

MSE_opt_pwrM3 = MSE_opt_pwrM3 ./ Mvec;
MSE_opt_maxM3 = MSE_opt_maxM3_it;
MSE_optM3 = MSE_optM3_it;

MSE_mh_polyH = MSE_mh_polyH_it;
MSE_mh_polyO3 = MSE_mh_polyO3_it;
MSE_mh_polyO5 = MSE_mh_polyO5_it;
MSE_mh_polyO20 = MSE_mh_polyO20_it;


FontSize = 14;
h=plot(1:Mc, 10*log10(MSE_mh), 1:Mc, 10*log10(MSE_opt), ...
     1:Mc, 10*log10(MSE_mh_pwrM3), 1:Mc, 10*log10(MSE_mh_maxM3), ...
     1:Mc, 10*log10(MSE_optM3), 1:Mc, 10*log10(MSE_mhM3), ...
     1:Mc, 10*log10(MSE_mh_polyO3), 1:Mc, 10*log10(MSE_mh_polyO5), ...
     1:Mc, 10*log10(MSE_mh_polyO20), 1:Mc, 10*log10(MSE_mh_polyH));
legend('MH', 'OPT', 'MH-pwrM3', 'MH-subM3', 'OPT-M3', 'MH-M3', ...
    'POLY-SDP3', 'POLY-SDP5', 'POLY-SDP20', 'POLY-HERM', 'location', 'best');
set(gca, 'FontSize', FontSize); 
set(h, 'linewidth', 2);
xlabel('{\it t}', 'FontSize', FontSize);
ylabel('MSE, dB', 'FontSize', FontSize);

epsilon = -100; % dB
Tave_mh_it = sum(double(10*log10(MSE_mh_it) >= epsilon), 1)' -1;
Tave_opt_it = sum(double(10*log10(MSE_opt_it) >= epsilon), 1)' -1;

Tave_mh_pwrM3_it = sum(double(10*log10(MSE_mh_pwrM3_it) >= epsilon), 1)' -1;
Tave_mh_maxM3_it = sum(double(10*log10(MSE_mh_maxM3_it) >= epsilon), 1)' -1;
Tave_mhM3_it = sum(double(10*log10(MSE_mhM3_it) >= epsilon), 1)' -1;

Tave_opt_pwrM3_it = sum(double(10*log10(MSE_opt_pwrM3_it) >= epsilon), 1)' -1;
Tave_opt_maxM3_it = sum(double(10*log10(MSE_opt_maxM3_it) >= epsilon), 1)' -1;
Tave_optM3_it = sum(double(10*log10(MSE_optM3_it) >= epsilon), 1)' -1;

Tave_mh_bound = log(sqrt(10^(epsilon/10))) ./ log(l2_mh);
Tave_opt_bound = log(sqrt(10^(epsilon/10))) ./ log(l2_opt);

Tave_mh_pwrM3_bound = log(sqrt(10^(epsilon/10))) ./ log(l2_mh_est_pwrM3);
Tave_mh_maxM3_bound = log(sqrt(10^(epsilon/10))) ./ log(l2_mh_est_maxM3);
Tave_mhM3_bound = log(sqrt(10^(epsilon/10))) ./ log(l2_mhM3);

Tave_opt_pwrM3_bound = log(sqrt(10^(epsilon/10))) ./ log(l2_opt_est_pwrM3);
Tave_opt_maxM3_bound = log(sqrt(10^(epsilon/10))) ./ log(l2_opt_est_maxM3);
Tave_optM3_bound = log(sqrt(10^(epsilon/10))) ./ log(l2_optM3);

MSE_mh_SHT_it = MSE_mh_it(Tsund_haj);
MSE_opt_SHT_it = MSE_opt_it(Tsund_haj);
MSE_mh_maxM3_SHT_it = MSE_mh_maxM3_it(Tsund_haj);
MSE_mhM3_SHT_it = MSE_mhM3_it(Tsund_haj);
MSE_opt_maxM3_SHT_it = MSE_opt_maxM3_it(Tsund_haj);
MSE_optM3_SHT_it = MSE_optM3_it(Tsund_haj);
MSE_mh_polyH_SHT_it = MSE_mh_polyH_it(Tsund_haj);
MSE_mh_polyO3_SHT_it = MSE_mh_polyO3_it(Tsund_haj);
MSE_mh_polyO5_SHT_it = MSE_mh_polyO5_it(Tsund_haj);
MSE_mh_polyO20_SHT_it = MSE_mh_polyO20_it(Tsund_haj);
for i = 1:max(Mvec)
    MSE_mh_pwrM3_SHT_it(i) = MSE_mh_pwrM3_it(Tsund_haj, i);
    MSE_opt_pwrM3_SHT_it(i) = MSE_opt_pwrM3_it(Tsund_haj, i);
end;

MSE_mh_SHT = mean(MSE_mh_SHT_it);
MSE_opt_SHT = mean(MSE_opt_SHT_it);
MSE_mh_pwrM3_SHT = mean(MSE_mh_pwrM3_SHT_it);
MSE_mh_maxM3_SHT = mean(MSE_mh_maxM3_SHT_it);
MSE_mhM3_SHT = mean(MSE_mhM3_SHT_it);
MSE_opt_pwrM3_SHT = mean(MSE_opt_pwrM3_SHT_it);
MSE_opt_maxM3_SHT = mean(MSE_opt_maxM3_SHT_it);
MSE_optM3_SHT = mean(MSE_optM3_SHT_it);
MSE_mh_polyH_SHT = mean(MSE_mh_polyH_SHT_it);
MSE_mh_polyO3_SHT = mean(MSE_mh_polyO3_SHT_it);
MSE_mh_polyO5_SHT = mean(MSE_mh_polyO5_SHT_it);
MSE_mh_polyO20_SHT = mean(MSE_mh_polyO20_SHT_it);

filename = strcat('cons_CIRC_Lin_POLY_N', num2str(Nvec));
save(filename);
