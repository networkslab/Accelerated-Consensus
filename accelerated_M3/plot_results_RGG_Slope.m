% Script for plotting the figures that depict the mean squred error
% performance comparison and the averaging time comparison for the setting
% where topology is modelled as a Random Geometric Graph and the
% initialization is "Slope". Plot all results, excluding polynomial filter
% and Sundaram-Hadjicostis.


clc;


Nsim = [25; 50; 100; 150; 200];
Nsim=[25; 50; 100; 150; 200];
% Nsim = 200;
in_folder = 'sim_results_2/cons_RG_Lin_POLY_N';
out_folder = 'fig3/cons_RG_Lin_poly_paper_';
FontSize_lab = 24;
FontSize_ax = 18;
MarkerSize = 8;

DFvec = [2; 2; 2; 3; 3];
DFbvec = [50; 200; 1000; 500; 1000];
shift_mag_vec = round(DFvec / 2);

mse_lb = -100;
epsilon = -100;

Tave_mh = zeros(size(Nsim));
Tave_mhM3 = zeros(size(Nsim));
Tave_mh_pwrM3 = zeros(size(Nsim));
Tave_opt = zeros(size(Nsim));
Tave_optM3 = zeros(size(Nsim));

for ii=1:length(Nsim)
    
    in_filename = strcat(in_folder, num2str(Nsim(ii)), '.mat');
    load(in_filename, 'Mc', 'MSE_mh_it', ...
                      'MSE_opt_it', 'MSE_mh_pwrM3_it', ...
                      'MSE_mh_pwrM3_it', 'MSE_mh_maxM3_it', ...
                      'MSE_mhM3_it', 'MSE_opt_pwrM3_it', ...
                      'MSE_opt_maxM3_it', 'MSE_optM3_it');
                  
    load(in_filename, 'Mc', 'MSE_mh', ...
                      'MSE_opt', 'MSE_mh_pwrM3', ...
                      'MSE_mh_pwrM3', 'MSE_mh_maxM3', ...
                      'MSE_mhM3', 'MSE_opt_pwrM3', ...
                      'MSE_opt_maxM3', 'MSE_optM3');

    x_ub = 100; % 3*sum(10*log10(MSE_mhM3) >= mse_lb);
    
    DF = DFvec(ii);
    DFb = DFbvec(ii);
    shift_mag = shift_mag_vec(ii);
 
    h=0;
    hold off;
    h(1)=semilogy(downsample((1:Mc)', round(DF*3)), ...
        downsample(sqrt(MSE_mh), round(DF*3)), '-^b');
    hold on;
    h(2)=semilogy(downsample((1:Mc)', round(DF*3), 3*shift_mag), ...
        downsample(sqrt(MSE_opt), round(DF*3), 3*shift_mag), '-+k');
    h(3)=semilogy(downsample((1:Mc)', 2*DF), ...
        downsample(sqrt(MSE_mhM3), 2*DF), '-dr');
    h(4)=semilogy(downsample((1:Mc)', 2*DF, 2*shift_mag), ...
        downsample(sqrt(MSE_mh_pwrM3), 2*DF, 2*shift_mag), '-.xg');
    h(5)=semilogy(downsample((1:Mc)', DF), ...
        downsample(sqrt(MSE_optM3), DF), '-sk');
    
    
    legend('MH', 'Opt', 'MH-Proposed', 'MH-ProposedEst', 'Opt-Proposed', ...
        'location', 'NorthEast');
    set(gca, 'FontSize', FontSize_ax); 
    set(h, 'LineWidth', 2);
    set(h, 'MarkerSize', MarkerSize);
    set(h(4), 'MarkerSize', 1.25*MarkerSize);
    set(h(2), 'MarkerSize', 1.25*MarkerSize);
    ylim([10^(mse_lb/20); 1]);
    xlim([1; x_ub]);
    xlabel('{\it t}', 'FontSize', FontSize_lab);
    ylabel('|| x(t) - \mu ||_2', 'FontSize', FontSize_lab);
    set(gca,'YTick', [1e-5; 1e-4; 1e-3; 1e-2; 1e-1; 1]);
    
    out_filename = strcat(out_folder, 'MSE_N', num2str(Nsim(ii)));
    saveas(gcf, out_filename, 'epsc2');
    delete(gcf);
    
    
    Tave_mh_it = sum(double(10*log10(MSE_mh_it) >= epsilon), 1)' -1;
    Tave_opt_it = sum(double(10*log10(MSE_opt_it) >= epsilon), 1)' -1;

    Tave_mh_pwrM3_it = sum(single(10*log10(MSE_mh_pwrM3_it) >= epsilon), 1)' -1;
    Tave_mh_maxM3_it = sum(single(10*log10(MSE_mh_maxM3_it) >= epsilon), 1)' -1;
    Tave_mhM3_it = sum(double(10*log10(MSE_mhM3_it) >= epsilon), 1)' -1;

    Tave_opt_pwrM3_it = sum(single(10*log10(MSE_opt_pwrM3_it) >= epsilon), 1)' -1;
    Tave_opt_maxM3_it = sum(single(10*log10(MSE_opt_maxM3_it) >= epsilon), 1)' -1;
    Tave_optM3_it = sum(double(10*log10(MSE_optM3_it) >= epsilon), 1)' -1;

    Tave_mh(ii) = mean(Tave_mh_it);    
    Tave_mhM3(ii) = mean(Tave_mhM3_it);
    Tave_mh_pwrM3(ii) = mean(Tave_mh_pwrM3_it);
    
    Tave_opt(ii) = mean(Tave_opt_it);
    Tave_optM3(ii) = mean(Tave_optM3_it);
    
    a=1;
end;
delete(gcf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=0;
hold off;
h(1)=plot(Nsim, Tave_mh, '-b^');
hold on;
h(2)=plot(Nsim, Tave_opt, '-+k');
h(3)=plot(Nsim, Tave_mhM3, '-dr');
h(4)=plot(Nsim, Tave_mh_pwrM3, '-.xg');
h(5)=plot(Nsim, Tave_optM3, '-sk');

legend('MH', 'Opt', 'MH-Proposed', 'MH-ProposedEst', 'Opt-Proposed', 'location', 'best');
set(gca, 'FontSize', FontSize_ax); 
set(h, 'LineWidth', 2);
set(h, 'MarkerSize', MarkerSize);
set(h(4), 'MarkerSize', 1.5*MarkerSize);
set(h(2), 'MarkerSize', 1.5*MarkerSize);
ylim([min(1); max(Tave_mh)]); % Tave_optM3
xlim([min(Nsim)-1; max(Nsim)]);
xlabel('{\it N}', 'FontSize', FontSize_lab);
ylabel('T_{ave}', 'FontSize', FontSize_lab);

out_filename = strcat(out_folder, 'TIME_TRUE');
saveas(gcf, out_filename, 'epsc2');
delete(gcf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=0;
hold off;
h(1)=plot(Nsim, Tave_mh./Tave_mhM3, '-dr');
hold on;
h(2)=plot(Nsim, Tave_mh./Tave_mh_pwrM3, '-.xg');
h(3)=plot(Nsim, Tave_mh./Tave_optM3, '-sk');
h(4)=plot(Nsim, Tave_mh./Tave_opt, '-+k');

legend('MH-ProposedEst', 'MH-Proposed', 'Opt-Proposed', 'Opt', 'location', 'best');
set(gca, 'FontSize', FontSize_ax); 
set(h, 'LineWidth', 2);
set(h, 'MarkerSize', MarkerSize);
set(h(2), 'MarkerSize', 1.5*MarkerSize);
set(h(4), 'MarkerSize', 1.5*MarkerSize);

ylim([1; max(Tave_mh./Tave_optM3)]);
xlim([min(Nsim)-0.25; max(Nsim)]);
xlabel('{\it N}', 'FontSize', FontSize_lab);
ylabel('T_{ave}(W) / T_{ave}(\Phi)', 'FontSize', FontSize_lab);

out_filename = strcat(out_folder, 'TIME_RATIO');
saveas(gcf, out_filename, 'epsc2');
delete(gcf);
