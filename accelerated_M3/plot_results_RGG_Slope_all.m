% Script for plotting the figures that depict the mean squred error
% performance comparison and the averaging time comparison for the setting
% where topology is modelled as a Random Geometric Graph and the
% initialization is "Linear". Plot all results.

clc;

Nsim = [25; 50; 100; 150; 200];
Nsim=[25; 50; 100; 150; 200];
% Nsim = 200;
in_folder = 'sim_results_2/cons_RG_Lin_POLY_N';
out_folder = 'fig3/cons_RG_Lin_poly_acc_paper_';
FontSize_lab = 24;
FontSize_ax = 18;
MarkerSize = 8;
DF = 10;
DFb = 50;
shift_mag = 5;

DFvec = [4; 4; 8; 8; 8];
DFbvec = [50; 200; 1000; 500; 1000];
shift_mag_vec = floor(DFvec / 2);

mse_lb = -100;
epsilon = -100; % dB

Tave_mh = zeros(size(Nsim));
Tave_mhM3 = zeros(size(Nsim));
Tave_mh_pwrM3 = zeros(size(Nsim));
Tave_mh_maxM3 = zeros(size(Nsim));
Tave_opt = zeros(size(Nsim));
Tave_mh_polyH = zeros(size(Nsim));
Tave_mh_polyO3 = zeros(size(Nsim));
Tave_mh_polyO5 = zeros(size(Nsim));
Tave_mh_polyO20 = zeros(size(Nsim));

MSE_opt_SHT_vec = zeros(size(Nsim));
MSE_mhM3_SHT_vec = zeros(size(Nsim));
MSE_mh_pwrM3_SHT_vec = zeros(size(Nsim));
MSE_mh_maxM3_SHT_vec = zeros(size(Nsim));
MSE_mh_polyH_SHT_vec = zeros(size(Nsim));
MSE_mh_polyO3_SHT_vec = zeros(size(Nsim));
MSE_mh_polyO5_SHT_vec = zeros(size(Nsim));
MSE_mh_polyO20_SHT_vec = zeros(size(Nsim));


for ii=1:length(Nsim)
    
    in_filename = strcat(in_folder, num2str(Nsim(ii)), '.mat');
                  
    load(in_filename, 'Mc', 'MSE_mh_it', ...
                      'MSE_opt_it', 'MSE_mhM3_it', ...
                      'MSE_mh_pwrM3_it', 'MSE_mh_maxM3_it', ...
                      'MSE_mh_polyH_it', 'MSE_mh_polyO3_it', ...
                      'MSE_mh_polyO5_it', 'MSE_mh_polyO20_it');
                  
%     load(in_filename, 'Mc', 'MSE_mh', ...
%                       'MSE_opt', 'MSE_mhM3', ...
%                       'MSE_mh_pwrM3', 'MSE_mh_maxM3', ...
%                       'MSE_mh_polyH', 'MSE_mh_polyO3', ...
%                       'MSE_mh_polyO5', 'MSE_mh_polyO20');
                  
    [MSE_mh, MSE_mh_it] = calc_mean_MSE(MSE_mh_it);
    [MSE_opt, MSE_opt_it] = calc_mean_MSE(MSE_opt_it);
    [MSE_mhM3, MSE_mhM3_it] = calc_mean_MSE(MSE_mhM3_it);
    [MSE_mh_pwrM3, MSE_mh_pwrM3_it] = calc_mean_MSE(MSE_mh_pwrM3_it);
    [MSE_mh_maxM3, MSE_mh_maxM3_it] = calc_mean_MSE(MSE_mh_maxM3_it);
    [MSE_mh_polyH, MSE_mh_polyH_it] = calc_mean_MSE(MSE_mh_polyH_it);
    [MSE_mh_polyO3, MSE_mh_polyO3_it] = calc_mean_MSE(MSE_mh_polyO3_it);
    [MSE_mh_polyO5, MSE_mh_polyO5_it] = calc_mean_MSE(MSE_mh_polyO5_it);
    [MSE_mh_polyO20, MSE_mh_polyO20_it] = calc_mean_MSE(MSE_mh_polyO20_it);
    
    load(in_filename, 'Mc', ...
                      'MSE_opt_SHT', 'MSE_mhM3_SHT', ...
                      'MSE_mh_pwrM3_SHT', 'MSE_mh_maxM3_SHT', ...
                      'MSE_mh_polyH_SHT', 'MSE_mh_polyO3_SHT', ...
                      'MSE_mh_polyO5_SHT', 'MSE_mh_polyO20_SHT');
                  
    MSE_opt_SHT_vec(ii) = MSE_opt_SHT;
    MSE_mhM3_SHT_vec(ii) = MSE_mhM3_SHT;
    MSE_mh_pwrM3_SHT_vec(ii) = MSE_mh_pwrM3_SHT;
    MSE_mh_maxM3_SHT_vec(ii) = MSE_mh_maxM3_SHT;
    MSE_mh_polyH_SHT_vec(ii) = MSE_mh_polyH_SHT;
    MSE_mh_polyO3_SHT_vec(ii) = MSE_mh_polyO3_SHT;
    MSE_mh_polyO5_SHT_vec(ii) = MSE_mh_polyO5_SHT;
    MSE_mh_polyO20_SHT_vec(ii) = MSE_mh_polyO20_SHT;
                         
    x_ub = 2*calc_Tave(MSE_mh_polyO20, mse_lb); % sum(10*log10(MSE_mh_polyO20) >= mse_lb);
    
    DF = DFvec(ii);
    DFb = DFbvec(ii);
    shift_mag = shift_mag_vec(ii);
 
    h=0;
    hold off;
    h(1)=semilogy(downsample((1:Mc)', DF, shift_mag), ...
        downsample(sqrt(MSE_opt), DF, shift_mag), '-k+');
    hold on;
    h(2)=semilogy(downsample((1:Mc)', DF, shift_mag), ...
        downsample(sqrt(MSE_mh_polyO3), DF, shift_mag), 'cv');
    h(3)=semilogy(downsample((1:Mc)', DF), ...
        downsample(sqrt(MSE_mhM3), DF), '-rd');
    h(4)=semilogy(downsample((1:Mc)', DF, shift_mag), ...
        downsample(sqrt(MSE_mh_pwrM3), DF, shift_mag), '-.xg');
    h(5)=semilogy(downsample((1:Mc)', DF, shift_mag), ...
        downsample(sqrt(MSE_mh_polyO20), DF, shift_mag), 'b>');
    
    legend('Opt', 'MH-PolyFilt3', 'MH-Proposed', 'MH-ProposedEst', 'MH-PolyFilt7', 'location', 'NorthEast');
    
    h(6)=semilogy((1:Mc)', sqrt(MSE_mh_polyO3), '-c');
    h(7)=semilogy((1:Mc)', sqrt(MSE_mh_polyO20), '-b');
    
    set(gca, 'FontSize', FontSize_ax); 
    set(h, 'LineWidth', 2);
    set(h, 'MarkerSize', MarkerSize);
    set(h(4), 'MarkerSize', 1.5*MarkerSize);
    set(h(1), 'MarkerSize', 1.25*MarkerSize);
    ylim([10^(mse_lb/20); 1]);
    xlim([1; x_ub]);
    xlabel('{\it t}', 'FontSize', FontSize_lab);
    ylabel('|| x(t) - \mu ||_2', 'FontSize', FontSize_lab);
    set(gca,'YTick', [1e-5; 1e-4; 1e-3; 1e-2; 1e-1; 1]);
    
    out_filename = strcat(out_folder, 'MSE_N', num2str(Nsim(ii)));
    saveas(gcf, out_filename, 'epsc2');
    delete(gcf);
    
    
    Tave_mh_it = calc_Tave(MSE_mh_it, epsilon);
    Tave_opt_it = calc_Tave(MSE_opt_it, epsilon);

    Tave_mh_pwrM3_it = calc_Tave(MSE_mh_pwrM3_it, epsilon);
    Tave_mh_maxM3_it = calc_Tave(MSE_mh_maxM3_it, epsilon);
    Tave_mhM3_it = calc_Tave(MSE_mhM3_it, epsilon);
        
    Tave_mh_polyH_it = calc_Tave(MSE_mh_polyH_it, epsilon);
    Tave_mh_polyO3_it = calc_Tave(MSE_mh_polyO3_it, epsilon);
    Tave_mh_polyO5_it = calc_Tave(MSE_mh_polyO5_it, epsilon);
    Tave_mh_polyO20_it = calc_Tave(MSE_mh_polyO20_it, epsilon);
    
    Tave_mh(ii) = mean(Tave_mh_it);
    Tave_opt(ii) = mean(Tave_opt_it);

    Tave_mh_pwrM3(ii) = mean(Tave_mh_pwrM3_it);
    Tave_mh_maxM3(ii) = mean(Tave_mh_maxM3_it);
    Tave_mhM3(ii) = mean(Tave_mhM3_it);
        
    Tave_mh_polyH(ii) = mean(Tave_mh_polyH_it);
    Tave_mh_polyO3(ii) = mean(Tave_mh_polyO3_it);
    Tave_mh_polyO5(ii) = mean(Tave_mh_polyO5_it);
    Tave_mh_polyO20(ii) = mean(Tave_mh_polyO20_it);
    
end;
delete(gcf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=0;
hold off;
h(1)=plot(Nsim, Tave_mh_polyO3, '-vc');
hold on;
h(2)=plot(Nsim, Tave_opt, '-k+');
h(3)=plot(Nsim, Tave_mh_polyO20, '-b>');
h(4)=plot(Nsim, Tave_mhM3, '-rd');

legend('MH-PolyFilt3', 'Opt', ...
       'MH-PolyFilt7', 'MH-Proposed', 'location', 'NorthWest');
set(gca, 'FontSize', FontSize_ax); 
set(h, 'LineWidth', 2);
set(h, 'MarkerSize', MarkerSize);
set(h(3), 'MarkerSize', 1.25*MarkerSize);
ylim([min(Tave_mh_polyO20); max(Tave_mh_polyO3)]);
xlim([min(Nsim)-1; max(Nsim)]);
xlabel('{\it N}', 'FontSize', FontSize_lab);
ylabel('T_{ave}', 'FontSize', FontSize_lab);

out_filename = strcat(out_folder, 'TIME_TRUE');
saveas(gcf, out_filename, 'epsc2');
delete(gcf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=0;
hold off;
h(1)=plot(Nsim, Tave_mh ./ Tave_mhM3, '-rd');
hold on;
h(2)=plot(Nsim, Tave_mh ./ Tave_mh_polyO20, '-b>');
h(3)=plot(Nsim, Tave_mh ./ Tave_opt, '-+k');
h(4)=plot(Nsim, Tave_mh ./ Tave_mh_polyO3, '-vc');
legend('MH-Proposed', 'MH-PolyFilt7', 'Opt', 'MH-PolyFilt3', 'location', 'NorthWest');
set(gca, 'FontSize', FontSize_ax); 
set(h, 'LineWidth', 2);
set(h, 'MarkerSize', MarkerSize);

ylim([1; max(Tave_mh./Tave_mhM3)+0.1]);
xlim([min(Nsim)-1; max(Nsim)]);
xlabel('{\it N}', 'FontSize', FontSize_lab);
ylabel('T_{ave}(W) / T_{ave}(\Phi)', 'FontSize', FontSize_lab);

out_filename = strcat(out_folder, 'TIME_RATIO');
saveas(gcf, out_filename, 'epsc2');
delete(gcf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=0;
hold off;
h(1)=semilogy(Nsim, sqrt(MSE_mh_polyO3_SHT_vec), '-vc');
hold on;
h(2)=semilogy(Nsim, sqrt(MSE_opt_SHT_vec), '-+k');
h(3)=semilogy(Nsim, sqrt(MSE_mh_polyO20_SHT_vec), '->b');
h(4)=semilogy(Nsim, sqrt(MSE_mhM3_SHT_vec), '-dr');

legend('MH-PolyFilt3', 'Opt', 'MH-PolyFilt7', ...
       'MH-Proposed', 'location', 'NorthEast');
set(gca, 'FontSize', FontSize_ax); 
set(h, 'LineWidth', 2);
set(h, 'MarkerSize', MarkerSize);
set(h(2), 'MarkerSize', 1.25*MarkerSize);

% ylim([sqrt(min(MSE_mhM3_SHT_vec)); sqrt(max(MSE_mh_polyO3_SHT_vec))]);
xlim([min(Nsim)-1; max(Nsim)]);
xlabel('{\it N}', 'FontSize', FontSize_lab);
ylabel('|| x - \mu ||_2', 'FontSize', FontSize_lab);

    
out_filename = strcat(out_folder, 'MSE_SHT');
saveas(gcf, out_filename, 'epsc2');
delete(gcf);
