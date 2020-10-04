clear
close all
clc
format compact
dbstop if error
disp('Example_3Wk.m')
disp('...')

fig_num = 0;
linewidth = 1.5;
markersize = 10;

load example.mat % Cargamos los vectores t, Q_meas y P_meas

Est_old_P = [38; 49; 0.32]; %[Pout; P0; tau] Estimación inicial
t_ini = 271; % consideramos el rango de valores de t=0.4 a t=0.76
t_end = length(t);
[Pout, P0, tau, Est_hist_Pout] = Estimate_Pout_P0_and_tau(t, t_ini, t_end, P_meas, Est_old_P); %Pout in mmHg, tau in s.

Est_old = [0.3; 4.0]; %Estimación inicial para R2 y C
Lim_t = [0 0.3; 0 5];
Lim   = [0 0.3; 0 5];
node = 1;
[R1, R2, C, Est_hist_Wk, Cost, fig_num] = Estimate_R2_and_C(t, P_meas, Q_meas, Pout, Est_old, Lim_t, Lim, node, fig_num);

[meanP_meas1, meanP_est1, fig_num] = Plot_Estimated_Windkessel(t, P_meas, Q_meas, R1, R2, C, Pout, fig_num);
[meanP_meas2, meanP_est2, fig_num] = Plot_Estimated_WindkesselB(t, P_meas, Q_meas, R1, R2, C, Pout, fig_num);

disp(' ')
disp('------- 3-element WK -----(Estimated: R2 and C)------------')
disp(['R1 = ' num2str(R1) ', R2 = ' num2str(R2) ', C = ' num2str(C) ', Pout = ' num2str(Pout) '.'])
disp(['Cost = +-' num2str(sqrt(Cost/length(t))) ' mmHg/measurement.'])


