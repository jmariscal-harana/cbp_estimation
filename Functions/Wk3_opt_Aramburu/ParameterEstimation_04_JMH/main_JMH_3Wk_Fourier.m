clear all
close all
clc
format compact
dbstop if error
disp('main_JMH_3Wk_Fourier.m')
disp('...')

fig_num = 0;
linewidth = 1.5;
markersize = 10;

load Subject_1.mat % Load vectors t (s), Q (ml/s) and P (mmHg)

%% 3-WK parameter estimation
Pout = 33; %(mmHg)
 
R2_0 = 0.7*0.9; C_0 = 1.77;
Est_old = [R2_0; C_0]; % Initial guess for R2 y C, which should be sufficiently close to the solution.
Lim_t = [0.5066 0.7066; 1 3];
Lim   = [0.5066 0.7066; 1 3];
node = 1;
[R1, R2, C, Est_hist_Wk, Cost, fig_num] = Estimate_R2_and_C(t, P, Q, Pout, Est_old, Lim_t, Lim, node, fig_num);
Rt = (mean(P)-Pout)/mean(Q);
R1_0 = Rt - R2_0;

[meanP_meas2, meanP_est2, fig_num] = Plot_Estimated_WindkesselB(t, P, Q, R1_0, R2_0, C_0, Pout, fig_num);
[meanP_meas3, meanP_est3, fig_num] = Plot_Estimated_WindkesselB(t, P, Q, R1, R2, C, Pout, fig_num);

disp(' ')
disp('------- 3-element WK -----(Estimated: R2 and C)------------')
disp(['R1 = ' num2str(R1) ', R2 = ' num2str(R2) ', C = ' num2str(C) ', Pout = ' num2str(Pout) '.'])
disp(['Cost = +-' num2str(sqrt(Cost/length(t))) ' mmHg/measurement.'])

%% Fourier analysis
Pfft = fft(P);
Pz = abs(Pfft(1:10));
Qfft = fft(Q);
Qz = abs(Qfft(1:10));
Modulus = Pz./Qz;
fig_num  = fig_num + 1;
subplot(2,1,1)
plot(Modulus,'-ob','MarkerFaceColor','b')
axis([0 10 0 2])
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10])
set(gca,'XTickLabel',{'0','1','2','3','4','5','6','7','8','9'})
set(gca,'fontsize',14)
ylabel('Modulus Z [mmHg s/mL]')
box off
%- Phase
Pp = angle(Pfft(1:10));
Qp = angle(Qfft(1:10));
Phase = Pp - Qp;
subplot(2,1,2)
plot(Phase,'-ob','MarkerFaceColor','b')
axis([0 10 -10 10])
set(gca,'XTick',[1 2 3 4 5 6 7 8 9 10])
set(gca,'XTickLabel',{'0','1','2','3','4','5','6','7','8','9'})
line([0 10],[0 0],'linestyle','- -','color','k')
set(gca,'fontsize',14)
xlabel('Frequency [Hz]','fontsize',14)
ylabel('Phase \Phi [rad]')
box off