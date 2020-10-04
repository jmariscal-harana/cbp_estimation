% Create_parametric_inflow.m
% Solve the system of equations describing the inlet flow as piecewise
% polynomial
%  Q1 (0-tsystole): 4th order
%  Q2 (tsystole - tstar): 2nd order
%  Q3 (tstar - tdiastole) (=dicrotic notch): 4the order
%  Q4 (tdiastole - tend):  constant=0

% ====================
% Marie Willemet
% KCL, October 2015
% ====================
%
% Modified by Jordi Alastruey
% January 2017
%

function [INFLOW] = Create_parametric_inflow(PATHS,INFLOW)

age = INFLOW.age;
dt	= INFLOW.dt;
QM	= INFLOW.CO;
HR	= INFLOW.HR;
T	= INFLOW.T;
SR	= round(1/dt);
ReverseFlow	= INFLOW.reverseflow;


%% Obtain parametric curve
Qmax	= INFLOW.Qmax;	  %Systolic peak flow
Qi		= INFLOW.Qi;      %Flow at inflection point
Qi2		= INFLOW.Qi2;     %Flow at inflection point in systolic decrease (old)
Qmin	= INFLOW.Qmin;    %Minimum flow during dicrotic notch

ti_0	= INFLOW.ti_0;   %Time of inflection point in systolic increase
tmax_0	= INFLOW.tmax_0; %Time of systolic peak
ti2_0	= INFLOW.ti2_0;  %Time of inflection point in systolic decrease (old only)
ts_0	= INFLOW.ts_0;   %Time of start of dicrotic notch
tmin_0	= INFLOW.tmin_0; %Time of minimum flow during dicrotic notch
td_0	= INFLOW.td_0;   %Time of start of diastole

%==== multiply systole period duration
mult_period=INFLOW.mult_period;

%=== multiply distribution of times of systole
mult_sysPeak=INFLOW.mult_sysPeak;
mult_sys=INFLOW.mult_sys;

%round ti,tmax,... to ensure right sizes of time vector
ti   = round(ti_0	*mult_sys*mult_period/dt).*dt;
tmax = round(tmax_0 *mult_sys*mult_period*mult_sysPeak/dt).*dt;
ti2  = round(ti2_0	*mult_sys*mult_period/dt).*dt;
ts   = round(ts_0	*mult_sys*mult_period/dt).*dt;
tmin = round(tmin_0	*mult_sys*mult_period/dt).*dt;
td   = round(td_0	*mult_sys*mult_period/dt).*dt;
T	 = round(T		*mult_sys*mult_period/dt).*dt;

%	Flow type: old, young, mix
Qtot_old   = Get_parametric_old(dt,Qi,Qmax,Qmin,Qi2,ti,tmax,tmin,ti2,ts,td,T);
Qtot_young = Get_parametric_young(dt,Qi,Qmax,Qmin,ti,tmax,tmin,ts,td,T);
Qtot_mix   = (Qtot_old+Qtot_young)./2;


%% Adjust the mean flow to QM
Qtot_old_QM = Qtot_old./mean(Qtot_old).*QM;
Qtot_young_QM = Qtot_young./mean(Qtot_young).*QM;
Qtot_mix_QM = Qtot_mix./mean(Qtot_mix).*QM;

%-- Secant method to adjust Qmax - Optimisation does not seem to work on old waveforms...
% Qm = 0.2; %mean(Qtot_old); %mean flow to obtain (parametric curve)
% Qmax_opt_o = secant_Qmax_old(Qm, dt, Qi, Qmin, Qi2, ti, tmax, tmin, ti2, ts, td);
% Qtot_old_opt= Get_parametric_old(dt, Qi, Qmax_opt_o, Qmin, Qi2, ti, tmax, tmin, ti2, ts, td);
% Qmax_opt_y = secant_Qmax_young(Qm, dt, Qi, Qmin, ti,tmax,tmin, ts, td);
% Qtot_young_opt = Get_parametric_young(dt, Qi, Qmax_opt_y, Qmin, ti,tmax,tmin, ts, td);


%% Eliminate reverse flow
if (ReverseFlow == 0)
	for j=1:length(Qtot_old_QM)
		if Qtot_old_QM(j)<0
			Qtot_old_QM(j) = 0;
			Qtot_young_QM(j) = 0;
			Qtot_mix_QM(j) = 0;
			
			Qtot_old(j) = 0;
			Qtot_young(j) = 0;
			Qtot_mix(j) = 0;
		else
			Qtot_old_QM(j) = Qtot_old_QM(j);
			Qtot_young_QM(j) = Qtot_young_QM(j);
			Qtot_mix_QM(j) = Qtot_mix_QM(j);
			
			Qtot_old(j) = Qtot_old(j);
			Qtot_young(j) = Qtot_young(j);
			Qtot_mix(j) = Qtot_mix(j);
		end
	end
end


%% Adjust the duration of the cardiac cycle to T=1/HR
% t = (0:dt:1)'*60/HR;
t = (0:dt:T)';
if (t(end)>T)
	Nend = round(abs(t(end)-T)/(dt*60/HR));
	tfinal = t(1:end-Nend);
	Qtot_young_QM = Qtot_young_QM(1:end-Nend);
	Qtot_old_QM = Qtot_old_QM(1:end-Nend);
	Qtot_mix_QM = Qtot_mix_QM(1:end-Nend);
else
	tfinal = t;
end

if (t(end)<T)
	Nend = round(abs(t(end)-T)/dt);
	tfinal = horzcat(t(1:end-1)',(t(end):dt:T))';
	Qtot_young_QM = horzcat(Qtot_young_QM,zeros(1,Nend));
	Qtot_old_QM = horzcat(Qtot_old_QM,zeros(1,Nend));
	Qtot_mix_QM = horzcat(Qtot_mix_QM,zeros(1,Nend));
end


%% Figures
%Comparison with O'Rourke waveforms
addpath(PATHS.Data_inflow)

old_tmp = load('orourke2009_old_modified5.dat');
old = old_tmp(2,:)'./max(old_tmp(2,:));
young_tmp = load('orourke2009_young_modified5.dat');
young = young_tmp(2,:)'./max(young_tmp(2,:));
%t = old_tmp(1,:)';

% % switch age
% %     case 'young'
% %         figure(1);
% %         set(figure(1),'Position',[0 400 560 420])
% %         set(gca,'fontsize',16)
% %         xlabel('t/tmax (/)')
% %         ylabel('Q/Qmax (/)')
% %         hold on
% %         plot(t/max(t)*HR/60,young,'r','linewidth',3);
% %         plot(t/max(t),Qtot_young,'k','linewidth',2)
% %         legend('O Rourke','parametric')
% %         plot(tmax/max(t)*60/HR,Qmax,'*k')
% %         plot(ti/max(t)*60/HR,Qi,'*k')
% %         plot(ts/max(t)*60/HR,Qtot_young(ts./dt),'*k')
% %         plot(td/max(t)*60/HR,0,'*k')
% %         title('young')
% %
% %     case 'old'
% %         figure(2);
% %         set(figure(2),'Position',[0 400 560 420])
% %         set(gca,'fontsize',16)
% %         xlabel('t/tmax (/)')
% %         ylabel('Q/Qmax (/)')
% %         hold on
% %         plot(t/max(t)*HR/60,old,'r','linewidth',3);
% %         plot(t/max(t),Qtot_old,'k','linewidth',2)
% %         legend('O Rourke','parametric')
% %         plot(tmax/max(t)*60/HR,Qmax,'*k')
% %         plot(ti/max(t)*60/HR,Qi,'*k')
% %         plot(ts/max(t)*60/HR,Qtot_old(ts./dt),'*k')
% %         plot(td/max(t)*60/HR,0,'*k')
% %         plot(ti2/max(t)*60/HR,Qi2,'*k')
% %         title('old')
% %
% %     case 'avg'
% %         figure(3);
% %         set(figure(3),'Position',[0 400 560 420])
% %         set(gca,'fontsize',16)
% %         hold on
% %         plot(t/max(t),Qtot_young_QM,'linewidth',3);
% %         plot(t/max(t),Qtot_old_QM,'linewidth',3);
% %         plot(t/max(t),Qtot_mix_QM,'linewidth',3);
% %         xlabel('t/tmax (/)')
% %         ylabel('Q/Qmax (/)')
% %         legend(['Young : ',num2str(mean(Qtot_young_QM))],['Old : ',num2str(mean(Qtot_old_QM))], ['Avg : ',num2str(mean(Qtot_mix_QM))])
% % end


if INFLOW.Plots == 1
	% Final inflow waveform in ml/s
	% figure
	% set(gca,'fontsize',16)
	% hold on
	switch age
		case 'young'
			plot(tfinal,Qtot_young_QM);
		case 'old'
			if strcmp(INFLOW.Inflow, 'Inflow_old_0000')
				plot(tfinal,Qtot_old_QM, 'k', 'LineWidth', 2);
			else
				plot(tfinal,Qtot_old_QM,'k');
			end
		case 'avg'
			plot(tfinal,Qtot_mix_QM);
	end
	% set(gca,'YTick',[0 500]); ylim([-100 600]);
	% set(gca,'XTick',[0:0.2:T]); xlim([0 T]);
	% ylabel('Q [mL/s]')
	% xlabel('t [s]')
	% title('Inflow','FontSize',16)
end


%%	Paths
OutputFile			= [PATHS.Wk_study,'Inflow/'];
OutputFile_freq		= [OutputFile,INFLOW.Inflow,'.bcs'];
OutputFile_time		= [OutputFile,INFLOW.Inflow,'_t.txt'];
OutputFile_time_m3s = [OutputFile,INFLOW.Inflow,'_t_m3s.txt'];

INFLOW.t			= tfinal;


%% Save waveforms to text file 'Inflow_{period}_{sysPeak}_age', e.g.: 'Inflow_125_80_y'
% Name = ['Inflow_',num2str(mult_period*100),'_',num2str(mult_sys*100)];
% cd([OutputPath,'_Qwaves_t/']);
% dlmwrite([Name,'_y_t.txt'],[t Qtot_young_QM'],' ');
% dlmwrite([Name,'_o_t.txt'],[t Qtot_old_QM'],' ');
% dlmwrite([Name,'_m_t.txt'],[t Qtot_mix_QM'],' ');

switch age
	case 'young'		
		INFLOW.Q	= Qtot_young_QM';
		dlmwrite(OutputFile_time,[INFLOW.t INFLOW.Q],' ');
		dlmwrite(OutputFile_time_m3s,[INFLOW.t INFLOW.Q*1e-6],' ');	%SI units
	case 'old'
		INFLOW.Q	= Qtot_old_QM';
		dlmwrite(OutputFile_time,[INFLOW.t INFLOW.Q],' ');
		dlmwrite(OutputFile_time_m3s,[INFLOW.t INFLOW.Q*1e-6],' ');	%SI units
	case 'avg'
		INFLOW.Q	= Qtot_mix_QM';
		dlmwrite(OutputFile_time,[INFLOW.t INFLOW.Q],' ');
		dlmwrite(OutputFile_time_m3s,[INFLOW.t INFLOW.Q*1e-6],' ');	%SI units
end


%% Save into frequency file
% [t_H_y, Q_H_y] = Harmonics_filter(t,Qtot_young_QM,1000,0, [OutputPath,'_Qwaves_f/', Name,'_y.bcs']);
% [t_H_o, Q_H_o] = Harmonics_filter(t,Qtot_old_QM,1000,0, [OutputPath,'_Qwaves_f/', Name,'_o.bcs']);
% [t_H_m, Q_H_m] = Harmonics_filter(t,Qtot_mix_QM,1000,0,   [OutputPath,'_Qwaves_f/', Name,'_m.bcs']);

switch age
	case 'young'
		[t_H_y, Q_H_y] = Harmonics_filter(t,Qtot_young_QM*1e-6,SR,0,OutputFile_freq);
	case 'old'
		[t_H_o, Q_H_o] = Harmonics_filter(t,Qtot_old_QM*1e-6,SR,0,OutputFile_freq);
	case 'avg'
		[t_H_m, Q_H_m] = Harmonics_filter(t,Qtot_mix_QM*1e-6,SR,0,OutputFile_freq);
end


% %% Save all files
% DirPath = [OutputPath,'/_Qwaves_t/'];
% d = dir(DirPath);
% str = {d.name};
% figure; hold on
%
% cd(DirPath)
% for i = 4:length(str)
%     dat = load(str{i});
%     plot(dat(:,1), dat(:,2))
% end


% Update inflow parameters with final values
INFLOW.T = INFLOW.t(end); % T [s]
INFLOW.HR = 60/INFLOW.T; % HR [bpm]
INFLOW.SV = trapz(INFLOW.Q)*dt;	%SV [mL]
INFLOW.CO = INFLOW.SV*INFLOW.HR/60;  % Cardiac output [mL/s]
INFLOW.LVET = td;

end