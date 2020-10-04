%%	Plot dataset characteristics
format compact;
close all;
clear;
clc;


%%	PATHS
PATHS.REF_data			= '~/Haemodynamic_Tools/Version6/VirtualDB/55_art/Nektar_outputFiles/Simulations/2018_12_06/';


%%	Load REF - 55_art
Population = '55_art';
load([PATHS.REF_data,Population,'_reference.mat'],'REF')


%%	Select patients to plot
%	Patients whose R_T (aortic root) < R_T_Wk (parallel sum of Windkessel resistances)
Patient_ID = find([REF.R_T]-[REF.R_T_Wk] < 0);


%%	Calculate relevant characteristics
for jj = 1:length(REF)
	%	R_T difference
	R_T_diff(jj) = REF(jj).R_T-REF(jj).R_T_Wk;
	
	%	Blood pressure values
	SBP(jj) = max(REF(jj).Asc_Ao.ONE_CYCLE.P);
	DBP(jj) = min(REF(jj).Asc_Ao.ONE_CYCLE.P);
	PP(jj) = SBP(jj) - DBP(jj);
	MBP_Asc(jj) = mean(REF(jj).Asc_Ao.ONE_CYCLE.P);
	MBP_Desc(jj) = mean(REF(jj).Desc_Ao.ONE_CYCLE.P);
end

for jj = 1:length(Patient_ID)
	%	Patient
	Patient{jj} = REF(Patient_ID(jj)).INFO.ID{1};
end


%%	Plots
%	Problematic patients
figure
hold on
for jj = Patient_ID
	plot(REF(jj).Asc_Ao.ONE_CYCLE.t,REF(jj).Asc_Ao.ONE_CYCLE.P)
	xlabel('Time [s]')
	ylabel('Pressure - ascending aorta [mmHg]')
	set(gca,'FontSize',20)
end
hold off
legend(num2str(Patient_ID'))


%	R_T difference
figure
hold on
plot([REF.R_T],'o')
plot([REF.R_T_Wk],'o')
plot(R_T_diff,'o')
plot(Patient_ID,R_T_diff(Patient_ID),'x','MarkerSize',20)
hold off
xlabel('Patient number')
ylabel('R_T [mmHg*s/mL]')
legend('R_T','R_{T,Wk}','R_T - R_{T,Wk}','R_T - R_{T,Wk} (problematic)')
set(gca,'FontSize',20)


%	C_T_arterial
figure
hold on
plot([REF.C_T_art],'o')
plot(Patient_ID,[REF(Patient_ID).C_T_art],'x','MarkerSize',20)
hold off
xlabel('Patient number')
ylabel('Arterial Compliance: post-simulation [mL/mmHg]')
set(gca,'FontSize',20)


%	PP
figure
hold on
plot(PP,'o')
plot(Patient_ID,PP(Patient_ID),'x','MarkerSize',20)
hold off
xlabel('Patient number')
ylabel('Pulse pressure [mmHg]')
set(gca,'FontSize',20)


%	MBP
figure
hold on
plot(Patient_ID,MBP_Asc(Patient_ID),'o')
plot(Patient_ID,MBP_Desc(Patient_ID),'x')
hold off
xlabel('Patient number')
ylabel('Mean pressure [mmHg]')
legend('Ascending aorta','Descending aorta')
set(gca,'FontSize',20)