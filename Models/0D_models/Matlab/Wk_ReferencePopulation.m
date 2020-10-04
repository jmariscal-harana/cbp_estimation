%%	Generate Wk-2 or Wk-3 patients using literature values for R, C, PWV, etc.
%	Input:
% 	-Arterial properties: 
%		aortic characteristic impedance (Z_0)
%		outflow pressure (P_out)
%		total arterial resistance (R_T)
% 		total arterial compliance (C_T)
%	-Cardiac properties:
%		stroke volume (SV)
%		heart rate (HR)
%	
%	Output:
%	-REF: reference data for each subject: CV parameters, flow, pressure...
%
%==========================================================================
%	Jorge Mariscal-Harana, King's College London (original author)
%	v1.0 (18/10/17)
%	v2.0 (05/12/17) - REF_out contains parameters and pressure and is  
%		saved as a structure with 10 fields
%	v3.0 (11/12/17) - Wk3 pressure calculation added
%	v3.1 (15/01/18) - tau_var added to keep tau equal for Wk2 and Wk3
%		populations
%
%==========================================================================
function [REF_out] = Wk_ReferencePopulation(PATHS,REF,Inflow_baseline,Plots_REF,Windkessel_REF)
%%	List of healthy parameter variations: mean +/- 1 SD
P_out_baseline		= REF.LIT.P_out_baseline;
P_out_SD			= REF.LIT.P_out_SD;
R_T_baseline		= REF.LIT.R_T_baseline;
R_T_SD				= REF.LIT.R_T_SD;
C_T_baseline		= REF.LIT.C_T_baseline;
C_T_SD				= REF.LIT.C_T_SD;

P_out_var			= [P_out_baseline - P_out_SD, P_out_baseline - 0.5*P_out_SD, P_out_baseline, P_out_baseline + 0.5*P_out_SD, P_out_baseline + P_out_SD];	% [-10*SD -5*SD 0 +5*SD +10*SD]
R_T_var				= [R_T_baseline - R_T_SD, R_T_baseline - 0.5*R_T_SD, R_T_baseline, R_T_baseline + 0.5*R_T_SD R_T_baseline + R_T_SD];	% [-SD -0.5*SD 0 +0.5*SD +SD]
C_T_var				= [C_T_baseline - C_T_SD, C_T_baseline - 0.5*C_T_SD, C_T_baseline, C_T_baseline + 0.5*C_T_SD C_T_baseline + C_T_SD];	% [-SD -0.5*SD 0 +0.5*SD +SD]

P_0_baseline		= REF.LIT.P_0_baseline;
P_0_var				= [0.5*P_0_baseline, P_0_baseline, 1.5*P_0_baseline];

SV_baseline			= REF.LIT.SV_baseline;
SV_SD				= REF.LIT.SV_SD;
HR_baseline			= REF.LIT.HR_baseline;
HR_SD				= REF.LIT.HR_SD;
 
SV_var				= [SV_baseline - SV_SD, SV_baseline - 0.5*SV_SD, SV_baseline, SV_baseline + 0.5*SV_SD SV_baseline + SV_SD];	% [-SD -0.5*SD 0 +0.5*SD +SD]
HR_var				= [HR_baseline - HR_SD, HR_baseline - 0.5*HR_SD, HR_baseline, HR_baseline + 0.5*HR_SD HR_baseline + HR_SD];	% [-SD -0.5*SD 0 +0.5*SD +SD]

PWV_baseline		= REF.LIT.PWV_baseline;
PWV_SD				= REF.LIT.PWV_SD;
Diam_baseline		= REF.LIT.Diam_baseline;
Diam_SD				= REF.LIT.Diam_SD;

PWV_var				= [PWV_baseline - PWV_SD, PWV_baseline - 0.5*PWV_SD, PWV_baseline, PWV_baseline + 0.5*PWV_SD PWV_baseline + PWV_SD];	% [-SD -0.5*SD 0 +0.5*SD +SD]
Diam_var			= [Diam_baseline - Diam_SD, Diam_baseline - 0.5*Diam_SD, Diam_baseline, Diam_baseline + 0.5*Diam_SD Diam_baseline + Diam_SD];	% [-SD -0.5*SD 0 +0.5*SD +SD]

Density				= REF.LIT.Density;

%	Negative correlation between cfPWV and diameter - 10.1038/hr.2014.101
Z_0_var				= Density*PWV_var./(pi.*(fliplr(Diam_var)/2).^2);	%[Pa*s/m3]
Z_0_var				= Z_0_var/(133.322368*10^6);						%[mmHg*s/mL]


%%	Variation codes for flow and pressure waveforms
Variation = ['A', 'B', '0', 'Y', 'Z'];	% [A; B; ... 0; ... Y; Z] [largest decrease; baseline value; largest increase]


%%	Generate inflow waveform

if Plots_REF == 1
	figure
	hold on
end

for cardiac1 = 1:5	% 0000-X000 - SV
	for cardiac2 = 3	% 0000-0X00 - HR
		for cardiac3 = 3:3	% 0000-00X0 - Available
			for cardiac4 = 3:3	% 0000-000X - Available				
				INFLOW.dt = 1e-3;

				INFLOW.Inflow	= [Inflow_baseline,'_',Variation([cardiac1,cardiac2,cardiac3,cardiac4])];
				
				INFLOW.SV		= SV_var(cardiac1);
				INFLOW.HR		= HR_var(cardiac2);
				
				addpath([PATHS.Root,'Others/GenerateInflow/'])
				INFLOW			= GenerateInflow(PATHS,INFLOW,Plots_REF);
				LVET(cardiac1,cardiac2) = INFLOW.LVET;
				
				clear INFLOW	
			end
		end
	end
end

if Plots_REF == 1
	hold off
	
	set(gca,'XTick',[0:0.2:1.2]); xlim([0 1.135]);
	set(gca,'YTick',[0:200:1000]); ylim([-100 800]);
	xlabel('Time [s]')
% 	ylabel('Flow [mL/s]')
	set(gca, 'FontSize', 30)
	
	set(gcf,'PaperOrientation','landscape');
	set(gcf,'PaperUnits','normalized');
	set(gcf,'PaperPosition', [0 0 1 1]);
	
	%	Inflow waveform plots for:
	if cardiac1*cardiac2 == 25
		%	A) ALL variations
		% 	title('SV and HR variations','FontSize',16)
		saveas(gcf, [PATHS.Figures,'Inflow_all'], 'png')
		saveas(gcf, [PATHS.Figures,'Inflow_all'], 'pdf')
		
	elseif cardiac1 == 5
		%	B) SV variations
		title('Stroke Volume Variations','FontSize',30)
		hline = findobj(gcf, 'type', 'line');
		set(hline(1),'LineStyle','--','Color','red','LineWidth',1)
		set(hline(2),'LineStyle','-.','Color','red','LineWidth',1)
		set(hline(4),'LineStyle','-.','Color','blue','LineWidth',1)
		set(hline(5),'LineStyle','--','Color','blue','LineWidth',1)
% 		legend('\mu - \sigma','\mu - 0.5\sigma','\mu','\mu + 0.5\sigma','\mu + \sigma')

		saveas(gcf, [PATHS.Figures,'Inflow_SV'], 'png')
		saveas(gcf, [PATHS.Figures,'Inflow_SV'], 'pdf')
					
	elseif cardiac2 == 5
		%	C) HR variations
		title('Heart Rate Variations','FontSize',30)
		hline = findobj(gcf, 'type', 'line');
		set(hline(1),'LineStyle','--','Color','red','LineWidth',1)
		set(hline(2),'LineStyle','-.','Color','red','LineWidth',1)
		set(hline(4),'LineStyle','-.','Color','blue','LineWidth',1)
		set(hline(5),'LineStyle','--','Color','blue','LineWidth',1)
		legend('\mu - \sigma','\mu - 0.5\sigma','\mu','\mu + 0.5\sigma','\mu + \sigma')
		
		saveas(gcf, [PATHS.Figures,'Inflow_HR'], 'png')
		saveas(gcf, [PATHS.Figures,'Inflow_HR'], 'pdf')
		
	end
end


%%	Calculate pressure waveform for the virtual patients
fileID = fopen([PATHS.Input,Windkessel_REF,'_patients_all.txt'],'w');

k = 1;
REF_out{k} = {};

switch Windkessel_REF
	case 'Wk2'
		Z_0 = 3;
	case 'Wk3'
		Z_0 = 1:5;
	otherwise
		error('Choose a valid 0-D model')
end

for arterial1 = Z_0	% X000-0000 - Z_0
	for arterial2 = 1:5	% 0X00-0000 - P_out
		for arterial3 = 1:5	% 00X0-0000 - R_T
			for arterial4 = 1:5	% 000X-0000 - C_T
				for cardiac1 = 1:5	% 0000-X000 - SV
					for cardiac2 = 1:5	% 0000-0X00 - HR
						for cardiac3 = 3:3	% 0000-00X0 - Available
							for cardiac4 = 3:3	% 0000-000X - Available
								REF.INFO.ID	= [Windkessel_REF,'_',Variation([arterial1,arterial2,arterial3,arterial4]),'-',...
									Variation([cardiac1,cardiac2,cardiac3,cardiac4])];
								REF.INFO.Inflow	= [Inflow_baseline,'_',Variation([cardiac1,cardiac2,cardiac3,cardiac4])];
								
								REF.P_0		= P_0_baseline;

								REF.Z_0		= Z_0_var(arterial1);
								REF.P_out	= P_out_var(arterial2);
								REF.R_T		= R_T_var(arterial3);
								REF.C_T		= C_T_var(arterial4);
								
								REF.SV		= SV_var(cardiac1);
								REF.HR		= HR_var(cardiac2);
								REF.LVET	= LVET(cardiac1,cardiac2);
								
								REF_out{k}	= Wk_ReferencePressure(PATHS,REF,Windkessel_REF);	% Calculate Wk patients' pressures
								
								REF_out{k}	= setstructfields(REF,REF_out{k});	% Keep parameters used to generate patients, update new ones
								ID{k}		= REF.INFO.ID;
								k = k + 1;
							end
						end
					end
				end
			end
		end
	end
end

REF_out = [REF_out{:}];

fprintf(fileID,'%s\n',ID{:});
fclose(fileID);


end