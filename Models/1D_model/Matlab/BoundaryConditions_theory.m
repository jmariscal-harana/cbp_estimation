%%	Boundary conditions (BCs) for all domains that can be prescribed in Nektar1D
%	Input
%	-PATHS
%	-REF
%	-PARAM
%	-MESH
%
%	Output
%	-BCS: boundary conditions
%
%==========================================================================
%	Jorge Mariscal-Harana, King's College London (original author)
%	v1.0 (08/05/19)
%
%==========================================================================

function [BCS] = BoundaryConditions_theory(PATHS,REF,PARAM,MESH,BCS)
%%	Scaling factors for pressure and flow
Scale_P = 133.322368;
Scale_Q = 1e6;


%%	Calculate optimal BCs for each terminal location
Artery_terminal_name	= {'Desc_Ao_I'; 'Brachio'; 'L_Carotid'; 'L_Subclavian'};
Artery_terminal_loc		= MESH.E_terminal;

addpath('~/Haemodynamic_Tools/Version6/Others//Wk3_opt_Aramburu/')
addpath('~/Haemodynamic_Tools/Version6/Others/PlotSave')

%	Optimised Wk3 calculation to get R_1, R_2, C_T
for jj = 1:length(Artery_terminal_loc)	
	%	Scale for numerical reasons within Aramburu's algorithm
	warning('Scaling units for numerical reasons within Optimise_Wk3.m')
	
	%Load P and Q waveforms
	eval(['Q = REF.',Artery_terminal_name{jj},'.ONE_CYCLE.Q*Scale_Q;']);	%[mL/s]
	eval(['P = REF.',Artery_terminal_name{jj},'.ONE_CYCLE.P/Scale_P;']);	%[mmHg]
	eval(['t = REF.',Artery_terminal_name{jj},'.ONE_CYCLE.t;']);			%[s]
	
	%Calculate R_T (= R_1 + R_2)
	P_out = REF.P_out/Scale_P;				%[mmHg]
	R_T	= (mean(P) - P_out) / mean(Q);		%[mmHg*s/mL]
	
	%Initial guess for R_1, R_2, and C_T
	R_1 = BCS.R1_Wk(Artery_terminal_loc(jj))/(Scale_P*Scale_Q);		%[mmHg s/mL]
	R_2 = R_T - R_1;
	if Artery_terminal_loc(jj) == 6
		C_T	= 0.5*BCS.C_Wk(Artery_terminal_loc(jj))*(Scale_P*Scale_Q);	%[mL/mmHg]
	else
		C_T	= BCS.C_Wk(Artery_terminal_loc(jj))*(Scale_P*Scale_Q);	%[mL/mmHg]
	end

	%Run Aramburu's optimisation
	[R_1_opt,R_2_opt,C_T_opt,P_est]	= Optimise_Wk3(t,Q,P,P_out,R_2,C_T,1,1);
	
	%Check that R_T = R_1 + R_2
	if 100*(R_T / (R_1_opt + R_2_opt) - 1) > 1, error('R_T and R_1_opt + R_2_opt differ by more than 1\%'), end
	
	%Print initial and optimal values
	fprintf('Initial R_1: %.3f; Optimised R_1: %.3f \n',R_1,R_1_opt)
	fprintf('Initial C_T: %.3f; Optimised C_T: %.3f \n',C_T,C_T_opt)
	
	%Save updated values of R_1, R_2, and C
	BCS.R_Wk_tot(Artery_terminal_loc(jj))	= (R_1_opt + R_2_opt)*(Scale_P*Scale_Q);
	BCS.C_Wk(Artery_terminal_loc(jj))		= C_T_opt/(Scale_P*Scale_Q);
	
	figure
	hold on
	plot(t,P,'k','linewidth',2)
	plot(t,P_est(end-length(t)+1:end),'b','linewidth',2)
	legend('P_{ref} [mmHg]','P_{est} [mmHg]')
	xlabel('t [s]')
	ylabel('P [mmHg]')
	set(gca,'FontSize',30)
	PlotSave(PATHS.Aortic1D_files,[PATHS.PatientName ,'_Wk3_',Artery_terminal_name{jj}])
	
end

%	Scaling factor for R1 - descending aorta Windkessel BC
if isfield(BCS,'Scale_R1')
	Scale_R1 = BCS.Scale_R1;
% 	Scale_R1 = ExtractScaleR1(REF);
	if length(Scale_R1) ~= 1, error('Length of Scale_R1 ~= 1'), end
	BCS.Scale_R1 = [1 1 1 Scale_R1 1 1 1];
	
	BCS.R2_Wk(MESH.E_terminal) = BCS.R_Wk_tot(MESH.E_terminal) - BCS.Scale_R1(MESH.E_terminal).*BCS.R1_Wk(MESH.E_terminal);
end

%	Terminal C for descending aorta Winkdessel BC
BCS.C_Wk(Artery_terminal_loc(1)) = (REF.C_T - BCS.C_Art_tot) * REF.R_T ./ BCS.R2_Wk(Artery_terminal_loc(1));

save([PATHS.Aortic1D_files,PATHS.PatientName,'_BCs'],'BCS')

end




