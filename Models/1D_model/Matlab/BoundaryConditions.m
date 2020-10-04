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
%	v1.0 (16/10/18)
%	v2.0 (09/11/18) - introduced scaling factor for R2
%
%==========================================================================

function [BCS] = BoundaryConditions(PATHS,REF,PARAM,MESH,BCS)
%%	Inflow boundary conditions - commonly at the inlet of the ascending aorta
addpath('~/Haemodynamic_Tools/Version6/Others/GenerateInflow_v2')
t = REF.Asc_Ao.ONE_CYCLE.t;
Q = REF.Asc_Ao.ONE_CYCLE.Q;
[t_H,Q_H,T,HR,SV] = oneDbio_GenerateInflowFile_v2([PATHS.REF_data,'Inflow/'],REF.INFO.ID,0,t,Q);

% figure
% plot(t_H,Q_H)
% title('Harmonics flow')

%	Copy boundary conditions to the folder for 1-D Aortic input files
CopyBcsFile(PATHS,REF)


%%	Outflow boundary conditions - commonly at the outlet of terminal domains
%	Arterial compliance and terminal R1
BCS = ArterialC_TerminalR1(PARAM,MESH,BCS);		
R1_Wk(MESH.E_terminal) = BCS.R1_Wk(MESH.E_terminal);	%remove non-terminal R1 values
BCS.R1_Wk = R1_Wk;

if ~isfield(REF.Desc_Ao_I,'ONE_CYCLE')	%NOTE: is this ever gonna happen?
	%	Terminal outlet area distribution - non-available descending aorta flow
	BCS = TerminalArea_MESH(MESH,BCS);
elseif ~isfield(REF.Brachio,'ONE_CYCLE') || ~isfield(REF.L_Carotid,'ONE_CYCLE') || ~isfield(REF.L_Subclavian,'ONE_CYCLE')
	%	Terminal outlet area distribution - non-available supra-aortic flow
	BCS = AreaDistribution_HIS(REF,MESH,BCS);
else
	% 	Terminal outflow distribution - available flows
	BCS = OutflowDistribution(REF,MESH,BCS);
end
BCS.R_Wk_tot(MESH.E_terminal) = REF.R_T ./ BCS.OutflowDistribution(MESH.E_terminal);

%	Terminal R2 for each Winkdessel BC
BCS.R2_Wk(MESH.E_terminal) = BCS.R_Wk_tot(MESH.E_terminal) - BCS.R1_Wk(MESH.E_terminal);

%	Scaling factor for R1 - descending aorta windkessel
if isfield(REF,'Scale_R1')
	Scale_R1 = REF.Scale_R1;
% 	Scale_R1 = ExtractScaleR1(REF);
	if length(Scale_R1) ~= 1, error('Length of Scale_R1 ~= 1'), end
	BCS.Scale_R1 = [1 1 1 Scale_R1 1 1 1];
	
	BCS.R2_Wk(MESH.E_terminal) = BCS.R_Wk_tot(MESH.E_terminal) - BCS.Scale_R1(MESH.E_terminal).*BCS.R1_Wk(MESH.E_terminal);
end

%	Scaling factor for R2 - every windkessel
if isfield(REF,'Scale_R2')
	BCS.Scale_R2 = REF.Scale_R2;
	BCS.R2_Wk(MESH.E_terminal) = BCS.Scale_R2 * BCS.R2_Wk(MESH.E_terminal);
	BCS.R_Wk_tot(MESH.E_terminal) = BCS.R2_Wk(MESH.E_terminal) + BCS.R1_Wk(MESH.E_terminal);
	tau = REF.R_T*REF.C_T;
	REF.R_T = 1/sum(1./BCS.R_Wk_tot(MESH.E_terminal));
	REF.C_T = tau/REF.R_T;
end

%	Terminal C for each Winkdessel BC
BCS.C_Art_tot = sum(BCS.C_Art);
if BCS.C_Art_tot > REF.C_T
	error('Arterial compliance > Total compliance')
end
BCS.C_Per_tot = (REF.C_T - BCS.C_Art_tot);
% BCS.C_Wk(MESH.E_terminal) =  BCS.C_Per_tot * BCS.OutflowDistribution(MESH.E_terminal) .*...
% 	(BCS.R1_Wk(MESH.E_terminal) + BCS.R2_Wk(MESH.E_terminal)) ./ BCS.R2_Wk(MESH.E_terminal);
BCS.C_Wk(MESH.E_terminal) =  BCS.C_Per_tot * REF.R_T ./ BCS.R2_Wk(MESH.E_terminal);


end


%%	Function definitions
%	Copy .BCS file to the right destination
function CopyBcsFile(PATHS,REF)

copyfile([PATHS.REF_data,'Inflow/',REF.INFO.ID,'_IN_1.bcs'],...
	[PATHS.Aortic1D_files,PATHS.PatientName,'_IN_1.bcs'])

end

function CopyBcsFile_AoCo(PATHS,REF)

copyfile([PATHS.REF_data,'Flow_Asc_Ao/MRI_inflow_',REF.INFO.ID{1},'_IN_1.bcs'],...
	[PATHS.Aortic1D_files,PATHS.PatientName,'_IN_1.bcs'])

end

%	Concatenate time and flow waveforms for the duration of the current simulation
function ExtendInflowBC(PATHS)
%	Original data
data		= dlmread([PATHS.Inflow,T_name]);
t			= data(:,1);
dt			= t(2) - t(1);
T			= t(end);
Q			= data(:,2);

if T ~= T_cycle
	error('Cardiac cycle durations do not coincide: check original .txt file')
end

%	Numerical parameters
Duration	= INPUT.duration;
N_cycles	= ceil(Duration/t(end));

%	Original data concatenated for the full duration
t_extended	= [t(1):dt:Duration]';
dummy		= [Q;repmat(Q(2:end),N_cycles-1,1)];
Q_extended	= dummy(1:length(t_extended));

%	Interpolated data
dt_sim		= INPUT.dt;
t_sim		= [t(1):dt_sim:t(end)]';
pp = spline(t, Q);    % creates a cubic-splines interpolator
Q_sim = ppval(pp, t_sim);

%	Interpolated data concatenated for the full duration
t_extended_sim	= [t_sim(1):dt_sim:Duration]';
dummy_sim		= [Q_sim;repmat(Q_sim(2:end),N_cycles-1,1)];
Q_extended_sim	= dummy_sim(1:length(t_extended_sim));

%	Interpolation plot
figure
hold on
plot(t_extended,Q_extended,'k')
plot(t_extended_sim,Q_extended_sim,'.b')
plot(t_extended_sim,Q_extended_sim,'.b','MarkerSize',20)
hold off

OutputFile	= [PATHS.Data, T_name];
% dlmwrite(OutputFile,[t_extended Q_extended],'delimiter',' ','precision',10);
dlmwrite(OutputFile,[t_extended_sim Q_extended_sim],'delimiter',' ','precision',10);
end

function [BCS] = ArterialC_TerminalR1(PARAM,MESH,BCS)
% Calculate average properties of each arterial segment and outlet R1
for d=1:MESH.nbSeg
%     BCS.c_prox(d) = sqrt(2/3/PARAM.Rho{1}*(0.1*(MESH.Eh_k1*exp(MESH.Eh_k2*100*MESH.Rad_prox(d))+MESH.Eh_k3)));
%     BCS.c_dist(d) = sqrt(2/3/PARAM.Rho{1}*(0.1*f(MESH.Eh_k1*exp(MESH.Eh_k2*100*MESH.Rad_dist(d))+MESH.Eh_k3)));
    BCS.c_avg(d) = 0.5*(MESH.c_prox(d)+MESH.c_dist(d));					%average of mean proximal and mean distal PWV
    BCS.A_avg(d) = 0.5*pi*(MESH.Rad_prox(d)^2+MESH.Rad_dist(d)^2);		%average of mean proximal and mean distal areas
    BCS.C_Art(d) = BCS.A_avg(d)*MESH.L(d)/PARAM.Rho{1}/BCS.c_avg(d)^2;	%average segment compliance
    BCS.R1_Wk(d) = PARAM.Rho{1}*MESH.c_dist(d)/pi/MESH.Rad_dist(d)^2;	%Wk characteristic impedance at terminal outlet
end

end

function [BCS] = TerminalArea_MESH(MESH,BCS)

A_tot = 0;  % initialise total area
for d=1:MESH.nbSeg
	if sum(MESH.E_terminal == d) == 1
		A_terminal(d) = pi*MESH.Rad_dist(d)^2;	%outlet area
		A_tot = A_tot + A_terminal(d);
	else
		A_terminal(d) = 0;
		A_tot = A_tot;
	end
end

% for d=1:MESH.nbSeg
% 	if sum(MESH.E_terminal == d) == 1
% 		A_terminal(d) = 0.5*pi*(MESH.Rad_prox(d)^2+MESH.Rad_dist(d)^2);	%average area
% 		A_tot = A_tot + A_terminal(d);
% 	else
% 		A_terminal(d) = 0;
% 		A_tot = A_tot;
% 	end
% end

BCS.A_tot					= A_tot;
BCS.A_terminal				= A_terminal;
BCS.OutflowDistribution		= A_terminal/A_tot;

end

function [BCS] = AreaDistribution_HIS(REF,MESH,BCS)

Q_mean_in = mean(REF.Asc_Ao.ONE_CYCLE.Q);				%mean cardiac inflow 
Q_mean_out_Desc_Ao = mean(REF.Desc_Ao_I.ONE_CYCLE.Q);	%mean descending aorta outflow

if Q_mean_in < Q_mean_out_Desc_Ao, error('Mean descending aorta outflow > mean inflow'), end

BranchFlow = 1 - Q_mean_out_Desc_Ao/Q_mean_in;

BranchOutletAreas = [mean(REF.Brachio.A_out)
	mean(REF.L_Carotid.A_out)
	mean(REF.L_Subclavian.A_out)];

BCS.OutflowDistribution(MESH.E_terminal(1)) = Q_mean_out_Desc_Ao/Q_mean_in;
BCS.OutflowDistribution(MESH.E_terminal(2:end)) = BranchFlow * (BranchOutletAreas ./ sum(BranchOutletAreas));

end

function [BCS] = OutflowDistribution(REF,MESH,BCS)
%%	Extract terminal mean flows from waveform data
Q_mean_in = mean(REF.Asc_Ao.ONE_CYCLE.Q);		%mean cardiac inflow 
Q_mean_out = [mean(REF.Desc_Ao_I.ONE_CYCLE.Q)
	mean(REF.Brachio.ONE_CYCLE.Q)
	mean(REF.L_Carotid.ONE_CYCLE.Q)
	mean(REF.L_Subclavian.ONE_CYCLE.Q)];		%individual mean branch outflow

if (Q_mean_in - sum(Q_mean_out))/Q_mean_in*100 > 5
	error('Difference between mean inflow and total mean outflow > 5%')
else
	BCS.OutflowDistribution(MESH.E_terminal) = Q_mean_out/sum(Q_mean_out);	%according to total mean outflow
end

end

function [Scale_R1] = ExtractScaleR1(REF)
S.v = REF.Pressure;
S.fs = round(1/(REF.Asc_Ao.ONE_CYCLE.t(2)-REF.Asc_Ao.ONE_CYCLE.t(1)));
options.do_plot = 0;

addpath('/Users/joh15/Haemodynamic_Tools/Version6/Others/PulseAnalyse/')
[~, fid_pts, peaks, filtered_signal] = PulseAnalyse5(S, options);

%	Normalise waveform
P_filtered = filtered_signal.v;
P_filtered = P_filtered - min(P_filtered);
P_filtered = P_filtered/max(P_filtered);
P_si_filtered = P_filtered(fid_pts.dia);

%	Calculate second peak pressure 
P_si = min(REF.Pressure) + P_si_filtered*(max(REF.Pressure) - min(REF.Pressure));
P_si = REF.Pressure(fid_pts.dia);

P_out = REF.P_out;
P_s = max(REF.Pressure);

%	Equation for Scale_R1 from?10.1098/rsif.2016.0073
Scale_R1 = (P_si - P_out) / (2*P_s - P_si - P_out);

end
