%%	Extract REF data from in silico database using info not available in vivo

function [REF] = Aortic1D_ExtractREF_theory(REF,Dataset)

switch Dataset
		
	case 'pwdb_6_reference'		
		%%	Arterial segments of the 116-artery model
		REF.Asc_Ao.Segment				= 1;
		REF.Ao_Arch_I.Segment			= 2;
		REF.Brachio.Segment				= 3;
		REF.Ao_Arch_II.Segment			= 14;
		REF.L_Carotid.Segment			= 15;
		REF.Desc_Ao_I.Segment			= 18;
		REF.L_Subclavian.Segment		= 19;
		REF.L_Brachial.Segment			= 21;
		
		
		%%	Extract geometry		
		REF = ExtractGeometry_InSilico(REF,Dataset);
		
		
		%%	Extract DBP
		REF.DBP = min(REF.Asc_Ao.ONE_CYCLE.P);

		
end

end

function [REF] = ExtractParameters_periodTEX(PATHS,REF)

Tex_File = fileread([PATHS.REF_data,REF.INFO.ID{1},'_period.tex']);

R_T_string = {'Peripheral vascular resistance: &', ' Pa s m$^{-3}$\\'};
C_T_string = {'Total compliance: &', ' m$^{3}$ Pa$^{-1}$ \\'};
tau_string = {'Time constant: &', ' s \\'};
PWV_string = {'$\rightarrow$', '&', '$' };
rho_string = {'Density (Kg/m$^{3}$): &', 'Blood viscosity (mPa s):  &'};

%	R_T	- sum (in parallel) of Windkessel resistances
R_T_start	= strfind(Tex_File, R_T_string{1}) + length(R_T_string{1});
R_T_end		= R_T_start + strfind(Tex_File(R_T_start:end), R_T_string{2}) - 1;
R_T_temp	= Tex_File(R_T_start:R_T_end);
REF.R_T		= str2double(R_T_temp);

%	C_T	- total compliance (arterial + impedance weighted peripheral compliance)
C_T_start	= strfind(Tex_File, C_T_string{1}) + length(C_T_string{1});
C_T_end		= C_T_start + strfind(Tex_File(C_T_start:end), C_T_string{2}) - 1;
C_T_temp	= Tex_File(C_T_start:C_T_end);
REF.C_T		= str2double(C_T_temp);

%	tau - time constant
tau_start	= strfind(Tex_File, tau_string{1}) + length(tau_string{1});
tau_end		= tau_start + strfind(Tex_File(tau_start:end), tau_string{2}) - 1;
tau_temp	= Tex_File(tau_start:tau_end);
REF.tau		= str2double(tau_temp);

%	PWV - mean PWV at ascending aorta (inlet)
PWV_locs	= strfind(Tex_File, PWV_string{1});
PWV_start	= PWV_locs(1) + length(PWV_string{1});
PWV_end		= PWV_locs(2);
PWV_temp	= Tex_File(PWV_start:PWV_end);
PWV_start	= strfind(PWV_temp, PWV_string{2}) + length(PWV_string{2});
PWV_end		= strfind(PWV_temp, PWV_string{3}) - length(PWV_string{3});
PWV_temp	= PWV_temp(PWV_start:PWV_end);
REF.PWV		= str2double(PWV_temp);

%	rho - density
rho_start	= strfind(Tex_File, rho_string{1}) + length(rho_string{1});
rho_end		= rho_start + strfind(Tex_File(rho_start:end), rho_string{2}) - 6;
rho_temp	= Tex_File(rho_start:rho_end);
REF.rho		= str2double(rho_temp);

end

function [REF] = ExtractGeometry_InSilico(REF,Dataset)

warning('All inputs should be SI units')

switch Dataset
		
	case 'pwdb_6_reference'
		%	Arterial segments of the 116-artery model
		Art_Domains_116_art = [REF.Asc_Ao.Segment
			REF.Ao_Arch_I.Segment
			REF.Ao_Arch_II.Segment
			REF.Desc_Ao_I.Segment
			REF.Brachio.Segment
			REF.L_Carotid.Segment
			REF.L_Subclavian.Segment]';
 		
 		%	Arterial segments of the aortic model (must be in the same order as the 55-artery model segments)
		Art_Domains_Aortic1D	= [1 2 3 4 5 6 7];
		Art_Domains_Aorta		= [1 2 3 4];
		
		%	Mesh connectivity for given aortic model segments
		REF.Nodes	= [1 2;
			2 3;
			3 4;
			4 5;
			2 6;
			3 7;
			4 8];	%[inlet, outlet]
		
		%	This should follow the same order as the mesh connectivity and geometry 
		Artery_name			= REF.INFO.Artery_geo_name;

		for jj = 1:length(Artery_name)
			eval(['Length(jj)	= REF.',Artery_name{jj},'.length;']);
			eval(['R_in(jj)		= REF.',Artery_name{jj},'.R_in_mesh;']);
			eval(['R_out(jj)	= REF.',Artery_name{jj},'.R_out_mesh;']);
		end
		
		
		%%	Mesh dimensions
		REF.L			= Length;	%[m]
		REF.Rad_prox	= R_in;		%mean radius; [m]
		REF.Rad_dist	= R_out;	%mean radius; [m]
		
		
		%%	Theoretical local PWV
		%	Extract local PWV using Olufsen's formula
		k				= REF.k;	%[g, s, cm]
		rho				= REF.rho;
		
		%	PWV at inlets
		r				= R_in;
		c_squared		= 2/(3*rho) * 0.1 * (k(1) * exp(100 * k(2) .* r) + k(3));	%Mynard 2015 [SI units]
		c				= sqrt(c_squared);
		c_in			= c;

		%	PWV at outlets
		r				= R_out;
		c_squared		= 2/(3*rho) * 0.1 * (k(1) * exp(100 * k(2) .* r) + k(3));	%Mynard 2015 [SI units]
		c				= sqrt(c_squared);
		c_out			= c;
		
		REF.c_prox		= c_in;
		REF.c_dist		= c_out;
		
		
end

end

function [REF] = ExtractArea_InSilico(PATHS,REF)
%%	Extract patient area for each desired arterial location
Artery_name = {'Asc_Ao', 'Ao_Arch_I', 'Ao_Arch_II', 'Desc_Ao', 'Brachio', 'L_Carotid', 'L_Subclavian'};
His_location = [3, 3, 3, 3, 3, 3, 3];

up.dir = PATHS.REF_data;
for jj = 1:length(Artery_name)
	eval(['Segment_no(jj) = REF.',Artery_name{jj},'.Segment;'])
	eval(['up.filename{jj} = ''',REF.INFO.ID{1},'_',num2str(Segment_no(jj)),'.his'';']);
end

Filename_modified =	REF.INFO.ID{1};
Filename_modified(strfind(Filename_modified,'-')) = '_';


%%	Extract last beat only (for comparison)
up.all_beats	= false;
data			= ConvertHistoryFiles(up);

for jj = 1:length(Artery_name)
	eval(['REF.',Artery_name{jj},'.ONE_CYCLE.A	= data.sim_',Filename_modified,'(jj).A(:,',num2str(His_location(jj)),');']);
	eval(['fs						= data.sim_',Filename_modified,'(jj).fs;'])
	eval(['REF.',Artery_name{jj},'.ONE_CYCLE.t	= [0:1/fs:(length(data.sim_',Filename_modified,'(jj).A)-1)/fs]'';'])
	eval(['REF.',Artery_name{jj},'.ONE_CYCLE.start	= data.sim_',Filename_modified,'(jj).start_sample(',num2str(His_location(jj)),')'';'])
end


%%	Extract all beats
up.all_beats	= true;
data			= ConvertHistoryFiles(up);

for jj = 1:length(Artery_name)
	eval(['REF.',Artery_name{jj},'.CYCLES.A	= data.sim_',Filename_modified,'(jj).A(:,',num2str(His_location(jj)),');']);
	eval(['fs						= data.sim_',Filename_modified,'(jj).fs;'])
	eval(['REF.',Artery_name{jj},'.CYCLES.t	= [0:1/fs:(length(data.sim_',Filename_modified,'(jj).A)-1)/fs]'';'])
end

end

function [REF] = ExtractFlow_InSilico(PATHS,REF)
%%	Extract patient flow for each desired arterial location
Artery_name = {'Asc_Ao', 'Desc_Ao', 'Brachio', 'L_Carotid', 'L_Subclavian'};
His_location = [1, 3, 3, 3, 3];

up.dir = PATHS.REF_data;
for jj = 1:length(Artery_name)
	eval(['Segment_no(jj) = REF.',Artery_name{jj},'.Segment;'])
	eval(['up.filename{jj}		= ''',REF.INFO.ID{1},'_',num2str(Segment_no(jj)),'.his'';']);
end

Filename_modified =	REF.INFO.ID{1};
Filename_modified(strfind(Filename_modified,'-')) = '_';


%%	Extract last beat only (for comparison)
up.all_beats	= false;
data			= ConvertHistoryFiles(up);

for jj = 1:length(Artery_name)
	eval(['REF.',Artery_name{jj},'.ONE_CYCLE.Q	= data.sim_',Filename_modified,'(jj).Q(:,',num2str(His_location(jj)),');']);
	eval(['fs						= data.sim_',Filename_modified,'(jj).fs;'])
	eval(['REF.',Artery_name{jj},'.ONE_CYCLE.t	= [0:1/fs:(length(data.sim_',Filename_modified,'(jj).Q)-1)/fs]'';'])
	eval(['REF.',Artery_name{jj},'.ONE_CYCLE.start	= data.sim_',Filename_modified,'(jj).start_sample(',num2str(His_location(jj)),')'';'])
end


%%	Extract all beats
up.all_beats	= true;
data			= ConvertHistoryFiles(up);

for jj = 1:length(Artery_name)
	eval(['REF.',Artery_name{jj},'.CYCLES.Q	= data.sim_',Filename_modified,'(jj).Q(:,',num2str(His_location(jj)),');']);
	eval(['fs						= data.sim_',Filename_modified,'(jj).fs;'])
	eval(['REF.',Artery_name{jj},'.CYCLES.t	= [0:1/fs:(length(data.sim_',Filename_modified,'(jj).Q)-1)/fs]'';'])
end

end

function [REF] = ExtractPressure_InSilico(PATHS,REF)
%%	Extract patient pressure for each desired arterial location
Artery_name = {'Asc_Ao', 'L_Brachial', 'L_Carotid'};
His_location = [1, 3, 2];

up.dir = PATHS.REF_data;
for jj = 1:length(Artery_name)
	eval(['Segment_no(jj) = REF.',Artery_name{jj},'.Segment;'])
	eval(['up.filename{jj}		= ''',REF.INFO.ID{1},'_',num2str(Segment_no(jj)),'.his'';']);
end

Filename_modified =	REF.INFO.ID{1};
Filename_modified(strfind(Filename_modified,'-')) = '_';


%%	Extract last beat only (for comparison)
up.all_beats	= false;
data			= ConvertHistoryFiles(up);

for jj = 1:length(Artery_name)
	eval(['REF.',Artery_name{jj},'.ONE_CYCLE.P	= data.sim_',Filename_modified,'(jj).P(:,',num2str(His_location(jj)),');']);
	eval(['fs						= data.sim_',Filename_modified,'(jj).fs;'])
	eval(['REF.',Artery_name{jj},'.ONE_CYCLE.t	= [0:1/fs:(length(data.sim_',Filename_modified,'(jj).P)-1)/fs]'';'])
	eval(['REF.',Artery_name{jj},'.ONE_CYCLE.start	= data.sim_',Filename_modified,'(jj).start_sample(',num2str(His_location(jj)),')'';'])
end


%%	Extract all beats
up.all_beats	= true;
data			= ConvertHistoryFiles(up);

for jj = 1:length(Artery_name)
	eval(['REF.',Artery_name{jj},'.CYCLES.P	= data.sim_',Filename_modified,'(jj).P(:,',num2str(His_location(jj)),');']);
	eval(['fs						= data.sim_',Filename_modified,'(jj).fs;'])
	eval(['REF.',Artery_name{jj},'.CYCLES.t	= [0:1/fs:(length(data.sim_',Filename_modified,'(jj).P)-1)/fs]'';'])
end

end

function [REF] = ExtractPWV_InSilico(REF)
addpath('~/Haemodynamic_Tools/Version6/Others/PWV_TT')

Wave_type = 2;	%flow
SamplingRate = round(1/(REF.Asc_Ao.ONE_CYCLE.t(2)-REF.Asc_Ao.ONE_CYCLE.t(1)));

%	Read waveform pair
Wave_in = REF.Asc_Ao.ONE_CYCLE.Q;
Wave_out = REF.Desc_Ao_I.ONE_CYCLE.Q;

%	Ensure equal vector length
[Wave_in, Wave_out] = CorrectSize(Wave_in, Wave_out);

%	Increase number of cycles
Wave_in = [Wave_in; repmat(Wave_in(2:end),4,1)];
Wave_out = [Wave_out; repmat(Wave_out(2:end),4,1)];

%	Identify starting time index for each waveform
Start_in = REF.Asc_Ao.ONE_CYCLE.start;
Start_out = REF.Desc_Ao_I.ONE_CYCLE.start;
while Start_out < Start_in
	Start_out = Start_out + length(REF.Asc_Ao.ONE_CYCLE.Q) - 1;
end

Wave_in = Wave_in(Start_out - Start_in + 1 : end);
[Wave_in, Wave_out] = CorrectSize(Wave_in, Wave_out);

%	Calculate transit time (TT)
TT = TTAlgorithm([Wave_in';Wave_out'],SamplingRate,1,Wave_type,1,0);

REF.PWV = REF.Aortic_L/mean(TT);

end

function [REF] = Extract_t_sys(REF,Pressure,time)

%	Method 1 - Aortic Q
% addpath('~/Haemodynamic_Tools/Version6/Others/t_sys/')
% REF.HAEMO.t_sys = t_sys_from_Q(REF.Asc_Ao.ONE_CYCLE.Q,time,0);

%	Method 2 - From P
addpath('~/Haemodynamic_Tools/Version6/Others/PulseAnalyse/')
Wave.v = Pressure;
Wave.fs = round((length(time) - 1) / time(end));	%(Number of data points - 1)/Period
Options.do_plot = false;
[~,fid_pts] = PulseAnalyse5(Wave,Options);
t_sys_index = fid_pts.dic;
% t_sys_index = fid_pts.dia;
REF.HAEMO.t_sys = time(t_sys_index);

end

function [REF] = ExtractParameters_IN(PATHS,REF)
%% P_out from .in
FileName_IN = [PATHS.REF_data,REF.INFO.ID{1},'.in'];

Text		= fileread(FileName_IN);

P_out_string = {'Viscosity', ' pinf'};	%text before and after the imposed P_out value

%	P_out
P_out_start	= strfind(Text, P_out_string{1}) + length(P_out_string{1});
P_out_end	= P_out_start + strfind(Text(P_out_start:end), P_out_string{2}) - 1;
P_out_temp	= Text(P_out_start:P_out_end);
REF.P_out	=  str2double(P_out_temp);

end

function [Wave_in, Wave_out] = CorrectSize(Wave_in, Wave_out)

if length(Wave_in) > length(Wave_out)
	Wave_in = Wave_in(1:length(Wave_out));
elseif length(Wave_in) < length(Wave_out)
	Wave_out = Wave_out(1:length(Wave_in));
end

end
