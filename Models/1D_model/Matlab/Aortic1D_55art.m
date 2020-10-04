%%	Specify required data to generate Nektar1D input files (.IN) for the Aortic1D model
%	Input
%	-PATHS: folders where I/O data is stored
%	-REF: reference parameters for the Aortic1D model
%
%	Output
%	-.IN and .BCS containing patient-specific files for Aortic1D simulations
%	-.MAT containing values required to generate .IN and .BCS files
%
%	For more information, open most recent READ_ME.DOCX file
%
%==========================================================================
%	Jorge Mariscal-Harana, King's College London (original author)
%	v1.0 (08/10/18)
%
%==========================================================================format compact;
format compact;
close all;
clear;
clc;


%%	PATHS
PATHS.Aortic1D			= '~/Haemodynamic_Tools/Version6/Aortic1D/';		%Root
PATHS.Matlab			= [PATHS.Aortic1D,'Matlab/'];						%Matlab functions
PATHS.Aortic1D_files	= [PATHS.Aortic1D,'Nektar_inputFiles/'];			%.BCS, .IN and .MAT aortic model files are saved here

addpath(PATHS.Matlab)


%%	Generate REF data for input file 
%	REF.INFO.ID should describe parameter variations for each patient. 
%	NOTE: Your inflow file (.BCS) should have the same ID.
%	EXAMPLE:
%	REF(1).INFO.ID = "0000"; will generate an .IN file named 
%	"Aortic_0000.in" (baseline subject), so you should name your .BCS file 
%	"Aortic_0000_IN_1.bcs"

nb = {'A', 'B', '0', 'Y', 'Z'};

Scale_R2 = [0.8, 0.9, 1.0, 1.1, 1.2]; 

%	Give suitable format to input data 
jj = 1;
for i1 = [1 3 5]%Scale_R2
	for i2 = [1 3 5]%Cardiac period
		for i3 = [1 3 5] %Stroke Volume
			for i4 = [1 3 5] %Duration of systole
				REF(jj).INFO.ID = {[nb{i1},nb{i2},nb{i3},nb{i4}]};
				REF(jj).Scale_R2 = Scale_R2(i1);
				
				REF_temp(jj) = BaselineData(REF(jj));
				jj = jj + 1;
			end
		end
	end
end

REF = REF_temp;
clear REF_temp


%%	Generate an aortic model for each patient
for jj = 1:length(REF)
	PATHS.PatientName = ['Aortic_',REF(jj).INFO.ID{1}];	
	Aortic1D_Input_File(PATHS,REF(jj))
end

restoredefaultpath
%%%%%%%%%%%%%%%%%%


%%	Function definitions
function [REF] = BaselineData(REF)
%%	Baseline geometry from 55-artery baseline subject
REF = BaselineGeometry(REF);


%%	Baseline mean flow for the calculation of flow distributions
REF.Asc_Ao.ONE_CYCLE.Q_mean = 1.016498141109771e-04;	%[m3/s]
REF.Desc_Ao.ONE_CYCLE.Q_mean = 7.24943368834285e-05;	%[m3/s]
REF.Brachiocephalic.ONE_CYCLE.Q_mean = 1.45721138505747e-05;	%[m3/s]
REF.L_Carotid.ONE_CYCLE.Q_mean = 6.30066960000000e-06;	%[m3/s]
REF.L_Subclavian.ONE_CYCLE.Q_mean = 7.83138482285714e-06;	%[m3/s]


%%	Baseline PWV
%	Regional (aorta)
REF.PWV = 5.961301423879362;	%[m/s]


%%	Baseline P_out, R_T, C_T
%	NOTE: R_T will change when changing R2
REF.P_out = 3.999671000000000e+03;	%[Pa]
REF.R_T = 7.572443529532567e+07;	%[Pa s/m3]
REF.C_T = 8.350193299450490e-09;	%[m3/Pa]


end

function [REF] = BaselineGeometry(REF)
%	Arterial segments of the aortic model
Art_Domains_Aortic1D	= [1 2 3 4 5 6 7];	%equivalent to [1 2 14 18 3 15 19] in the 55-artery model
Art_Domains_Aorta		= [1 2 3 4];		%required to calculate aortic length
Art_Domains_Branches	= [5 6 7];			%not required atm

%	Mesh connectivity for given aortic model segments
REF.nodes		= [	1 2;
					2 3;
					3 4;
					4 5;
					2 6;
					3 7;
					4 8];	%[inlet, outlet]

%	Mesh dimensions
REF.L			= [0.0580000000000000,0.0230000000000000,0.0450000000000000,0.0600000000000000,0.0390000000000000,0.160000000000000,0.0390000000000000];	%[m]
REF.Rad_prox	= [0.0184000000000000,0.0154000000000000,0.0130000000000000,0.0120000000000000,0.0122000000000000,0.00439000000000000,0.00634000000000000];	%mean radius; [m]
REF.Rad_dist	= [0.0183000000000000,0.0147000000000000,0.0126000000000000,0.0113000000000000,0.0108000000000000,0.00370000000000000,0.00482000000000000];	%mean radius; [m]
REF.c_prox		= [5.66000000000000,5.89000000000000,6.13000000000000,6.24000000000000,6.23000000000000,8.02000000000000,7.30000000000000];			%PWV at mean radius; [m/s]
REF.c_dist		= [5.67000000000000,5.96000000000000,6.18000000000000,6.33000000000000,6.41000000000000,8.39000000000000,7.83000000000000];		%PWV at mean radius; [m/s]

%	Aortic length for f2f PWV calculation between aortic root and descending aorta
REF.Aortic_L	= sum(REF.L(Art_Domains_Aorta));


end



