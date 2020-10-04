%%	Generates Nektar1D input files (.IN and .BCS) from reference data (in silico or clinical)
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
%==========================================================================

function Aortic1D_Input_File(PATHS,REF)
%%	1. PARAMETER LIST
%	Required
PARAM.EQTYPE		= 0;									%0: nonlinear formulation
															%1: linear formulation
UP.Duration			= 15.0;									%[s]
PARAM.DT			= 5e-5;									%[s]
UP.His_dt			= PARAM.DT*50;							%[s]
PARAM.Rho			= REF.rho;								%Blood density [Kg/m3]
PARAM.Viscosity		= 2.5e-3;								%Blood viscosity [Pa s] - Pete

%	Optional
PARAM.INTTYPE		= 2;									%Integration order of numerical scheme: 1, 2 (DEFAULT) or 3
lastCycle			= floor(UP.Duration/REF.HAEMO.T);
PARAM.T_initial		= (lastCycle-1).*REF.HAEMO.T;
PARAM.T_final		= lastCycle.*REF.HAEMO.T;
UP.mcf				= 9;									%NOT IN USE - momentum correction factor associated with the velocity profile

PARAM				= ParameterList(REF,PARAM,UP);


%%	2. MESH DEFINITION
MESH.Viscoeslastic	= 0;			%0: elastic simulations
									%1: viscoelastic simulations
MESH.MechProps		= 'Beta';		%'Beta', 'Eh' or 'Empirical', as explained in the Reference Manual	
MESH.PWVType		= 'f2f';		%'Prescribed' from .TEX or 'f2f' as average aortic PWV
MESH.p_order		= 3;
MESH.q_order		= 3;
MESH.elmt			= 0.02;			%Length of a single element in spatial discretisation [m]

MESH				= MeshDefinition(REF,PARAM,MESH);


%%	3. BOUNDARY CONDITIONS
BCS.InflowType		= 'F0';		%'F0' = harmonics
BCS.WkType			= 'W';		%Reflective terminal BCs - 3-element Windkessel
% BCS.WkType			= 'R';		%Absorbent terminal BCs - purely resistance

BCS = BoundaryConditions(PATHS,REF,PARAM,MESH,BCS);

%	Theoretical study
% BCS = BoundaryConditions_theory(PATHS,REF,PARAM,MESH,BCS);
% load([PATHS.Aortic1D_files,PATHS.PatientName,'_BCs'],'BCS')


%%	4. INITIAL CONDITIONS
ICS.Type			= 'A0 mesh';	%'A0 mesh': mesh area; 'A0 DBP': (sqrt(A_DBP - DBP / beta)^2
ICS.Pdiastolic		= REF.DBP;


%%	5. HISTORY POINTS
HIS.ArteryOutput	= [1 2 3 4 5 6 7];
% HIS.ArteryOutput	= [1];


%%	Generate .IN file
WriteInFile(PATHS,PARAM,MESH,BCS,ICS,HIS)


%%	Save data to .MAT file
save([PATHS.Aortic1D_files,PATHS.PatientName,'.mat']);


%%	END OF CODE		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end
