%%	General parameters of the simulation
%	Input
%	-REF: reference haemodynamic data
%	-PARAM: required to specify or calculate some parameters
%
%	Output
%	-PARAM: parameter list values and data types
%
%==========================================================================
%	Jorge Mariscal-Harana, King's College London (original author)
%	v1.0 (16/10/18)
%	v1.1 (14/03/19) - took out hardcoded parameter values
%==========================================================================

function [PARAM] = ParameterList(REF,PARAM,UP)
%%	For any simulation, the following parameters must always be included:
PARAM.EQTYPE		= {PARAM.EQTYPE, '%d'};						%Linear/nonlinear formulation
PARAM.DT			= {PARAM.DT, '%.3e'};						%Time step [s]
PARAM.NSTEPS		= {UP.Duration/PARAM.DT{1}, '%.3e'};		%Number of time steps for the simulation 
PARAM.HISSTEP		= {round(UP.His_dt/PARAM.DT{1}), '%d'};	%Number of time steps for the .his file
PARAM.Rho			= {PARAM.Rho, '%.1f'};						%Blood density [Kg/m3]
PARAM.Viscosity		= {PARAM.Viscosity, '%.4f'};				%Blood viscosity [Pa s]


%%	The following parameters are optional (include only those to be printed to the .IN file):
PARAM.INTTYPE		= {PARAM.INTTYPE,	'%d'};								%Integration order of numerical scheme: 1, 2 (DEFAULT) or 3
% PARAM.IOSTEP		= {0,	'%.3e'};							%Number of time steps until the next step solution is dumped in the .out file. (DEFAULT = 0)
% PARAM.Beta			= {0,	'%.3e'};	
% PARAM.Gamma			= {0,	'%.3e'};
% PARAM.Ao			= {0,	'%.3e'};
PARAM.pinf			= {REF.P_out,	'%.3f'};					%Pressure at 3-element Windkessel outflow (DEFAULT = 0) [Pa]
PARAM.Pext			= {REF.DBP,	'%.3f'};						%Pressure to the vessel wall (DEFAULT = 0) [Pa]
% PARAM.Junc_losses	= {0,	'%d'};
% PARAM.Bifurc_losses	= {0,	'%d'};
% PARAM.Periodic		= {0,	'%.5f'};
PARAM.T_initial		= {PARAM.T_initial,	'%.5f'};
PARAM.T_final		= {PARAM.T_final,	'%.5f'};
% PARAM.SCAL_F		= {0,	'%.4f'};							%Scaling factor for the inflow waveform
PARAM.Alpha			= {(2 + UP.mcf)/(1 + UP.mcf), '%.3f'};		%Velocity profile parameter
PARAM.Alpha			= {1, '%.3f'};								%1 for inviscid flow
PARAM.Alpha			= {4/3, '%.3f'};							%4/3 for Poiseuille flow - Pete
% PARAM.Factor		= {0,		'%d'};		% Scaling factor for sensitivity analysis
% 											% 0: not in use
% PARAM.WallViscosity	= {0,		'%.4f'};		% Wall viscosity (Pa s) - 0: not used
% PARAM.hD			= {0.07,	'%.4f'};		% Thickness to diameter ratio
% PARAM.PI			= {pi,	'%2.3e'};		% value of pi

end