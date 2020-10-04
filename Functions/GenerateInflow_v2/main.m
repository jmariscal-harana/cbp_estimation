%%	Creates a harmonic .BCS file from a time-flow file

InVivo.Data     = '1D-inflow-file.txt';	%Time-flow file
InVivo.Inflow   = '1D-input-file';		%Harmonics file

InVivo.SR               = 1000;		% Sampling rate of the interpolated waveforms
InVivo.Q.NHarm          = 20;		% Number of harmonics used to define the flow waveforms
InVivo.A.NHarm          = 20;		% Number of harmonics used to define the area waveforms
InVivo.P.NHarm          = 20;		% Number of harmonics used to define the pressure waveforms
InVivo.ShowHarmonics    = 0;		% 0: Do not show harmonics info;
									% 1: Show harmonics info

InVivo.Tcorrection      = 0;        % Calculate the cardiac period as
									% 0: given by the inflow waveform
									% 1: the mean of the periods of all flow waveforms

InVivo.Q.Root           = 100;      % Inflow waveform number
InVivo.Plots			= 1;		% 0: do not show in vivo plots
									% 1: show all in vivo plots
%	Scaling (improve this)
InVivo.SCAL_Q			= 1e-6;		
InVivo.SCAL_C			= 1e-6;		% Convert compliance into m^3/Pa
InVivo.SCAL_R			= 1e6;		% Convert resistance into Pa/(m^3/s)
InVivo.SCAL_A			= 1e-4;		% Convert measured area into m^2
InVivo.SCAL_t			= 1e-3;		% Convert time of measured flow rate into s

InVivo.Q.tpeak			= 0.125;

%%	Generate the .bcs inflow file
InVivo = oneDbio_GenerateInflowFile(InVivo);