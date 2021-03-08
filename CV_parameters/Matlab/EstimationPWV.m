%%	Estimate PWV from reference data

function [EST] = EstimationPWV(PATHS,REF,EST)

switch EST.MethodPWV
		
	case {1, 2, 3, 4}	%TT algorithm
		addpath([PATHS.Functions,'PWV_TT'])
		
		SR = round(1/(REF.t_in(2) - REF.t_in(1)));

		switch EST.MethodPWV
			case {1, 3} %Ascending-descending aorta flow waveforms
				Wave_type = 2;	%flow
				
				%	Read waveform pair
				Wave_in = REF.Asc_Ao.ONE_CYCLE.Q;
				Wave_out = REF.Desc_Ao_I.ONE_CYCLE.Q;
				%	Ensure equal vector length
				[Wave_in, Wave_out] = CorrectSize(Wave_in, Wave_out);
				%	Increase number of cycles 
				Wave_in = [Wave_in; repmat(Wave_in(2:end),3,1)];
				Wave_out = [Wave_out; repmat(Wave_out(2:end),3,1)];
				%	Identify starting time index for each waveform
				Start_in = REF.Asc_Ao.ONE_CYCLE.start;
				Start_out = REF.Desc_Ao_I.ONE_CYCLE.start;
				if Start_in == Start_out, error('Both waves cannot start at same time') , end
				while Start_out < Start_in
					Start_out = Start_out + length(REF.Asc_Ao.ONE_CYCLE.Q) - 1;
				end
				%	Determine arterial path
				if ~isfield(REF,'TT_length')
					%	55/116-art locations [1 2 14 18]
					Aortic_length = sum([REF.Asc_Ao.length, REF.Ao_Arch_I.length, REF.Ao_Arch_II.length, REF.Desc_Ao_I.length]);	%[m]
					REF.TT_length = Aortic_length;
					warning('hardcoded TT length for foot-to-foot algorithm')
				end
				
			case {2, 4}	%Carotid-femoral pressure waveforms
				Wave_type = 1;	%pressure
				
				%	Read waveform pair
				Wave_in = REF.L_Carotid.ONE_CYCLE.P;
				Wave_out = REF.L_Femoral.ONE_CYCLE.P;
				%	Ensure equal vector length
				[Wave_in, Wave_out] = CorrectSize(Wave_in, Wave_out);
				%	Increase number of cycles 
				Wave_in = [Wave_in; repmat(Wave_in(2:end),3,1)];
				Wave_out = [Wave_out; repmat(Wave_out(2:end),3,1)];
				
				%	Identify starting time index for each waveform
				Start_in = REF.L_Carotid.ONE_CYCLE.start;
				Start_out = REF.L_Femoral.ONE_CYCLE.start;
				if Start_in == Start_out, error('Both waves cannot start at same time') , end
				while Start_out < Start_in
					Start_out = Start_out + length(REF.L_Carotid.ONE_CYCLE.P) - 1;
				end
				%	Determine arterial path
				if ~isfield(REF,'TT_length')
					%	116/55-art locations [14 18 27 28 35 37 39 41 43 50]
					Carotid_length = REF.L_Carotid.ONE_CYCLE.artery_dist;		%[m]
					Femoral_length = REF.L_Femoral.ONE_CYCLE.artery_dist;		%[m]
					Aortic_length = REF.Asc_Ao_L_Femoral.dist;	%[m]
					REF.TT_length = Aortic_length + Femoral_length - Carotid_length;
					warning('hardcoded TT length for foot-to-foot algorithm')
				end
		end
		
		Wave_in = Wave_in(Start_out - Start_in + 1 : end);
		[Wave_in, Wave_out] = CorrectSize(Wave_in, Wave_out);
		
		switch EST.MethodPWV
			case {1, 2}
				TT = TTAlgorithm([Wave_in';Wave_out'],SR,1,Wave_type,1,EST.Plots);
			case {3, 4}
				TT = TTAlgorithm([Wave_in';Wave_out'],SR,3,Wave_type,1,EST.Plots);
		end
		
		if EST.Plots == 1
% 			warning('Temporary: thesis configuration only')
% 			title('Foot-to-foot method')
% 			xlabel('Time [s]')
% 			ylabel('Pulse wave [AU]')
% 			xlim([0,2.5])
% 			ylim([-0.0001,0.0005])
% 			set(gca,'XTick',[0:0.5:2.5])
% 			set(gca,'YTick',[])
% 			box on
% 			set(gca,'FontSize',40)
% 			addpath([PATHS.Functions,'PlotSave/'])
% 			PATHS.Figures = '/Users/joh15/PhD/THESIS - PhD/v6/Latex/Figures/';
% 			PlotSave(PATHS.Figures,'f2f_method')
% 			PlotSave(PATHS.Figures,'least_squares_method')
		end
		
		EST.PWV = REF.TT_length/mean(TT);
		
	case 5	%sum of squares
		addpath([PATHS.Functions,'SumOfSquares/'])
		addpath([PATHS.Functions,'InterpolateSpline/'])
		
		P = REF.Pressure;				%[Pa]
		U = REF.Q_in/REF.Asc_Ao.A_in;	%[m/s]
		dens = REF.rho;					%[kg/m3]
		K = 5;
		F = 11;
		
		[P,U] = InterpolateSpline(REF.t_in,P,U);
		
		EST.PWV = c_SumSqrt_sgolay(P,U,dens,K,F);
		
end

end


%%	Function definition
function [Wave_in, Wave_out] = CorrectSize(Wave_in, Wave_out)

if length(Wave_in) > length(Wave_out)
	Wave_in = Wave_in(1:length(Wave_out));
elseif length(Wave_in) < length(Wave_out)
	Wave_out = Wave_out(1:length(Wave_in));
end

end
