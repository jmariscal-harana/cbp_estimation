%%	Shaded plot showing max/min P, Q, etc. values for multiple datasets
%	Input
%	-Time: x-axis, time vectors
%	-Waves: extract max/min from these waveforms
%
%	Output
%	-Shaded plot
%
%==========================================================================
%	Jorge Mariscal-Harana, King's College London (original author)
%	v1.0 (10/05/18)
%	v1.1 (19/04/19) - improved zero value handling
%
%==========================================================================

function [t_max] = ShadedPlots(Time,Waves,NaN_on)
%%	Rearrange data to extract max time, and min/max waveform
N_waves = length(Waves);

for jj = 1:N_waves
	dt(jj) = Time{jj}(2) - Time{jj}(1);
end

dt_fine = min(dt);

for jj = 1:N_waves
	Time_fine{jj} = Time{jj}(1):dt_fine:Time{jj}(end);
	Waves_spline{jj} = spline(Time{jj},Waves{jj},Time_fine{jj});
end


for jj = 1:N_waves
	Wave_length(jj) = length(Waves_spline{jj});
end

Wave_length_max		= max(Wave_length);
Wave_length_min		= min(Wave_length);
Wave_zeros			= zeros(N_waves,Wave_length_max);

for jj = 1:N_waves
	dummy				= Waves_spline{jj};
	N					= Wave_length(jj)+1;
	if N <= Wave_length_max
		dummy(N:Wave_length_max) = 0;
	end
	Wave_zeros(jj,:)	= dummy;
end

DATA = Wave_zeros;

%	Remove zeros arising from autofilling
if NaN_on == 1
	for jj = 1:N_waves
		Detect_zeros = (DATA(jj,:) == 0);
		Detect_zeros(1:Wave_length(jj)) = 0;
		DATA(jj,Detect_zeros) = NaN;
	end

end

Upper	= max(DATA);
Lower	= min(DATA);

clear DATA

[~, t_index]	= max(Wave_length);
t_max			= Time_fine{t_index};

if size(t_max,1) < size(t_max,2)
else
	t_max = t_max';
end

jbfill(t_max,Upper,Lower,[0.5 0.5 0.5],[0.5 0.5 0.5],0,0.25);

end