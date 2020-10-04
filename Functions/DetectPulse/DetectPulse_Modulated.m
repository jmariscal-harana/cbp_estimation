%%	Extract single cycle from multiple cycles of modulated data
%	Input
%	-Waveform
%
%	Output
%	-Start: cycle start index 
%	-Duration: average cycle duration (no units)
% 
%==========================================================================
%	Jorge Mariscal-Harana, King's College London (original author)
%	v1.0 (07/05/18)
%
%==========================================================================
function [Start,Duration] = DetectPulse_Modulated(Waveform)
N = length(Waveform);
kk = 1;

Start = [];

for jj = 1:N-1
	if Waveform(jj) == 0 && Waveform(jj+1) ~= 0
		Start(kk) = jj;
		kk = kk + 1;
	end
end

if isempty(Start)
	error('I could not find the start of a cycle (where Waveform(n) = 0). Check the input waveform data')
end

Duration = mean(Start(2:end) - Start(1:end-1));

%	Plot first cycle
plot(Waveform([Start(1):Start(2)]))
pause(0.1)

end



	