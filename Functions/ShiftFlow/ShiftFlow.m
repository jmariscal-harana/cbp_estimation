function [Wave_shift2,shift_factor] = ShiftFlow(Wave)

Wave = Wave(:);

%	Analyse Q from start to global max
[~, max_loc] = max(Wave);
Wave_temp = Wave(1:max_loc);

%	Find locations where Q < 0 and make them = 0
Wave_negative = find(Wave_temp < 0);
if ~isempty(Wave_negative)
	Wave_negative = Wave_negative(end);
end

%	Shift waveform so that cycle starts at last early systole 0 value
if Wave_negative == 1
	Wave_shift1 = [0; Wave(2:end)];	%first datapoint is negative; make it = 0
elseif Wave_negative > 1	
	Wave_temp = Wave;
% 	Wave_shift1 = [0; Wave_temp(Wave_negative+1:end); zeros(Wave_negative-1,1)];	%pad with zeros
	Wave_shift1 = [0; Wave_temp(Wave_negative+1:end); Wave_temp(end)*ones(Wave_negative-1,1)];	%pad with last value
else
	Wave_shift1 = Wave;
end

%	Shift waveform according to early-systole max(dQ) and x-intercept
Wave_temp_dQ = diff(Wave_shift1);
[slope, ind_max] = max(Wave_temp_dQ(1:ceil(length(Wave_temp_dQ)/2)));	%find largest positive slope for 50% of the cycle
temp = Wave_shift1(ind_max) + slope*[-1000:0]; if sum(temp < 0) == 0, error('x-intercept required'), end
temp = find(temp>=0);
shift_factor = ind_max - length(temp);
if shift_factor < 0 
	error('Shift factor < 0')	%not physically possible
elseif shift_factor <= 1
	Wave_shift2 = Wave_shift1;	%no shifting required
else 
% 	Wave_shift2 = [0; Wave_shift1(shift_factor+1:end); zeros(shift_factor-1,1)];	%pad with zeros
	Wave_shift2 = [0; Wave_shift1(shift_factor+1:end); Wave_shift1(end)*ones(shift_factor-1,1)];	%pad with last value
end

% figure, hold on,
% plot(Wave,'--b')
% plot(Wave_shift1,'bo')
% % plot([-1000+ind_max:ind_max], [temp], 'r');
% plot(Wave_shift2,'rx')
% ylabel('VTI [m/s]')
% xlabel('Datapoint ID')
% legend('Cubic spline interpolation','Shift: remove Q_{upstroke} < 0', 'Shift: max(dQ) backprojection')
% set(gca,'fontsize',20)

end