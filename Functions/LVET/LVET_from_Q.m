%%	Extract LVET from flow and time data
%
%==========================================================================
%	Jorge Mariscal-Harana, King's College London (original author)
%	v1.0 (2018)
%
%==========================================================================

function [LVET] = LVET_from_Q(PATHS,Q_in,t_in,Plots)

Q_in = Q_in(:);
t_in = linspace(t_in(1),t_in(end),length(Q_in))';

[gl_max, gl_max_loc] = max(Q_in);	%find global max

[gl_min, gl_min_loc] = min(Q_in(gl_max_loc:round(end/2)));	%find global min after global max for the first half of the cycle
gl_min_loc = gl_max_loc + gl_min_loc - 1;

Q_in_dias = Q_in(gl_min_loc:end);		%Q for the diastolic portion
Q_in_temp = Q_in_dias(1:round(end/2));	%analyse the first half of the diastolic portion

Diastolic_threshold = 0.01*max(Q_in);

if sum(abs(Q_in_temp) > Diastolic_threshold) == 0	%all Q values after diastolic min are smaller than 1% of Q_max
	LVET = t_in(gl_min_loc);
	
else
	sign_change = find(sign(Q_in_temp(1:end-1)) == -1 & (sign(Q_in_temp(2:end)) == 1 | sign(Q_in_temp(2:end)) == 0));	%find neg -> pos after global min
	sign_change = sign_change + gl_min_loc;
	
	[~, max_peaks_loc] = findpeaks([Q_in_temp; -gl_max]);	%find first local max after global min
	max_peaks_loc = max_peaks_loc + gl_min_loc - 1;
	
	zero_locs = find(Q_in_temp == 0);	%find Q = 0 after global min	
	zero_locs = zero_locs + gl_min_loc - 1;
	
	if Plots == 1
		Plot_Q_after_min(t_in,Q_in,gl_max_loc,gl_max,gl_min_loc,gl_min,max_peaks_loc,sign_change)
		addpath([PATHS.Root,'Others/PlotSave/'])
		PlotSave(PATHS.Figures,'LVET_Q_analysis')
	end
	
	if ~isempty(sign_change)
		LVET_index = min(sign_change(1),max_peaks_loc(1));
	elseif ~isempty(zero_locs) && ~isempty(max_peaks_loc)
		LVET_index = min(zero_locs(1),max_peaks_loc(1));
	elseif ~isempty(zero_locs)
		LVET_index = zero_locs(1);
	elseif ~isempty(max_peaks_loc)
		LVET_index = max_peaks_loc(1);
	else
		LVET_index = round(0.37/sqrt(t_in(end))*length(t_in));
		warning('LVET automatically set as 0.37*sqrt(T); check your flow wave')
	end
	
	LVET = t_in(LVET_index);
end

end


%%	Function definitions
function Plot_Q_after_min(t,Q,gl_max_loc,gl_max,gl_min_loc,gl_min,max_peaks_loc,sign_change)
Scale_Q = 10^6;

figure
hold on
plot(t,Q*Scale_Q,'k','LineWidth',2)
plot(t(gl_max_loc),gl_max*Scale_Q,'^b','MarkerSize',15,'MarkerFaceColor','b')
plot(t(gl_min_loc),gl_min*Scale_Q,'vb','MarkerSize',15,'MarkerFaceColor','b')
plot(t(max_peaks_loc(1)),Q(max_peaks_loc(1))*Scale_Q,'^r','MarkerSize',15)
if ~isempty(sign_change)
	plot(t(sign_change(1)),Q(sign_change(1))*Scale_Q,'og','MarkerSize',10,'MarkerFaceColor','g')
end
hold off

if ~isempty(sign_change)
	legend('Q','global max','global min','local max','sign change')
else
	legend('Q','global max','global min','local max')
end

% title('t_{sys} (flow): 0-crossing/local max after min')
xlabel('Time [s]')
ylabel('Flow [mL/s]')
xlim([0,t(end)])
% ylim([-0.0001,0.0005])
set(gca,'XTick',[0:0.2:t(end)+0.2])
set(gca,'YTick',[0:100:400])
box on
legend boxoff
set(gca,'FontSize',40)

end