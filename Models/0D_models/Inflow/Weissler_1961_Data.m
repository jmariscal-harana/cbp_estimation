%%	Ejection time, heart rate and stroke volume data for 17 normal subjects in supine position from Weissler 1961 (10.1016/0002-8703(61)90403-3)
ET_normal = [.27, .275, .235, .265, .29, .2, .315, .32, .31, .285, .305, .26, .3, .27, .305, .285, .25];	% Ejection time (s)
HR_normal = [68, 68, 85, 72, 49, 120, 56, 60, 66, 65, 56, 72, 62, 68, 57, 65, 80];	% Heart rate (bpm)
SV_normal = [82, 81, 62, 95, 98, 44, 106, 102, 96, 82, 92, 74, 107, 93, 101, 109, 79];	% Stroke volume for 17 normal subjects (mL)

ET_from_HR = .266 - .0021*(HR_normal - 73);	% Regression equations extracted from data
ET_from_SV = .266 + .0017*(SV_normal - 82);	% Regression equations extracted from data
ET_from_HR_SV = .266 + .0011*(SV_normal - 82) - .0009*(HR_normal - 73);	% Regression equations extracted from data

%	Comparison betweeen regression equations for ejection time
figure
hold on
plot(ET_normal, ET_from_HR, 'o', 'LineWidth', 2)
plot(ET_normal, ET_from_SV, '*', 'LineWidth', 2)
plot(ET_normal, ET_from_HR_SV, 'd', 'LineWidth', 2)
plot([0 1], [0 1], 'k--')
hold off
xlim([0.15 0.35])
ylim([0.15 0.35])
xlabel('ET (s)')
ylabel('ET (s)')
legend('ET = .266 - .0021(HR - 73)', 'ET = .266 + .0017(SV - 82)', ...
	char({'ET = .266 + .0011(SV - 82)','                 - .0009(HR - 73)'}), 'Location', 'NorthWest')
set(gca, 'FontSize', 18')

% saveas(gcf, [PATHS.Figures, 'ET_regression_HR_SV'], 'pdf')
% saveas(gcf, [PATHS.Figures, 'ET_regression_HR_SV'], 'png')
