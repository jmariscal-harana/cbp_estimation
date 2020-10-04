function [] = PlotLinearCorrelation(REF,EST,PLOTS)

%%	Scatter plots: systolic, mean, diastolic, PP, etc.
figure

hold on
for jj = 1:length(EST)
	plot(REF,EST{jj},PLOTS.Markers{jj},'MarkerSize',PLOTS.MarkerSize,'LineWidth',PLOTS.LineWidth)
end

%	PLOT IDENTITY LINE
plot([-500,500],[-500,500],'k--','LineWidth',PLOTS.LineWidth)
hold off


%%	Format plot
if PLOTS.Format == 1
	addpath('/Users/joh15/Haemodynamic_Tools/Version6/Others/PlotFormat/')
	FormatPlots(PLOTS)
end

end