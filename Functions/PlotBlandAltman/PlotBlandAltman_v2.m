%%	Function to generate Bland Altman plots
%	Barry Greene, September 2008. Bland, J.M., Altman, D.G. 'Statistical 
%	methods for assessing agreement between two methods of clinical 
%	measurement'(1986) Lancet, 1 (8476), pp. 307-310.
%
%	INPUT
%	-REF:	reference data - single vector
%	-EST:	estimation data - can be a single vector or a cell of vectors
%	-PLOTS:	options for advanced plot formatting
%
%	OUTPUT
%	Bland-Altman plot with bias and limits of agreement
%
%==========================================================================
%	Jorge Mariscal-Harana, King's College London
%	v2.0 (30/01/19) clean version for multiple comparisons
%
%==========================================================================

function [md,sd] = PlotBlandAltman_v2(REF,EST,PLOTS)
%%	Ensure that REF is a row vector
REF = REF(:)';


%%	Determine if REF is compared against multiple EST and ensure that EST is a row vector
if iscell(EST)
	for jj = 1:length(EST)
		EST_temp(jj,:) = EST{jj};
	end
	EST = EST_temp;
else
	EST = EST(:)';
end

if(length(REF)~=length(EST)), error('REF and EST must be the same length'), end


%%	Calculate differences and SD from data
for jj = 1:size(EST,1)
% 	data_mean(jj,:) = mean([REF;EST(jj,:)]);		% Average of measurements from each instrument
	data_mean(jj,:) = REF;							% Reference measurements
	data_diff(jj,:) = EST(jj,:) - REF;				% Difference between data from each instrument
	md(jj) = mean(data_diff(jj,:));	% Mean of difference between instruments (bias)
	sd(jj) = std(data_diff(jj,:));	% SD of differewnce between instruments
end


%%	Bland-Altman plots
%	Plot bias and limits of agreement for (first estimation) - (reference)
figure, hold on
plot([0 10*max(data_mean(jj,:))],md(1)*ones(1,2),'k-','LineWidth',PLOTS.LineWidth)					% Mean difference line
plot([0 10*max(data_mean(jj,:))],(md(1) + 1.96*sd(1))*ones(1,2),'k--','LineWidth',PLOTS.LineWidth)	% Mean + 1.96*SD line
plot([0 10*max(data_mean(jj,:))],(md(1) - 1.96*sd(1))*ones(1,2),'k--','LineWidth',PLOTS.LineWidth)	% Mean - 1.96*SD line
for jj = 1:size(EST,1)
	plot(data_mean(jj,:),data_diff(jj,:),PLOTS.Markers{jj},'MarkerSize',PLOTS.MarkerSize,'LineWidth',PLOTS.LineWidth);	% Bland Altman plot
end
hold off


%%	Format plot
if PLOTS.Format == 1
	PLOTS.md = md(1);
	PLOTS.sd = sd(1);
	
	addpath('/Users/joh15/Haemodynamic_Tools/Version6/Others/PlotFormat/')
	FormatPlots(PLOTS)
end

end

