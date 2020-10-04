%%	Function to generate Bland Altman plots
%	Barry Greene, September 2008. Bland, J.M., Altman, D.G. 'Statistical 
%	methods for assessing agreement between two methods of clinical 
%	measurement'(1986) Lancet, 1 (8476), pp. 307-310.
%
%	INPUT
%	-Pressure_REF: gold standard data
%	-Pressure_EST{j}: estimation data 
%
%	OUTPUT
%	Produces Bland Altman plot with mean difference and mean difference +/-
%	2*SD difference lines.

%==========================================================================
%	Jorge Mariscal-Harana, King's College London
%	v1.0 (14/06/17) solved array size issue; improved plot format
%	v2.0 (26/06/17) xlim improved; title on/off added
%	v2.1 (27/06/17) xlim can now be an input
%	v3.0 (28/09/17) works for multiple datasets; PLOTS parameters added
%	v4.0 (15/11/18) generalised for any dataset comparison
%
%==========================================================================

function [md,sd] = PlotBlandAltman(REF,EST,PLOTS)
%%	Calculate differences and SD from data
[m,n] = size(REF); if(n>m), REF = REF'; end
[m,n] = size(EST); if(n>m), EST = EST'; end
if(size(REF)~=size(EST)), error('Data matrices must be the same size'), end
	
data_mean = mean([REF,EST],2);		% Average of measurements from each instrument
data_diff = EST - REF;				% Difference between data from each instrument
md = round(mean(data_diff),2);	% Mean of difference between instruments (bias)
sd = round(std(data_diff),2);	% SD of difference between instruments


%%	Bland-Altman plots
%	Plot bias and limits of agreement (LoA) for: (first estimation) - (reference)
figure
hold on
plot([0 10*max(data_mean)],md*ones(1,2),'k-','LineWidth',PLOTS.LineWidthM)				% Mean difference line  
plot([0 10*max(data_mean)],(md+1.96 *sd(1))*ones(1,2),'k--','LineWidth',PLOTS.LineWidthM)	% Mean plus 2*SD line  
plot([0 10*max(data_mean)],(md-1.96 *sd(1))*ones(1,2),'k--','LineWidth',PLOTS.LineWidthM)	% Mean minus 2*SD line   
plot(data_mean,data_diff,PLOTS.Markers,'MarkerSize',PLOTS.FontSizeL,'LineWidth',PLOTS.LineWidthS);	% Bland Altman plot
hold off


%%	Format plot
if PLOTS.Format == 1
	addpath('/Users/joh15/Haemodynamic_Tools/Version6/Others/PlotFormat/')
	FormatPlots(PLOTS)
end

end

