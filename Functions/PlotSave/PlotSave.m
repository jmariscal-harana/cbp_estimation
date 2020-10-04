function PlotSave(Path,Name,Number)

% saveas(gcf,[Path,Name],'epsc')
% saveas(gcf,[Path,Name],'png')
% savefig([Path,Name])

if nargin == 2
	set(gcf,'PaperOrientation','landscape');
	set(gcf,'PaperUnits','normalized');
	set(gcf,'PaperPosition', [0 0 1 1]);
	print(gcf, '-dpdf', [Path,Name]);
elseif nargin == 3
	set(Number,'PaperOrientation','landscape');
	set(Number,'PaperUnits','normalized');
	set(Number,'PaperPosition', [0 0 1 1]);
	print(Number, '-dpdf', [Path,Name]);
else
	error('Incorrect number of function inputs')
end

end