%%	Generates the incidence matrix using the following inputs
% 	-Number of domains
% 	-Bifurcation nodes
% 	-Outlet nodes

%==========================================================================
%	Jorge Mariscal-Harana, King's College London
%	v1.0 (05/16)
%==========================================================================

function [MeshProp] = IncidenceMatrix(MeshProp)

Domains			= MeshProp.Art_Domains;
Nodes			= Domains + 1;
Bifurcations	= MeshProp.N_Bifurcations;
Outlets			= MeshProp.N_Outlets;

Incidence		= zeros(Nodes,Domains);

k = 1;

for i = 1:Nodes
	for j = 1:Domains
		if i == j && sum(i == Outlets) == 0
			Incidence(i,j) = -1;
		end
		if sum(i == Bifurcations) == 1 && j == Outlets(k)
			Incidence(i,j) = -1;
		end
		if i == j+1
			Incidence(i,j) = 1;	% this is always the case
		end
	end
	if k < length(Outlets) && i == Bifurcations(k)
		k = k+1;
	end
end

MeshProp.Incidence = Incidence;

end