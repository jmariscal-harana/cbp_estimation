%%	Geometrical and mechanical properties of each arterial segment (called a "domain" in Nektar1D)
%	Input
%	-NUM: required to specify or calculate some parameters
%
%	Output
%	-PARAM: parameter list values and data types
%
%==========================================================================
%	Jorge Mariscal-Harana, King's College London (original author)
%	v1.0 (16/10/18)
%
%==========================================================================

function [MESH] = MeshDefinition(REF,PARAM,MESH)
%%	Arterial geometry
%	Copy from REF to MESH
MESH.L				= REF.L;
MESH.Rad_prox		= REF.Rad_prox;
MESH.Rad_dist		= REF.Rad_dist;

MESH.nbSeg			= length(MESH.L);
MESH.nbElemBySeg	= ceil(MESH.L./MESH.elmt);		%number of elements in each segment

%	Polynomial order for normal segments
MESH.p = MESH.p_order.*ones(MESH.nbSeg,1);		%vector with polynomial orders
MESH.q = MESH.q_order.*ones(MESH.nbSeg,1);		%vector with quadrature orders

%	Polynomial order for short elements
MESH.elmt_shortL_reduced = 0.011;	%Reduces complexity if element is less than this length [m] (Marie's)
% MESH.elmt_shortL_reduced = 0.019;	%Reduces complexity if element is less than this length [m] (Pete's)
MESH.elmt_shortL_vw = 0.011; %Original setting [m] (Marie's)
% MESH.elmt_shortL_vw = 0.009; %Ensures that all the arteries in the path of the digital artery are visco-elastic [m] (Pete's)

MESH.p_order_reduced = 2;					%Reduced order for short segments
MESH.q_order_reduced = 2;					%Reduced order for short segments

MESH = ReducedElements(MESH);

% Radius(x) = (Rad_dist - Rad_prox)*x/L + Rad_prox;
%           = params.mesh.Rad_coef(:,1)*x + params.mesh.Rad_coef(:,2);
% Area(x) = pi*Radius(x)*Radius(x) 
MESH.Rad_coef(:,1) = (MESH.Rad_dist - MESH.Rad_prox)./MESH.L;
MESH.Rad_coef(:,2) = MESH.Rad_prox;


%%	Mesh connectivity
MESH.NodesAtEdges	= REF.Nodes;
MESH.node1			= MESH.NodesAtEdges(:,1);
MESH.node2			= MESH.NodesAtEdges(:,2);
MESH.E				= MESH.nbSeg;								%Number of edges
MESH.N				= max([max(MESH.node1) max(MESH.node2)]);	%Number of nodes

if MESH.E ~= MESH.N - 1, error('Number of Segments ~= Number of Nodes - 1'), end

MESH				= MeshConnectivity(MESH);


%%	Mesh wall properties
switch MESH.MechProps
	case 'Beta'
		MESH.beta_coef(1:MESH.nbSeg) = 2*PARAM.Rho{1}*REF.PWV^2;
	case 'Eh'
		%	Empirical parameters for arterial stiffness
		MESH.Eh_k1 = REF.k(1);
		MESH.Eh_k2 = REF.k(2);
		MESH.Eh_k3 = REF.k(3);
		MESH.gamma_b0 = REF.gamma_b0;
		MESH.gamma_b1 = REF.gamma_b1;
	case 'Empirical'
		MESH.PWV_a(1:MESH.nbSeg) = 14.3;	%only required for 'Empirical' wall properties
end

switch MESH.PWVType
	case 'f2f'
		MESH.c_prox = REF.PWV.*ones(length(MESH.L),1);
		MESH.c_dist = REF.PWV.*ones(length(MESH.L),1);
	case 'Prescribed'
		MESH.c_prox = REF.c_prox;
		MESH.c_dist = REF.c_dist;
end


end

%%	Function definitions
%	If the length of the element is < MESH.elmt_shortL_reduced, 
% - group the two last elements of the domain (if there are more than 1 element in the domain)
% - or reduce the polynomial and quadrature order to 2 (if there is only one element in the domain)
function [MESH] = ReducedElements(MESH)
displayed_red_order = 0;
% cycle through each element
for i=1:MESH.nbSeg
	
	% identify distances of each element. NB: this elaborate scheme avoids the
	% rounding errors sometimes produced by:      temp = [[0:UP.elmt:MESH.L(i)],MESH.L(i)];
	no_els = ceil(MESH.L(i)/MESH.elmt);
	no_complete_els = floor(MESH.L(i)/MESH.elmt);
	temp = 0;
	for el_no = 1 : no_complete_els
		temp = [temp, el_no*MESH.elmt];
	end
	if no_els>no_complete_els
		temp(end+1) = MESH.L(i);
	end
	
	%ensure you don't repeat the last point twice (if L=multiple of MESH.elmt)
	if (temp(end) == temp(end-1))
		temp=temp(1:end-1);
	end
	
	% calculate shortened length of the last element of this segment
	shortL = temp(end)-temp(end-1);
	
	% if the length of this last element is less than a threshold value, then
	if (shortL <MESH.elmt_shortL_reduced)
		% if there is more than one element in this segment
		if(not(MESH.nbElemBySeg(i) == 1))
			% then join the second last and last segment into one
			temp = [temp(1:end-2), temp(end)];
			MESH.nbElemBySeg(i) = MESH.nbElemBySeg(i)-1;
		else
			% use reduced order for this segment
			MESH.p(i) = MESH.p_order_reduced;
			MESH.q(i) = MESH.q_order_reduced;
			if ~displayed_red_order
				fprintf(['NB: Reduced order at element:\n',num2str(i),]);
				displayed_red_order = 1;
			else
				fprintf([num2str(i),'\n'])
			end
		end
	end
	clear shortL
	
	% store the distances the elements in this segment
	MESH.x_elem{i} = temp; clear temp
end

end

function [MESH] = MeshConnectivity(MESH)
%	MESH contains the variables
%   -MESH.N : number of nodes in network
%   -MESH.E : number of edges in network
%   -MESH.NodesAtEdges (E x 2): Start and end nodes for each edge

%	The function constructs the following geometry matrices:
%   -MESH.B            (NxE)  : Incidence matrix
%   -MESH.EdgesByNode  (N x 3): All edges connected to each node
%   -MESH.degree       (1xN)  : Degree of each node
%   -MESH.BifType      (1xN)  : Type of bifurcation for each node

%--- Incidence Matrix B
MESH.B = zeros(MESH.N, MESH.E);
MESH.BifType = zeros(MESH.N, 1);

for i = 1:MESH.E
    %start of node
    MESH.B(MESH.NodesAtEdges(i,1),i) = -1;
    %end of node
    MESH.B(MESH.NodesAtEdges(i,2),i) = 1;
end
    

%-- EdgesByNodes & degree of nodes
MESH.EdgesByNode(1:MESH.N,1:3) = 0;

for i=1:MESH.N
    % Get the degree of node i
    MESH.degree(i) = sum(abs(MESH.B(i,:)));
    if (MESH.degree(i)==3) 
        MESH.BifType(i) = sum(MESH.B(i,:));
    end
    % Get the edges connected to the internal node i
    clear e Y;
    s=1; 
    for j=1:MESH.E  
        if abs(MESH.B(i,j))>0
            e(s)=j;
            s=s+1;
        end
    end
    MESH.EdgesByNode(i,:) = [e zeros(1,4-s)];
end

k = 1;
for i = 1:MESH.E
	if sum(MESH.B(i+1,:)) == 1
		MESH.E_terminal(k) = i;
		k = k + 1;
	end
end

end

