%%	Extract REF data from in silico database as if it were in vivo data

function [REF] = Aortic1D_ExtractREF_v2(PATHS,REF,Dataset,Scenario)

warning('All inputs should be SI units')


%%	Extract geometry
if sum(strcmp(Dataset,{'pwdb_4374_reference';'pwdb_4374_reference_v2';'pwdb_6_reference'})) == 1
	REF								= ExtractGeometry_InSilico(REF);
elseif strcmp(Dataset,'Isra_AoCoarctation_reference')
	REF								= ExtractGeometry_InVivo(PATHS,REF);
end


%%	Extract waveform data required for parameter estimation
REF.t_in						= REF.Asc_Ao.ONE_CYCLE.t;
REF.Q_in						= REF.Asc_Ao.ONE_CYCLE.Q;

dt								= REF.t_in(2) - REF.t_in(1);
REF.SV							= REF.HAEMO.SV;

if sum(strcmp(Dataset,{'pwdb_4374_reference';'pwdb_4374_reference_v2';'pwdb_6_reference'})) == 1
	switch Scenario
		case 'Waveform'
			REF.Pressure					= REF.L_Carotid.ONE_CYCLE.P;
		case 'BP'
			REF.Pressure					= REF.L_Brachial.ONE_CYCLE.P;
	end
elseif strcmp(Dataset,'Isra_AoCoarctation_reference')
	REF.Pressure					= REF.Desc_Ao_I.ONE_CYCLE.P;
end

REF.DBP							= min(REF.Pressure);
REF.SBP							= max(REF.Pressure);

REF								= ExtractReferenceBP(Scenario,REF);
switch Scenario
	case 'BP'
		REF		= rmfield(REF,'Pressure');
end


%%	Additional parameters
REF.rho = 1060;	%[kg/m3]

		
end

function [REF] = ExtractGeometry_InSilico(REF)

%	Arterial segments of the 116-artery model
REF.Asc_Ao.Segment				= 1;
REF.Ao_Arch_I.Segment			= 2;
REF.Brachio.Segment				= 3;
REF.Ao_Arch_II.Segment			= 14;
REF.L_Carotid.Segment			= 15;
REF.Desc_Ao_I.Segment			= 18;
REF.L_Subclavian.Segment		= 19;
REF.L_Brachial.Segment			= 21;

Art_Domains_116_art = [REF.Asc_Ao.Segment
	REF.Ao_Arch_I.Segment
	REF.Ao_Arch_II.Segment
	REF.Desc_Ao_I.Segment
	REF.Brachio.Segment
	REF.L_Carotid.Segment
	REF.L_Subclavian.Segment]';

%	Arterial segments of the aortic model (must be in the same order as the 116-artery model segments)
Art_Domains_Aortic1D	= [1 2 3 4 5 6 7];
Art_Domains_Aorta		= [1 2 3 4];

%	Mesh connectivity for given aortic model segments:
%	Aortic root		-> Thoracic aorta
%	1st bifurcation	-> Brachiocephalic
%	2nd bifurcation	-> Carotid
%	3rd bifurcation	-> Subclavian
REF.Nodes	= [1 2;
	2 3;
	3 4;
	4 5;
	2 6;
	3 7;
	4 8];	%[inlet, outlet]

%	This should follow the same order as the mesh connectivity and geometry
Artery_name			= REF.INFO.Artery_geo_name;

for jj = 1:length(Artery_name)
	eval(['Length(jj)	= REF.',Artery_name{jj},'.length;']);
	eval(['R_in(jj)		= REF.',Artery_name{jj},'.R_in;']);
	eval(['R_out(jj)	= REF.',Artery_name{jj},'.R_out;']);
end


%%	Mesh dimensions
REF.L			= Length;	%[m]
REF.Rad_prox	= R_in;		%diastolic radius; [m]
REF.Rad_dist	= R_out;	%diastolic radius; [m]


%%	Aortic length for PWV calculation
REF.Aortic_L	= sum(REF.L(Art_Domains_Aorta));


%%	Visualise geometry
% VisualiseGeometry(REF)
% VisualiseGeometry_v2(REF)


end

function [REF] = ExtractGeometry_InVivo(PATHS,REF)
%	Connectivity matrix
Matrix_File = [PATHS.REF_data,REF.INFO.ID,'/Geometry_Arna/matrix.txt'];
Matrix = dlmread(Matrix_File);

Nodes_N = max(size(Matrix));
Segments_N = min(size(Matrix));

%	Visualise geometry
% VisualiseConnectivity(Matrix,Nodes_N,Segments_N)

%	Inlet and outlet nodes - All segments
for jj = 1:Segments_N
	Node_Inlet(jj) = find(Matrix(:,jj) == -1);
	Node_Outlet(jj) = find(Matrix(:,jj) == 1);
end

REF.Nodes					= [Node_Inlet(:), Node_Outlet(:)];	%[inlet, outlet]

%	Segments (domains) - aorta and supra-aortic arteries
Temp_1						= find(Node_Inlet ~= [1:Segments_N]);	%first segment for each supra-aortic artery
Segments_Aorta_N			= Temp_1(1) - 1;
Bifurcation_Nodes			= Node_Inlet(Temp_1);
[~,Bifurcation_Order]		= ismember(Bifurcation_Nodes,unique(Bifurcation_Nodes));

Temp_2						= find(Node_Outlet(1:end-1) ~= Node_Inlet(2:end));
Terminal_Nodes				= [Temp_2 + 1, Node_Outlet(end)];

if sum(Bifurcation_Nodes ~= sort(Bifurcation_Nodes)) ~= 0, warning('Order of segmentation is different'), end

Art_Segments_Aorta			= [1:Segments_Aorta_N];					%required to calculate aortic length
Art_Segments_SupraAortic	= [Segments_Aorta_N + 1:Segments_N];	%not required atm

%	Aorta nodes (inlet, bifurcations, outlet)
Aorta_Nodes					= Node_Inlet(Node_Inlet ~= [1:length(Node_Inlet)]);						%internal aortic nodes
Aorta_Nodes					= [Art_Segments_Aorta(1), Aorta_Nodes, Art_Segments_Aorta(end) + 1];	%all aortic nodes
Aorta_Nodes					= sort(Aorta_Nodes);

%	Segments - aorta
REF.Aortic_1.Segments		= [Aorta_Nodes(1):Aorta_Nodes(2)-1];
REF.Aortic_2.Segments		= [Aorta_Nodes(2):Aorta_Nodes(3)-1];
REF.Aortic_3.Segments		= [Aorta_Nodes(3):Aorta_Nodes(4)-1];
REF.Aortic_4.Segments		= [Aorta_Nodes(4):Aorta_Nodes(5)-1];

%	Segments - supra-aortic arteries
Supra{1}					= [Terminal_Nodes(1):Terminal_Nodes(2) - 1];
Supra{2}					= [Terminal_Nodes(2):Terminal_Nodes(3) - 1];
Supra{3}					= [Terminal_Nodes(3):Terminal_Nodes(4) - 1];
REF.Brachio.Segments		= Supra{Bifurcation_Order(1)};
REF.L_Carotid.Segments		= Supra{Bifurcation_Order(2)};
REF.L_Subclavian.Segments	= Supra{Bifurcation_Order(3)};

%	Aortic geometry file
Geometry_File = [PATHS.REF_data,REF.INFO.ID,'/Geometry_Arna/Area.txt'];

Geometry = ReadGeometryFile(Geometry_File);
Geometry = [Geometry{:}];

Data_NaN = find(isnan(Geometry));
Data_Num = find(~isnan(Geometry));

Number_of_Domains	= Geometry(Data_NaN(1)+1);
Length				= Geometry(Data_NaN(2)+1 : Data_NaN(3)-1);
Proximal_Area		= Geometry(Data_NaN(3)+1 : Data_NaN(4)-1);
Distal_Area			= Geometry(Data_NaN(4)+1 : Data_NaN(5)-1);
REF.PWV				= Geometry(Data_NaN(5)+1 : Data_NaN(6)-1); warning('PWV from Arna f2f extraction; please make automatic')	%TEMPORARY
% SBP					= Geometry(Data_NaN(6)+1 : Data_NaN(7)-1);
% DBP					= Geometry(Data_NaN(7)+1);

if Number_of_Domains ~= Segments_N
	error('N_domains in Matrix and Geometry files differ')
end

%	Save geometry information as REF
REF.L			= Length*1e-3;	%[m]
REF.Rad_prox	= sqrt(Proximal_Area/pi)*1e-3;	%mean radius; [m]
REF.Rad_dist	= sqrt(Distal_Area/pi)*1e-3;	%mean radius; [m]
REF.c_prox		= REF.PWV;					%PWV at mean radius; [m/s]
REF.c_dist		= REF.PWV;					%PWV at mean radius; [m/s]

%	Aortic length for f2f PWV calculation between aortic root and descending aorta
REF.TT_length	= sum(REF.L(Art_Segments_Aorta));	%[m]

%	Visualise geometry
% VisualiseGeometry(REF)

%	Combine segments belonging to a single artery (inlet node to outlet node)
Dom_1 = REF.Aortic_1.Segments;
Dom_2 = REF.Aortic_2.Segments;
Dom_3 = REF.Aortic_3.Segments;
Dom_4 = REF.Aortic_4.Segments;
Dom_5 = REF.Brachio.Segments;
Dom_6 = REF.L_Carotid.Segments;
Dom_7 = REF.L_Subclavian.Segments;

Length = [sum(Length(Dom_1)), sum(Length(Dom_2)), sum(Length(Dom_3)), ...
	sum(Length(Dom_4)), sum(Length(Dom_5)), sum(Length(Dom_6)), sum(Length(Dom_7))];
Proximal_Area = [Proximal_Area(Dom_1(1)), Proximal_Area(Dom_2(1)), Proximal_Area(Dom_3(1)), ...
	Proximal_Area(Dom_4(1)), Proximal_Area(Dom_5(1)), Proximal_Area(Dom_6(1)), Proximal_Area(Dom_7(1))];
Distal_Area = [Distal_Area(Dom_1(end)), Distal_Area(Dom_2(end)), Distal_Area(Dom_3(end)), ...
	Distal_Area(Dom_4(end)), Distal_Area(Dom_5(end)), Distal_Area(Dom_6(end)), Distal_Area(Dom_7(end))];

%	Save geometry information as REF
REF.Aortic_1.Segments = 1;
REF.Aortic_2.Segments = 2;
REF.Aortic_3.Segments = 3;
REF.Aortic_4.Segments = 4;
REF.Brachio.Segments = 5;
REF.L_Carotid.Segments = 6;
REF.L_Subclavian.Segments = 7;

REF.Nodes		= [1 2 3 4 2 3 4; 2:8]';	%[inlet, outlet]
REF.L			= Length*1e-3;	%[m]
REF.Rad_prox	= sqrt(Proximal_Area/pi)*1e-3;	%mean radius; [m]
REF.Rad_dist	= sqrt(Distal_Area/pi)*1e-3;	%mean radius; [m]
REF.c_prox		= REF.PWV;						%PWV at mean radius; [m/s]
REF.c_dist		= REF.PWV;						%PWV at mean radius; [m/s]

REF.Brachio.A_out		= Distal_Area(REF.Brachio.Segments)*1e-6;		%[m2]
REF.L_Carotid.A_out		= Distal_Area(REF.L_Carotid.Segments)*1e-6;		%[m2]
REF.L_Subclavian.A_out	= Distal_Area(REF.L_Subclavian.Segments)*1e-6;	%[m2]

%	Visualise geometry
% VisualiseGeometry(REF)
% VisualiseGeometry_v2(REF)

end

function [Geometry] = ReadGeometryFile(Geometry_File)
%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
delimiter = {' '};
formatSpec = '%s%[^\n\r]';

%% Open the text file.
fileID = fopen(Geometry_File,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string',  'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
	raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

% Converts text in the input cell array to numbers. Replaced non-numeric
% text with NaN.
rawData = dataArray{1};
for row=1:size(rawData, 1)
	% Create a regular expression to detect and remove non-numeric prefixes and
	% suffixes.
	regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
	try
		result = regexp(rawData(row), regexstr, 'names');
		numbers = result.numbers;
		
		% Detected commas in non-thousand locations.
		invalidThousandsSeparator = false;
		if numbers.contains(',')
			thousandsRegExp = '^[-/+]*\d+?(\,\d{3})*\.{0,1}\d*$';
			if isempty(regexp(numbers, thousandsRegExp, 'once'))
				numbers = NaN;
				invalidThousandsSeparator = true;
			end
		end
		% Convert numeric text to numbers.
		if ~invalidThousandsSeparator
			numbers = textscan(char(strrep(numbers, ',', '')), '%f');
			numericData(row, 1) = numbers{1};
			raw{row, 1} = numbers{1};
		end
	catch
		raw{row, 1} = rawData{row};
	end
end


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
Geometry = raw;

end

function VisualiseConnectivity(Matrix,Nodes_N,Segments_N)
figure, hold on, grid on, title('Connectivity matrix')
for jj = 1:Segments_N
	for kk = 1:Nodes_N
		if Matrix(kk,jj) == 0
			plot(jj,kk,'ko')
		elseif Matrix(kk,jj) == -1
			plot(jj,kk,'bv')
		elseif Matrix(kk,jj) == 1
			plot(jj,kk,'r+')
		end
	end
end
xlabel('Segments'), xticks([1:Segments_N])
ylabel('Nodes'), yticks([1:Nodes_N])
set(gca, 'FontSize', 10)
end

function VisualiseGeometry(REF)

%	Arterial lengths
figure, hold on
plot(REF.Aortic_1.Segments, REF.L(REF.Aortic_1.Segments),'ko','MarkerSize',12)
plot(REF.Aortic_2.Segments, REF.L(REF.Aortic_2.Segments),'ko','MarkerSize',12)
plot(REF.Aortic_3.Segments, REF.L(REF.Aortic_3.Segments),'ko','MarkerSize',12)
plot(REF.Aortic_4.Segments, REF.L(REF.Aortic_4.Segments),'ko','MarkerSize',12)
plot(REF.Brachio.Segments, REF.L(REF.Brachio.Segments),'x','MarkerSize',12)
plot(REF.L_Carotid.Segments, REF.L(REF.L_Carotid.Segments),'x','MarkerSize',12)
plot(REF.L_Subclavian.Segments, REF.L(REF.L_Subclavian.Segments),'x','MarkerSize',12)
hold off

legend('Aorta 1','Aorta 2','Aorta 3','Aorta 4','Brachiocephalic','Carotid','Subclavian','Location','northwest')
set(gca,'FontSize',20)
end

function VisualiseGeometry_v2(REF)
Diam_prox_max = 2*max(REF.Rad_prox);

figure, hold on
for jj = 1:length(REF.L)
x{1,1} = [0, REF.L(jj)];
x{2,1} = [0, REF.L(jj)];
y{1,1} = [jj + REF.Rad_prox(jj)/Diam_prox_max, jj + REF.Rad_dist(jj)/Diam_prox_max]; 
y{2,1} = [jj - REF.Rad_prox(jj)/Diam_prox_max, jj - REF.Rad_dist(jj)/Diam_prox_max]; 

ShadedPlots(x,y,1);
end

xlabel('Segment length [m]')
names = {'Asc Aorta','Aortic arch I','Aortic arch II','Desc Aorta','Brachiocephalic','Carotid','Subclavian'};
set(gca,'ytick',[1:7],'yticklabel',names)
end
