%%	Calculate quartiles
% if there are 4n datapoints, divide the data in 4 equal subsets
% if there are (4n+1) datapoints, ignore the median
% if there are (4n+2) datapoints, include Q1 and Q3 in the 1st/2nd and
% 3rd/4th quartiles, respectively
% if there are (4n+3) datapoints, include the median twice
%
%	Input
%	-Data_in: column/row vector containing data to be split in quartiles
%	-Data_in: 2xN/Nx2 matrix containing (1) data to be split in quartiles,
%	(2) corresponding data to be plotted instead
%
%	Output
%	-Data_out: 4 equal-sized subsets
%	-ID_out: indices for datapoints within each subset
%
%==========================================================================
%	Jorge Mariscal-Harana, King's College London (original author)
%	v1.0 (20/09/19)
%
%==========================================================================

function [Data_out,ID_out] = Quartiles(Data_in,Plots)

if min(size(Data_in)) == 1
	Data_1 = Data_in;
elseif min(size(Data_in)) == 2
	if size(Data_in,1) == 2
		Data_1 = Data_in(1,:);
		Data_2 = Data_in(2,:);
	elseif size(Data_in,2) == 2
		Data_1 = Data_in(:,1);
		Data_2 = Data_in(:,2); 
	end
	Data_2 = Data_2(:);
else
	error('The size of the smallest matrix dimension can only be 1 or 2')
end
Data_1 = Data_1(:);

Data = sort(Data_1);
N = length(Data);

Q(1) = round(N/4);
Q(2) = round(N/2);
Q(3) = round(3*N/4);

if mod(N,4) == 0
	Data_0 = Data(1:Q(1));
	Data_25 = Data(Q(1)+1:Q(2));
	Data_50 = Data(Q(2)+1:Q(3));
	Data_75 = Data(Q(3)+1:end);
	
elseif mod(length(Data),4) == 1
	Data_0 = Data(1:Q(1));
	Data_25 = Data(Q(1)+1:Q(2)-1);
	Data_50 = Data(Q(2)+1:Q(3));
	Data_75 = Data(Q(3)+1:end);
	
elseif mod(length(Data),4) == 2
	Data_0 = Data(1:Q(1));
	Data_25 = Data(Q(1):Q(2));
	Data_50 = Data(Q(2)+1:Q(3));
	Data_75 = Data(Q(3):end);

elseif mod(length(Data),4) == 3
	Data_0 = Data(1:Q(1));
	Data_25 = Data(Q(1)+1:Q(2));
	Data_50 = Data(Q(2):Q(3));
	Data_75 = Data(Q(3)+1:end);
	
else
	error('Check input data')
end

ID_0 = find (Data_1 <= Data_0(end));
ID_25 = find (Data_1 >= Data_25(1) & Data_1 <= Data_25(end));
ID_50 = find (Data_1 >= Data_50(1) & Data_1 <= Data_50(end));
ID_75 = find (Data_1 >= Data_75(1));

if Plots == 1
	figure
	Data_0_x = 1*ones(length(Data_0),1);
	Data_25_x = 2*ones(length(Data_25),1);
	Data_50_x = 3*ones(length(Data_50),1);
	Data_75_x = 4*ones(length(Data_75),1);
	plot([Data_0_x,Data_25_x,Data_50_x,Data_75_x],[Data_0,Data_25,Data_50,Data_75],'o')
elseif Plots == 2
	figure
	boxplot([Data_0,Data_25,Data_50,Data_75],'Labels',{'0-25%','25-50%','50-75%','75-100%'})
elseif Plots == 3
	figure
	boxplot([Data_2(ID_0),Data_2(ID_25),Data_2(ID_50),Data_2(ID_75)],'Labels',{'0-25%','25-50%','50-75%','75-100%'})
end

Data_out = [Data_0,Data_25,Data_50,Data_75];
ID_out = [ID_0;ID_25;ID_50;ID_75]';

end