%%	Read .his files for each subject
function [Data] = ReadHis(SimulationFolder)
fileID		= fopen(SimulationFolder);
Data		= textscan(fileID,'%f %f %f %f %f %f %f','CommentStyle','#');
Data		= [Data{:}];

if 	Data(1,7) == 1
elseif Data(1,6) == 1
	Data = Data(:,1:6);
else
	error('.his format not supported')
end

end