%%	Move each virtual patient simulation to a folder with the name M116art_XXXX-XXXX
clear all

CurrentDir = cd;
SimulationDir = [cd,'/Nektar_outputFiles/Simulations/v1/opt/'];
addpath(SimulationDir);
fileID = fopen('Nektar_Success.txt');
SimulationName = fgetl(fileID);

while ischar(SimulationName)	%while fgetl returns a simulation name
	FullFolder = fullfile(SimulationDir,SimulationName); %full path to the simulation folder
	if exist(FullFolder,'dir')	%check if folder exists
		FolderContent = dir(FullFolder);	%check folder content
		if numel(FolderContent) >= 1
			%folder exists and is not empty
			warning(['Folder ',SimulationName,' is not empty: files will '...
				'not be copied to avoid overwriting'])
		else
			% 			folder exists and is empty
			movefile([SimulationDir,SimulationName,'.*'],[SimulationDir,SimulationName])	%move other files
			movefile([SimulationDir,SimulationName,'_*'],[SimulationDir,SimulationName])	%move .his files
			movefile([SimulationDir,SimulationName(1:end-3),'BC_opt.*'],[SimulationDir,SimulationName])	%move other files
		end
	else
		%if folder does not exist
		mkdir([SimulationDir,SimulationName])	%create folder
		movefile([SimulationDir,SimulationName,'.*'],[SimulationDir,SimulationName])	%move other files
		movefile([SimulationDir,SimulationName,'_*'],[SimulationDir,SimulationName])	%move .his files
		movefile([SimulationDir,SimulationName(1:end-3),'BC_opt.*'],[SimulationDir,SimulationName])	%move other files
	end
	SimulationName = fgetl(fileID);	%read next line from Nektar_Success.txt
end

fclose(fileID);
