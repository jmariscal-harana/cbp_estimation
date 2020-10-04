RootFolder = '~/Haemodynamic_Tools/Version6/Aortic1D/Nektar_outputFiles/Simulations/2019_09_04/Sc2/';
addpath('~/Haemodynamic_Tools/Version6/Others/CopyFiles/')

Number_of_folders = 4374;

for jj = 1:Number_of_folders
	SourceFolder = [RootFolder,'pwdb_',num2str(jj),'/'];

	FileName = '*_1.his*';
	TargetFolder = RootFolder;
	
	CopyFiles(SourceFolder,FileName,TargetFolder)
	
	fprintf('Copying folder %i of %i\n',jj,Number_of_folders)
end

restoredefaultpath