RootFolder = '~/Simulations/';
SourceFolder = [RootFolder,'2019_08_30/'];
addpath('~/Haemodynamic_Tools/Version6/Others/CopyFiles/')

Number_of_folders = 4374;

for jj = 1:Number_of_folders
    FileName = ['pwdb_',num2str(jj),'_*'];
	TargetFolder = [RootFolder,'pwdb_',num2str(jj),'/'];
	
	CopyFiles(SourceFolder,FileName,TargetFolder)
	
	fprintf('Copying folder %i of %i\n',jj,Number_of_folders)
end

restoredefaultpath