function pressure2csv(filepath)

% get filepath to the folder. than this function will read the pressures via the pressure extraction kit and
% subsequently imports it into a csv file with column 1:3 being ECG signals I II III and 4 and 5 being the two
% pressure files.pressure. B_R


if nargin < 1
    filepath = '*.000';
end

files = dir(filepath);

for j=1:length(files)
    file = char(files(j).name);
    
    if contains(filepath,'*')
        path = fileparts(filepath);
        fullname = strcat(path,file);
        cd(path)
    elseif length(files)==1 % assume the input was a filename instead of a folder, may not be true but the fileparts gets confused by the . in the filenames...
        fullname = filepath;
        path = fileparts(fullname);
        cd(path);
    else
        path = filepath;
        fullname = strcat(path,file);
        cd(path);
    end
    cmd = sprintf('python ~/Haemodynamic_Tools/Version6/Others/Extracttocsv/Cardiotek_pressureExtraction.py "%s"',fullname);
    [status,result] = system(cmd);
    fprintf('file converted, status = %d\n',status);
    
    if status == 1
        continue
    end
    
    outfiles = {'ECG_0.txt','ECG_1.txt','ECG_2.txt','Pressure_3.txt','Pressure_4.txt'};
    data = [];
    for i=1:length(outfiles)
        fid = fopen(outfiles{i});
        tline = '';
        k=0;
        while ischar(tline)
            k=k+1;
            tline = fgetl(fid);
            data_i(k,i) = (str2double(tline) - 32768)*0.01;
        end
        
        %data_i = (importdata(outfiles{i}) - 32768)*0.01;
        %data(:,i) = data_i;
%         figure(j); subplot(2,3,i); plot(data_i(:,i)); title(file);
	end
% 	delete (outfiles{:})
	
    fout = [file '.csv'];
    csvwrite(fout,data_i);
    fprintf('wrote %s\n',fout);
end
