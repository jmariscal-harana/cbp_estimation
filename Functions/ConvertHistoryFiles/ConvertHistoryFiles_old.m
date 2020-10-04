function data = ConvertHistoryFiles(up)
% Convert .his files created by Nektar into Matlab format.
%    Inputs:    up.dir      - the directory containing .his files. If this
%                               is not specified then the user is prompted
%                               to select the directory manually.
%               up.filename - a cell containing the filename(s) of .his
%                               files to be imported. If this is not
%                               specified then all .his files within the
%                               chosen directory are imported.
%               up.all_beats - a logical indicating whether data should be
%                               extracted from all beats (= true), or just
%                               the final beat of the simulation (= false).
%               all inputs are optional
%
%    Outputs:   simulation_data.m - a single Matlab file containing the
%                               data from all of the imported .his files,
%                               saved in the chosen directory.
%
% Peter H Charlton, King's College London, January 2018

% setup universal parameters
if nargin == 0
    up = struct;
end
up = setup_up(up);

% order according to simulation number
counter = 0;
for s = 1 : length(up.filename)
    if strcmp(up.filename{s}(1:4), 'sim_')
        counter = counter+1;
    end
end
dont_use = true;
sim_nos = nan(length(up.filename),1);
if counter == length(up.filename)
    dont_use = false;
    for s = 1 : length(up.filename)
        temp = strfind(up.filename{s}, '_');
        if length(temp) ~=2
            dont_use = true;
        else
            sim_nos(s) = str2double(up.filename{s}(temp(1)+1:temp(2)-1));
        end
        clear temp
    end
end
clear counter s

if ~dont_use
    [~, order] = sort(sim_nos);
    up.filename = up.filename(order);
end
clear dont_use order sim_nos

% convert files into matlab data
for file_no = 1 : length(up.filename)
    curr_file = up.filename{file_no};
    fprintf(['\n - ' curr_file])
    curr_file_path =[up.dir, curr_file];
    data_from_files{file_no} = ReadHistoryFiles(curr_file_path, up); 
    if ~up.all_beats
        data_from_files{file_no} = Extract_data_for_single_beat(data_from_files{file_no});
    end
    temp = strfind(curr_file, '_');
    filename_root{file_no,1} = curr_file(1:temp(end)-1); clear temp curr_file curr_file_path
end

% separate data according to filename root:
filename_root = strrep(filename_root,'-','_');
filename_roots = unique(filename_root);

char_filename_roots = filename_roots;
for s = 1 : length(filename_roots)
    if ~isstrprop(filename_roots{s}(1),'alpha')
        char_filename_roots{s} = ['sim_' filename_roots{s}];
    end    
end
clear s

for filename_root_no = 1 : length(filename_roots)
    rel_files = find(strcmp(filename_root, filename_roots{filename_root_no}));
    
    for rel_file_no = 1 : length(rel_files)
    
        % store domain no
        domain_no = up.filename{rel_files(rel_file_no)};
        temp = strfind(domain_no, '_');
        temp2 = strfind(domain_no, '.his');
        domain_no = str2double(domain_no(temp(end)+1:temp2-1));
        eval(['data.' char_filename_roots{filename_root_no} '(rel_file_no).domain_no = domain_no;'])
        clear domain_no temp temp2
            
        % store each variable's data in turn (including units)
        temp = data_from_files{rel_files(rel_file_no)};
        vars = fieldnames(temp);
        units = struct;
        for var_no = 1 : length(vars)
            eval(['data.' char_filename_roots{filename_root_no} '(rel_file_no).' vars{var_no} ' = temp.' vars{var_no} ';'])
            eval(['units.' vars{var_no} ' = find_units(''' vars{var_no} ''');']);
        end
        eval(['data.' char_filename_roots{filename_root_no} '(rel_file_no).units = units;']);
        clear units vars var_no temp
        
    end
    clear rel_file_no rel_files
    
    % sort data according to domain no
    eval(['rel_data = data.' char_filename_roots{filename_root_no} ';'])
    domain_nos = extractfield(rel_data, 'domain_no');
    [~, rel_order] = sort(domain_nos);
    rel_data = rel_data(rel_order);
    eval(['data.' char_filename_roots{filename_root_no} ' = rel_data;'])
    clear rel_data rel_order domain_nos
    
end
clear filename_root_no filename_roots filename_root data_from_files char_filename_roots

% Eliminate unwanted data
sims = fieldnames(data);
if ~up.all_data
    
    for sim_no = 1 : length(sims)
        curr_sim = sims{sim_no};
        eval(['sim_data = data.' curr_sim ';'])
        
        % eliminate unwanted signals
        curr_signals = fieldnames(sim_data);
        curr_signals = setxor(curr_signals, {'domain_no', 'distances', 'fs', 'units', 'start_sample'});
        unwanted_signals = setxor(curr_signals, up.required_signals); clear curr_signals
        sim_data = rmfield(sim_data, unwanted_signals); clear unwanted_signals
        
        % eliminated unwanted domains
        curr_domains = extractfield(sim_data, 'domain_no');
        [~,req_rows,~] = intersect(curr_domains, up.required_domains); clear curr_domains
        sim_data = sim_data(req_rows); clear req_rows
        
        % eliminate unwanted measurement points
        for domain_el = 1 : length(sim_data)
            curr_domain = sim_data(domain_el).domain_no;
            req_measurement_points = up.required_distance_els{up.required_domains == curr_domain}; clear curr_domain
            % Remove unwanted signal values
            for signal_no = 1 : length(up.required_signals)
                curr_signal = up.required_signals{signal_no};
                eval(['sim_data(domain_el).' curr_signal ' = sim_data(domain_el).' curr_signal '(:,req_measurement_points);'])
                clear curr_signal
            end
            clear signal_no
            % Remove unwanted distances
            sim_data(domain_el).distances = sim_data(domain_el).distances(req_measurement_points);
            % Remove unwanted start samples
            sim_data(domain_el).start_sample = sim_data(domain_el).start_sample(req_measurement_points);
            clear req_measurement_points
        end
        clear domain_el
        
        % store reduced data
        eval(['data.' curr_sim ' = sim_data;'])
        
    end
    
end

% save data in matlab file
savepath = [up.save_dir, 'history_files_data'];
save(savepath, 'data')

end

function up = setup_up(up)

% Update dir field if necessary
if sum(strcmp(fieldnames(up), 'dir')) && ~strcmp(up.dir(end), filesep)
    up.dir = [up.dir, filesep];
end

% Identify missing fields
required_fields = {'dir', 'filename', 'all_beats', 'all_data', 'required_domains', 'required_distance_els', 'required_signals', 'save_dir'};
current_fields = fieldnames(up);
missing_fields = setxor(required_fields, current_fields);

% Fill in missing fields
for field_no = 1 : length(missing_fields)
    curr_missing_field = missing_fields{field_no};
    switch curr_missing_field
        case 'dir'
            up.dir = uigetdir('', 'Please select the directory containing files for analysis');
            if ~up.dir
                error('No directory selected')
            end
            if ~strcmp(up.dir(end), filesep)
                up.dir = [up.dir, filesep];
            end
        case 'filename'
            temp = dir([up.dir, '*.his']); 
            if ~isempty(temp)
                up.filename = extractfield(temp, 'name');
            else
                error('There aren''t any .his files in this directory')
            end
        case 'all_beats'
            up.all_beats = true;
        case 'all_data'
            up.all_data = true;
        case 'required_domains'
            up.required_domains =      [1, 15, 21, 22, 42, 46, 49, 84, 87, 112];
        case 'required_distance_els'
            up.required_distance_els = {1, 2,  2,  3,  2,  1,  3,  3,  2,  3};
        case 'required_signals'
            up.required_signals = {'P', 'U', 'Q', 'A'}; % {'P', 'Pe', 'U', 'Q', 'A'};
        case 'save_dir'
            up.save_dir = up.dir;
    end
end

end

function history_file_data = ReadHistoryFiles(curr_file_path, up)

% Identify header lines and measurement points:
fid = fopen(curr_file_path);
header_line_log = true; line_no = 0; history_file_data.distances = [];
while header_line_log
    curr_line = fgetl(fid);
    line_no = line_no + 1;
    if ~strcmp(curr_line(1), '#')
        header_line_log = false;
    end
    if strcmp(curr_line(3:8), 'Point ')
        history_file_data.distances(end+1) = str2double(curr_line(16:23));
    end
    if strcmp(curr_line(1:3), '# t')
        header_text = strrep(curr_line, '(x,t)', '');
        header_text = strrep(header_text, ' ', '');
        header_text = strrep(header_text, '#', '');
        headers = textscan(header_text, '%s', 'Delimiter', ',');
        headers = headers{1}; clear header_text
    end
end
fclose all; 
no_header_lines = line_no - 1;
clear curr_line header_line_log fid line_no


% Import Nektar data:
raw = importdata(curr_file_path, ' ', no_header_lines);
raw = raw.data;
for col_no = 1 : length(headers)
    eval(['ext_data.' headers{col_no} ' = raw(:,col_no);']);    
end
history_file_data.fs = 1/median(diff(unique(ext_data.t)));
history_file_data.fs = round(1e6*history_file_data.fs)/1e6;
t_col = strcmp(headers, 't');
point_col = strcmp(headers, 'point');
points = raw(:, point_col);
raw = raw(:, ~t_col & ~point_col);
headers = headers(~t_col & ~point_col);

% Separate data according to each measurement location:
for header_no = 1 : length(headers)
    curr_header = headers{header_no};
    temp = raw(:,header_no);
    header_data = nan(ceil(length(temp)/length(history_file_data.distances)), length(history_file_data.distances));
    for point_no = 1 : length(history_file_data.distances)
        point_temp = temp(points == point_no);
        header_data(1:length(point_temp),point_no) = point_temp;
        clear point_temp
    end
    clear temp
    eval(['history_file_data.' curr_header ' = header_data;']); clear header_data
end

end

function units = find_units(var_name)

switch var_name
    case 'P'
        units = 'Pa';
    case 'Pe'
        units = 'Pa';
    case 'Pext'
        units = 'Pa';
    case 'U'
        units = 'm/s';
    case 'Q'
        units = 'm3/s';
    case 'A'
        units = 'm2';
    case 'distances'
        units = 'm';
    case 'fs'
        units = 'Hz';
    case 'start_sample'
        units = 'no samples';
end

end

function new_data_from_file = Extract_data_for_single_beat(data_from_file)

% cycle through each measurement point
for measurement_pt_no = 1 : length(data_from_file.distances)
    
    % identify start and end of final complete beat
    
    % Find minima in pressure (P)
    a = data_from_file.P(:,measurement_pt_no);
    minima = 3 + find( ...
        (a(3:end-2)>a(2:end-3) & a(1:end-4)>a(2:end-3)) | ...
        (a(5:end)>a(4:end-1) & a(1:end-4)>a(2:end-3)) ...
		);
    
    % Eliminate repeated minima (or those which are very close to each other)
    % - repeated
    keep_minima = true(size(minima));
    for minima_no = 2:length(minima)
        if minima(minima_no) == minima(minima_no-1)+1
            keep_minima(minima_no) = false;
        end
    end
    minima = minima(keep_minima);
    % - close to each other
    keep_minima = true(size(minima));
    for minima_no = 1:length(minima)-1
        diffs = minima - minima(minima_no);
        tol_samps = data_from_file.fs * 0.0001;  % within 0.0001 s
        if sum(diffs > 0 & diffs < tol_samps)
            keep_minima(minima_no) = false;
        end
    end
    minima = minima(keep_minima);
    
    % determine how high each minimum is in relation to the previous 2 s of
    % data
    time_int = 1; % in secs
    heights = nan(length(minima),1);
    for minimum_no = 1 : length(minima)
        prev_time_int.deb = minima(minimum_no) - round(time_int*data_from_file.fs);
        prev_time_int.fin = minima(minimum_no);
        if prev_time_int.deb <= 0
            heights(minimum_no) = nan;
        else
            heights(minimum_no) = (a(minima(minimum_no))-min(a(prev_time_int.deb:prev_time_int.fin)))/range(a(prev_time_int.deb:prev_time_int.fin));
        end
    end
    
    % identify reliable beat onsets according to amplitude
    rel_els = heights<0.125;
    beat_onsets = minima(rel_els);
    new_heights = heights(rel_els); clear rel_els
    % identify reliable beat onsets according to timing
    thresh_secs = 0.3;  % threshold in secs
    fs = data_from_file.fs;
    thresh_samps = thresh_secs*fs;
    beat_onsets_to_eliminate = [];
    for beat_onset_no = 1 : length(beat_onsets)-1
        curr_beat_onsets = beat_onsets(beat_onset_no:beat_onset_no+1);
        if range(curr_beat_onsets) < thresh_samps
            [~, temp] = min(curr_beat_onsets);
            rel_beat_onset_el = temp+beat_onset_no-1;
            beat_onsets_to_eliminate = [beat_onsets_to_eliminate, rel_beat_onset_el];
        end
    end
    rel_els = setxor(1:length(beat_onsets), beat_onsets_to_eliminate);
    beat_onsets = beat_onsets(rel_els);
    new_heights = new_heights(rel_els); clear rel_els
    
    % determine how high the max value is between each minimum and the
    % next one in relation to the prior signal
    time_int = 3; % in secs
    max_heights = nan(length(beat_onsets),1);
    for beat_onset_no = 1 : length(beat_onsets)-1
        current_beat = a(beat_onsets(beat_onset_no):beat_onsets(beat_onset_no+1));
        prev_time_int.deb = beat_onsets(beat_onset_no) - round(time_int*data_from_file.fs);
        prev_time_int.deb = max([1, prev_time_int.deb]);
        prev_time_int.fin = beat_onsets(beat_onset_no);
        max_heights(beat_onset_no) = (max(current_beat)-min(current_beat))/(max(a(prev_time_int.deb:prev_time_int.fin)) - min(a(prev_time_int.deb:prev_time_int.fin)));
    end
    
    % identify reliable beat onsets according to max height between
    rel_els = max_heights > 0.5 | isnan(max_heights);
    beat_onsets = beat_onsets(rel_els);
    
    % check that the last two beat onsets are a similar time apart
    if range(beat_onsets(end-1:end)) < 0.95*range(beat_onsets(end-2:end-1))
        beat_onsets = beat_onsets(1:end-1);
    end
    
    % check that there is some upslope after the final beat onset
    if a(beat_onsets(end)) == max(a(beat_onsets(end):end))
        beat_onsets = beat_onsets(1:end-1);
    end
    
    % check that there is sufficient signal in the relevant data for this
    % last complete beat
    if measurement_pt_no ~= 1
        curr_relevant_inds = beat_onsets(end-1):beat_onsets(end-1)+pulse_wave_duration_in_samples-1;
        if curr_relevant_inds(end) > length(a)
            beat_onsets = beat_onsets(1:end-1);
        end
        clear curr_relevant_inds
    end
            
    % identify two last beat onsets:
    if ~isempty(beat_onsets)
        last_beat_onsets = beat_onsets(end-1:end);
    else
        last_beat_onsets = [1, length(a)];
        fprintf('\n ---- Couldn''t find a reliable complete beat, so outputting all data')
    end
    
    % extract relevant data for this last complete beat
    if measurement_pt_no == 1
        relevant_inds = last_beat_onsets(1):last_beat_onsets(2)-1;
        pulse_wave_duration_in_samples = length(relevant_inds);
    else
        relevant_inds = last_beat_onsets(1):last_beat_onsets(1)+pulse_wave_duration_in_samples-1;
    end
    vars = fieldnames(data_from_file);
    for var_no = 1 : length(vars)
        curr_var = vars{var_no};
        if strcmp(curr_var, 'fs') || strcmp(curr_var, 'distances')
            eval(['new_data_from_file.' curr_var ' = data_from_file.' curr_var ';']);
            continue
        end
        eval(['new_data_from_file.' curr_var '(:,measurement_pt_no) = data_from_file.' curr_var '(relevant_inds,measurement_pt_no);']);
    end
    
    % Note the time of this beat onset
    new_data_from_file.start_sample(measurement_pt_no) = relevant_inds(1);
    
end

end