%% ABP beat detection (accounting for respiration)
% P H Charlton, KCL, 13th June 2017
%
% This script detects the beats and breaths in the ABP signals, and then calculates
% systolic, mean and diastolic BPs.
%
% It is based on Sun's implementation of Zong's ABP beat detector. The
% original article is at:
%   http://ecg.mit.edu/george/publications/abp-cinc-2003.pdf
% Sun's implementation is at:
%   https://physionet.org/physiotools/cardiac-output/code/2analyze/wabp.m
%
% To run the script use the command:
%
%    [final_data] = perform_abp_beat_detection_resp;
%
% It will ask you where the "curves_ext.xlsx" file is. Once you have selected
% the file then it should show you some plots. Just close each plot to move
% on to the next one.
%
% Note that the BPs can be shown by looking at e.g. final_data.abp_1
%
% Licence: this script uses some code from the RRest toolbox, which is
% covered by the GNU GPL.

function [final_data] = perform_abp_beat_detection_resp

% setup universal parameters (contains settings)
up = setup_up;

% load data
data = load_raw_data(up);

% downsample data to 125 Hz
data_ds = downsample_data(data, up);

% detect beats (including plotting pulse onsets and peaks)
data_with_beats = detect_beats(data_ds, up);

% detect resp cycles
data_with_breaths = detect_resp_cycles(data_with_beats, up);

% find BPs
final_data = calc_BPs(data_with_breaths, up);

% Note that the BPs can be shown by looking at e.g. final_data.abp_1

end

function up = setup_up

[FileName,PathName,~] = uigetfile('curves.xlsx', 'Please select the "curves.xlsx" file');
up.paths.data_file = [PathName, FileName]; % 'C:\Documents\Data\Arna_ABP\curves_ext.xlsx'
%up.paths.data_file = 'C:\Documents\Data\Arna_ABP\curves_ext.xlsx';

% sampling frequency required by Zong's ABP beat detector
up.downsample_fs = 125;  

% tolerance in secs for beat onset detection (searches in a region of twice this width either side of the detected onset)
up.beat_detection.onset_tolerance = 0.1;  

% set to 1 to do plots, or 0 to not do them
up.do_plots = 1;  

% choose how to calculate the average BPs (mean or median)
up.bp_ave_method = 'mean';

% Resp filtering

% Filter characteristics: Eliminate HFs (above resp freqs)
%up.paramSet.elim_hf.Fpass = 1.2;  % in Hz
%up.paramSet.elim_hf.Fstop = 0.8899;  % in Hz     (1.2 and 0.8899 provide a -3dB cutoff of 1 Hz)
up.paramSet.elim_hf.Fpass = 0.8;  % in Hz
up.paramSet.elim_hf.Fstop = 0.4;  % in Hz     (0.661 and 0.567 provide a -3dB cutoff of 0.6 Hz)
up.paramSet.elim_hf.Dpass = 0.05;
up.paramSet.elim_hf.Dstop = 0.01;
% Filter characteristics: Eliminate VLFs (below resp freqs): For 4bpm cutoff
up.paramSet.elim_vlf.Fpass = 0.157;  % in Hz
up.paramSet.elim_vlf.Fstop = 0.02;   % in Hz     (0.157 and 0.02 provide a - 3dB cutoff of 0.0665 Hz)
up.paramSet.elim_vlf.Dpass = 0.05;
up.paramSet.elim_vlf.Dstop = 0.01;
% Filter characteristics: duration of Tukey window taper in secs
up.paramSet.tukey_win_duration_taper = 2;

end

function data = load_raw_data(up)

% load raw data
[num,txt,raw] = xlsread(up.paths.data_file);

% extract each of the four signals
col_offset = 5;
for s = 1 : 4
    temp = num(:,col_offset+s);
    eval(['data.abp_' num2str(s) ' = temp(~isnan(temp));']);
end

% select relevant portion of each signal
for s = 1 : 4
    
    % extract current data
    eval(['curr_data = data.abp_' num2str(s) ';']);
    
    % extract relevant data
    if s == 1
        rel_els = 25680:length(curr_data);
    elseif s == 2
        rel_els = 1:length(curr_data);
    elseif s == 3
        rel_els = 1:36580;
    elseif s == 4
        rel_els = 1:22670;
    end
    rel_data = curr_data(rel_els); clear curr_data rel_els
    
    % store relevant data
    eval(['data.abp_' num2str(s) ' = rel_data;']); clear rel_data
    
end

end

function data_ds = downsample_data(data, up)

% assume sampling freq is 1000 Hz
initial_fs = 1000;

no_sigs = length(fieldnames(data));
for s = 1 : no_sigs
    
    % load this signal
    eval(['curr_sig.v = data.abp_' num2str(s) ';'])
    curr_sig.fs = initial_fs;
    curr_sig.t = [0:length(curr_sig.v)-1]/curr_sig.fs;
    
    % downsample
    decimation_factor = curr_sig.fs / up.downsample_fs;
    ds_sig.fs = up.downsample_fs;
    ds_sig.v = decimate(curr_sig.v(:)', decimation_factor);
    ds_sig.t = [0:length(ds_sig.v)-1]/ds_sig.fs;
    
    % store data
    eval(['data_ds.abp_' num2str(s) ' = ds_sig;'])
    clear ds_sig
    
end

end

function data_with_beats = detect_beats(data_ds, up)

% detect beats in each signal

% cycle through each signal
no_sigs = length(fieldnames(data_ds));
for sig_no = 1 : no_sigs
    
    % extract data for this signal
    eval(['curr_sig = data_ds.abp_' num2str(sig_no) ';']);
    
    % detect pulse onsets in this signal
    candidate_onsets = wabp_pc(curr_sig.v);
    
    % identify the tolerance (in samples) of the onset detection
    tol = floor(up.beat_detection.onset_tolerance*curr_sig.fs);
    
    % refine pulse onsets (as they are not always right at the foot of the pulse)
    refined_onsets = nan(length(candidate_onsets),1);
    for beat_no = 1 : length(candidate_onsets)
        curr_onset_index = candidate_onsets(beat_no);
        possible_indices.i = curr_onset_index - tol : curr_onset_index + tol; clear curr_onset_index
        possible_indices.v = curr_sig.v(possible_indices.i);
        [~, rel_el] = min(possible_indices.v);
        refined_onset_index = possible_indices.i(rel_el); clear rel_el possible_indices
        refined_onsets(beat_no) = refined_onset_index; clear refined_onset_index
    end
    clear tol beat_no candidate_onsets
    
    % find pulse peaks
    pulse_peaks = nan(length(refined_onsets)-1,1);
    for beat_no = 1 : length(refined_onsets)-1
        [~, el] = max(curr_sig.v(refined_onsets(beat_no):refined_onsets(beat_no+1)));
        pulse_peaks(beat_no) = refined_onsets(beat_no) + el -1; clear el
    end
    clear beat_no

    % store data for this signal
    curr_sig.onsets = refined_onsets; clear refined_onsets
    curr_sig.peaks = pulse_peaks; clear pulse_onsets
    
    % do plot if desired
    if up.do_plots
        plot(curr_sig.t, curr_sig.v), hold on, plot(curr_sig.t(curr_sig.peaks), curr_sig.v(curr_sig.peaks), 'or'), plot(curr_sig.t(curr_sig.onsets), curr_sig.v(curr_sig.onsets), 'ok')
        title('Beat detection: onsets (black), peaks (red)')
        uiwait(gcf)
    end
    
    % store data
    eval(['data_with_beats.abp_' num2str(sig_no) ' = curr_sig;']); clear curr_sig
end

end

function data_with_breaths = detect_resp_cycles(data_with_beats, up)

% cycle through each signal
no_sigs = length(fieldnames(data_with_beats));
for sig_no = 1 : no_sigs
    
    % extract data for this signal
    eval(['curr_sig = data_with_beats.abp_' num2str(sig_no) ';']);
    
    % extract breathing signal
    s = curr_sig;
    breathing = bpf_signal_to_remove_non_resp_freqs(s, s.fs, up);
    
    % detect individual breaths
    [~,breaths] = CtO(breathing, up);
    
    % extract data for this signal
    curr_sig.breathing = breathing; clear breathing
    curr_sig.breaths = breaths; clear breaths
    %curr_sig.rr = rr; clear rr
    eval(['data_with_breaths.abp_' num2str(sig_no) ' = curr_sig;']);
    
    % Make plot
    breath_offset = 0;  % in secs
    breath_ht = max(curr_sig.breathing.v)+4; % in mmHg
    if up.do_plots
        plot(curr_sig.t, curr_sig.v-mean(curr_sig.v)), hold on,
        plot(curr_sig.breathing.t, curr_sig.breathing.v, 'r', 'LineWidth', 2),
        for breath_no = 1 : length(curr_sig.breaths.start)
            t = [curr_sig.breathing.t(curr_sig.breaths.start(breath_no)) + breath_offset, ...
                curr_sig.breathing.t(curr_sig.breaths.end(breath_no)) - breath_offset];
            v = breath_ht*ones(1,2);
            plot(t, v, 'o-k', 'LineWidth', 2)
        end
        title('Breath detection: breaths (black) detected from breathing (red)')
        uiwait(gcf)
    end
    clear curr_sig
    
end


end

function s_filt = bpf_signal_to_remove_non_resp_freqs(s, Fs, up)
% Filter pre-processed signal to remove freqs outside of the range of plausible resp freqs

s.fs = Fs;

%% Window signal to reduce edge effects
s_win = tukey_win_data(s,up);

%% HPF to eliminate VLFs
s_evlfs = s_win;
%s_evlfs.v = elim_vlfs(s_win, up);

%% LPF to remove freqs above resp
s_filt.t = s_evlfs.t;
s_filt.v = lp_filter_signal_to_remove_freqs_above_resp(s_evlfs.v, s_evlfs.fs, up);
s_filt.fs = s_evlfs.fs;

end

function d_s_win = tukey_win_data(d_s, up)
%TUKEY_WIN_DATA extracts respiratory activity using a BPF (PC's implementation)
%
%	d_s_win = tukey_win_data(data, up)
%
%	Inputs:
%       data            raw signal data
%       up              universal parameters structure
%
%	Outputs:
%       d_s_win             a signal after Tukey windowing
%

%% Window Signal to reduce edge effects

duration_of_signal = d_s.t(end) - d_s.t(1);
prop_of_win_in_outer_regions = 2*up.paramSet.tukey_win_duration_taper/duration_of_signal;
tukey_win = tukeywin(length(d_s.v), prop_of_win_in_outer_regions);
d_s_win = d_s;    % copy time and fs
d_s_win.v = detrend(d_s.v(:)).*tukey_win(:);

end

function s_filt = lp_filter_signal_to_remove_freqs_above_resp(s, Fs, up)
%% Filter pre-processed signal to remove freqs above resp

% parameters for the low-pass filter to be used
flag  = 'scale';
Dpass = up.paramSet.elim_hf.Dpass;
Dstop = up.paramSet.elim_hf.Dstop;
Fstop = up.paramSet.elim_hf.Fstop;
Fpass = up.paramSet.elim_hf.Fpass;

% create filter
[N,Wn,BETA,TYPE] = kaiserord([Fstop Fpass]/(Fs/2), [1 0], [Dstop Dpass]);
b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);
AMfilter = dfilt.dffir(b);

%% Check frequency response
% Gives a -3 dB cutoff at ? Hz, using:
% freqz(AMfilter.Numerator)
% norm_cutoff_freq = 0.3998;    % insert freq here from plot
% cutoff_freq = norm_cutoff_freq*(Fs/2);

% Prepare signal
s_dt=detrend(s);
s_filt = filtfilt(AMfilter.numerator, 1, s_dt);
end

function [rr, breaths] = CtO(data, up)
%CtO estimates the RR using an implementation of Count Orig, as specified in:
% A. Schäfer et al., “Estimation of breathing rate from respiratory sinus arrhythmia: comparison of various methods,” Ann. Biomed. Eng., vol. 36, no. 3, pp. 476–85, Mar. 2008.
%	            CtO(option, up)
%
%	Inputs:
%       rel_data    .t  vector of times
%                   .v  vector of resp Sig values
%                   .fs sampling freq
%       up              universal parameters structure
%       wins        .t  vector of start times
%
%	Outputs:
%       rr          .t  vector of times of estimated RRs
%                   .v  vector of estimated RRs
%
% modified on 13 June 2017 by Pete Charlton to ignore windows and store breaths

%% identify peaks
diffs_on_left_of_pt = diff(data.v); diffs_on_left_of_pt = diffs_on_left_of_pt(1:(end-1)); diffs_on_left_of_pt = logical(diffs_on_left_of_pt>0);
diffs_on_right_of_pt = diff(data.v); diffs_on_right_of_pt = diffs_on_right_of_pt(2:end); diffs_on_right_of_pt = logical(diffs_on_right_of_pt<0);
peaks = find(diffs_on_left_of_pt & diffs_on_right_of_pt)+1;

%% identify troughs
diffs_on_left_of_pt = diff(data.v); diffs_on_left_of_pt = diffs_on_left_of_pt(1:(end-1)); diffs_on_left_of_pt = logical(diffs_on_left_of_pt<0);
diffs_on_right_of_pt = diff(data.v); diffs_on_right_of_pt = diffs_on_right_of_pt(2:end); diffs_on_right_of_pt = logical(diffs_on_right_of_pt>0);
troughs = find(diffs_on_left_of_pt & diffs_on_right_of_pt)+1;

%% define threshold
% q3 = quantile(data.v(peaks), 0.75);
q3 = prctile(data.v(peaks), 75);
thresh = 0.2*q3;

%% find relevant peaks and troughs
extrema = sort([peaks(:); troughs(:)]);
rel_peaks = peaks(data.v(peaks) > thresh);
rel_troughs = troughs(data.v(troughs) < 0);

%% find valid breathing cycles
% valid cycles start with a peak:
valid_cycles = zeros(length(rel_peaks)-1,1);
cycle_durations = nan(length(rel_peaks)-1,1);
counter_no = 1;
for peak_no = 1 : (length(rel_peaks)-1)
    
    % valid if there is only one rel_trough between this peak and the
    % next
    cycle_rel_troughs = rel_troughs(rel_troughs > rel_peaks(peak_no) & rel_troughs < rel_peaks(peak_no+1));
    if length(cycle_rel_troughs) == 1
        valid_cycles(peak_no) = 1;
        breaths.start(counter_no) = rel_peaks(peak_no);
        breaths.end(counter_no) = rel_peaks(peak_no+1); counter_no = counter_no+1;
        cycle_durations(peak_no) = data.t(rel_peaks(peak_no+1)) - data.t(rel_peaks(peak_no));
    end
    
end

%% Calc RR
if sum(valid_cycles) == 0
    rr = nan;
else
    % Using average breath length
    ave_breath_duration = nanmean(cycle_durations);
    rr = 60/ave_breath_duration;
end

% To check breaths:
%plot(data.t, data.v), hold on, plot(data.t(breaths.start), data.v(breaths.start), '.r'), plot(data.t(breaths.end), data.v(breaths.end), 'ok')

end

function final_data = calc_BPs(data_with_beats, up)

% cycle through each signal
no_sigs = length(fieldnames(data_with_beats));
for sig_no = 1 : no_sigs
    
    % extract data for this signal
    eval(['curr_sig = data_with_beats.abp_' num2str(sig_no) ';']);
    
    % identify those times which are within valid breaths
    valid_i = false(length(curr_sig.t),1);
    for breath_no = 1 : length(curr_sig.breaths.start)
        valid_i(curr_sig.breaths.start(breath_no) : curr_sig.breaths.end(breath_no)) = true;        
    end
    
    % extract BPs (only from during valid breaths)
    [dbp,sbp,mbp] = deal([]);
    counter_no = 0;
    used_beat = false(length(curr_sig.peaks),1);
    for beat_no = 1 : length(curr_sig.peaks)
        curr_onset = curr_sig.onsets(beat_no);
        curr_peak = curr_sig.peaks(beat_no);
        curr_end = curr_sig.onsets(beat_no+1)-1;
        
        if mean(valid_i(curr_onset:curr_end)) == 1
            counter_no = counter_no+1;
            sbp(counter_no) = curr_sig.v(curr_peak);
            dbp(counter_no) = curr_sig.v(curr_end);
            mbp(counter_no) = mean(curr_sig.v(curr_onset:curr_end));
            used_beat(beat_no) = true;
        end
    end
    
    % store data
    if strcmp(up.bp_ave_method, 'mean')
        curr_sig.dbp = mean(dbp);
        curr_sig.mbp = mean(mbp);
        curr_sig.sbp = mean(sbp);
    elseif strcmp(up.bp_ave_method, 'median')
        curr_sig.dbp = median(dbp);
        curr_sig.mbp = median(mbp);
        curr_sig.sbp = median(sbp);
    end
     
    % do plot if desired
    if up.do_plots
        plot(curr_sig.t, curr_sig.v), hold on,
        title('Calculated BPs: DBP (red), MBP (black), and SBP (green)')
        rel_els = curr_sig.onsets(find(used_beat)+1); % DBP
        plot(curr_sig.t(rel_els), curr_sig.v(rel_els), 'or')
        rel_els = curr_sig.peaks(find(used_beat)); % SBP
        plot(curr_sig.t(rel_els), curr_sig.v(rel_els), 'og')
        rel_els = curr_sig.peaks(find(used_beat)); % MBP
        plot(curr_sig.t(rel_els), mbp, 'ok')
        plot(curr_sig.t([1,end]), ones(2,1)*curr_sig.mbp, 'k', 'LineWidth', 2)
        plot(curr_sig.t([1,end]), ones(2,1)*curr_sig.sbp, 'g', 'LineWidth', 2)
        plot(curr_sig.t([1,end]), ones(2,1)*curr_sig.dbp, 'r', 'LineWidth', 2)
        uiwait(gcf)
    end
    
    % store data
    eval(['final_data.abp_' num2str(sig_no) ' = curr_sig;']); clear curr_sig
end

end

function r = wabp_pc(abp)
% WABP  ABP waveform onset detector.
%   r = WABP(ABP) obtains the onset time (in samples) 
%       of each beat in the ABP waveform.
%
%   In:   ABP (125Hz sampled)
%   Out:  Onset sample time
% 
%   Usage:
%   - ABP waveform must have units of mmHg
%
%   Written by James Sun (xinsun@mit.edu) on Nov 19, 2005.  This ABP onset
%   detector is adapted from Dr. Wei Zong's wabp.c.
%
%   Modified by PC on 29-Oct-2015 as follows:
%   - Removed input checks
%   - Added step to ensure correct orientation of input vector
%   - Added/modified comments to try to make the code clearer to the non-expert
%   - Changed period of initial threshold calculation to 10 s to be in keeping with the paper.
%   - Removed a "+5" in the initial "if" clause to determine if SSF has exceeded threshold.
%   - Changed window length around the SSF threshold-crossing point to be in keeping with the paper.
%   - Added check to see if there is at least 10 secs of data.

%% Setup

% ensure correct orientation of input vector
abp = abp(:);

% scale physiologic ABP
offset   = 1600;
scale    = 20;
Araw     = abp*scale-offset;

%% LPF
A = filter([1 0 0 0 0 -2 0 0 0 0 1],[1 -2 1],Araw)/24+30;
A = (A(4:end)+offset)/scale; % Takes care of 4 sample group delay

%% Slope-sum function
dypos          = diff(A);
dypos(dypos<0) = 0;
ssf            = [0; 0; conv(ones(16,1),dypos)];   % Produce analysing window of 128 ms

%% Decision rule
if length(ssf)<(125*10)
    % less than 10 secs of data so don't perform beat detection
    r = [];
    return
end

avg0       = mean(ssf(1:(125*10)));   % average of 1st 10 seconds of SSF (changed by PC)
Threshold0 = 3*avg0;                  % initial decision threshold

% ignoring "learning period" for now
lockout    = 0;    % lockout >0 means we are in refractory
timer      = 0;    % Timer used for counting time since previous ABP pulse
z          = zeros(100000,1);   % Stores the indices of the detected beats
counter    = 0;    % Counter of number of beats detected

for t = 50:length(ssf)-19       % PC changed to 19 to be in keeping with the change below (which satisfies 150 ms)
    lockout = lockout - 1;      % increase lockout by 1 sample
    timer   = timer   + 1;      % Timer used for counting time after previous ABP pulse

    if (lockout<1) && (ssf(t)>avg0)  % Not in refractory and SSF has exceeded threshold here (pc changed from "(lockout<1) & (ssf(t)>avg0)")
        timer = 0;
        maxSSF = max(ssf(t:t+19));  % Find local max of SSF
        minSSF = min(ssf(t-19:t));  % Find local min of SSF  - pc changed to 19 to fit in with the 150 ms in the paper
        if maxSSF > (minSSF + 10)
            onset = 0.01*maxSSF ;  % Onset is at the time in which local SSF just exceeds 0.01*maxSSF

            tt       = t-16:t;
            dssf     = ssf(tt) - ssf(tt-1);
            BeatTime = find(dssf<onset,1,'last')+t-17;
            counter  = counter+1;

            if isempty(BeatTime)
                counter = counter-1;
            else
                z(counter) = BeatTime;
            end
            Threshold0 = Threshold0 + 0.1*(maxSSF - Threshold0);  % adjust threshold
            avg0 = Threshold0 / 3;        % adjust avg

            lockout = 32;   % lock so prevent sensing right after detection (refractory period)
        end
    end

    if timer > 312  % Lower threshold if no pulse detection for a while
        Threshold0 = Threshold0 - 1;
        avg0       = Threshold0/3;
    end
end
r = z(find(z))-2;

end