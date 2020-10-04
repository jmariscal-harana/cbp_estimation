%%	Returns an ensemble average cycle of a raw periodic signal
%
%	INPUT
%	-signal: filtered signal to ensemble average
%	-freq: frequency of the signal (Hz)
%	-detect = 1: ensemble average is based on peak detection 
%			= 2: ensemble average is based on foot detection
%
%	OUTPUT
%	-MeanCycle: ensemble average cycle
%	-Sig_int: interpolated signal for a given frequency
%
%==========================================================================
% Marie Willemet (original author)
% August 2015
% KCL, London
%
%	Jorge Mariscal-Harana, King's College London
%	v1.0 (14/01/19) adapted from original author
%
%==========================================================================

function [MeanCycle,Sig_int] = EnsembleAverage_main(waveform, freq, filterStatus, savePath, saveName)
%% Pre-processing and fixed variables
Yor = waveform;
dt0 = 1/freq;

%	1. Apply filtering 
% ff=figure;
% plot(Yor);
% filterStatus = questdlg('Apply filtering to signal?', 'Filter?', 'Savitzky-Golay', 'Low-pass', 'None', 'None');                     
switch filterStatus
    case 'Savitzky-Golay'
        %Savitzky Golay is a low-pass filter: the lower the order or window
		%length, the lower the cut-off frequency
		addpath('~/Haemodynamic_Tools/Version6/Others/SignalFilter/')
        filtered = ApplyFiltering(Yor, 1, [2 83], 1);
        disp('Filtering with Savitzky Golay filtering (Order 2, 83 points)')
    case 'Low-pass'
        %Low-pass filter (fcut = 50, freq =1000)
		addpath('~/Haemodynamic_Tools/Version6/Others/SignalFilter/')
        filtered = ApplyFiltering(Yor, 2, [50 freq], 1);
        disp('Filtering with Low-pass filter (fcut = 50)')
    case 'None'
        filtered = Yor;
        disp('No filter applied');
    otherwise
        filtered = Yor;
        disp('No filter applied');        
end
% close(ff);

%	2. Ensemble average if signal longer than 2 seconds
if length(Yor)>2.0*freq
%     cycle = EnsembleAverage(filtered,freq, 1); %ensemble averaging based on peaks of waves
    MeanCycle = EnsembleAverage(filtered,freq, 2); %ensemble averaging based on feet of waves
else
    MeanCycle = filtered;
end

%-- Fixed data for cardiac cycle
tfinal0 = length(MeanCycle)*dt0;
t0 = [0:dt0:tfinal0-dt0]';

%-- Save processed cycle
% cd(PPGPath)
if nargin == 5
% 	dlmwrite([savePath,'/',saveName,'_ensemble.txt'],[t0 MeanCycle]);
	dlmwrite([savePath,'/',saveName,'_ensemble.txt'],[MeanCycle]);
end

%-- Offset the signal to start at onset of systole - assumes a periodic signal
% [minSig, I_minSig] = min(MeanCycle);
% ff2=figure;
% hold on
% plot(cycle);
% plot(I_minSig,minSig,'*');
% title('Select the start of cycle')
% [x,y] = ginput(1);
% x = I_minSig;
% if(round(x)<0) x=1;
% elseif(round(x)>length(MeanCycle)) x=length(MeanCycle);
% end
% Sig0_onset = [MeanCycle(round(x):end); MeanCycle(1:round(x)-1)];
% close(ff2)

%Sig0_onset = cycle;

%-- Interpolate to a spline of sampling=fsampl Hz
% fsampl = 200; %Hz
% dt=1/fsampl;
% ti = [0:dt:tfinal0-dt]';
% Sig_int = interp1(t0,Sig0_onset,ti,'spline');
% disp(['Interpolating signal to sampling of ', num2str(fsampl),' Hz'])

end