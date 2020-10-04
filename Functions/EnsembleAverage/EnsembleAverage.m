function MeanCycle = EnsembleAverage(signal, freq, detect)
%
% Returns an ensemble average cycle of a filtered period signal 
% signal: filtered signal to ensemble average
% freq: frequency of the signal (Hz)
% detect = 1: ensemble average is based on peak detection 
%        = 2: ensemble average is based on foot detection
%
% Marie Willemet
% August 2015
% KCL, London

%
%	Jorge Mariscal-Harana, King's College London
%	v1.0 (22/08/19) adapted from original author - Remove cycles which are 
%	25% longer or shorther than the median cycle;


% === Example
% load('example/PPG_finger_filtered.txt');
% MeanCycle = EnsembleAverage(PPG_finger_filtered, 1000,1);

TTpath = '~/Haemodynamic_Tools/Version6/Others/PWV_TT/';

%-- Select window of data to analyse
%	Manual selection
% figure(22)
% hold on 
% plot(signal,'color',[0.5 0.5 0.5])
% disp('Ensemble averaging: click to select start and end of data to include in ensemble average')
% [xcut,ycut] = ginput(2);
%	Entire signal
xcut = [1, length(signal)];


if(round(xcut(1))<0) xcut(1)=1;
elseif(round(xcut(1))>length(signal)) xcut(1)=length(signal);end
if(round(xcut(2))<0) xcut(2)=1;
elseif(round(xcut(2))>length(signal)) xcut(2)=length(signal);end

Y = signal(round(xcut(1)):round(xcut(2)));
% close(figure(22))

% figure(99)
% title('Signal used to compute ensemble average')
% plot(Y)

%-- Detect landmarks of waveform (feet, max) using TTAlgorithm_b (modified output)
Y1 = Y(1:end-100);
Y2 = Y(101:end);

addpath(TTpath)
algo=1;	%foot-to-foot
wavetype=1;	%pressure
[TT, t_foot, i_max] =  TTAlgorithm_b([Y1, Y2],freq,algo,wavetype,1,0);
% TTAlgorithm_b(signal,f,algorithm,waveform,continuous,show)

%flip vectors to be in increasing time
t_foot = fliplr(t_foot);
i_max = fliplr(i_max);

%-- Compute mean cycle 
% Wave = Y1; Ix = 2; (second waveform with offset)
Wave = Y2; Ix = 1; 

% Average based on the peak/feet of wave
if detect == 1 %peak detection
    detect_mark = 'peaks';
    delay = 130; %ms
    index = i_max(Ix,:);
    disp(['Ensemble average based on peaks detection, delay from peak = ', num2str(delay/freq*1000),' ms'])
elseif detect == 2 %foot detection
    detect_mark = 'feet';
    delay = 30;	%[ms]
    index = round(t_foot(Ix,:).*freq);
    disp(['Ensemble average based on feet detection, delay from foot = ', num2str(delay/freq*1000),' ms'])
else
    errordlg('Select if detection if based on peaks (detect=1), or on feet of wave (detect=2)','detect parameter')
end

% % add the end index of last cycle analysed ==> this includes a cycle that is not a period long
% index(end+1) = min(index(end)+round(mean(diff(index))), length(Wave));

nb=size(t_foot,2)-1; %number of cycles

% Extract individual cycles
figure(21), hold on
for k = 1:nb
    SingleCycle{k} = Wave((index(k)-delay):(index(k+1)-delay)); %Raw individual cycles
	WaveLength(k) = length(SingleCycle{k});
	plot(SingleCycle{k},'color', [0.5 0.5 0.5]);
end 

%Remove cycles which are 25% longer or shorther than the median cycle
warning('Removing cycles 25% longer or shorther than the median cycle')
RemoveCycles = [~(WaveLength > 1.25*median(WaveLength) | WaveLength < 0.75*median(WaveLength))];
SingleCycle = {SingleCycle{RemoveCycles}};
nb = length(SingleCycle);

% Update lengths
for k = 1:nb
	WaveLength(k) = length(SingleCycle{k});
end

%Truncate to a given reference cycle
TruncateLength = min(WaveLength);

figure(22), hold on
j=1;
for k = 1:nb
	t = [0:1/freq:(length(SingleCycle{k}) - 1)*1/freq];
	plot(t,SingleCycle{k},'k','LineWidth',0.5);
    UniqueSingleCycle(:,j) = SingleCycle{k}(1:TruncateLength); 
    j=j+1;
end

%Ensemble average:
MeanCycle = mean(UniqueSingleCycle,2);
SDCycle = std(UniqueSingleCycle,0,2); %standard deviation
t = [0:1/freq:(length(MeanCycle) - 1)*1/freq];
% plot(t,SDCycle,'-r')
plot(t,MeanCycle,'b', 'linewidth',2)
disp(['   Ensemble average of ',num2str(nb),' cycles']);
% set(gca, 'fontsize',24)
% title(['Detection of ', detect_mark]) 


%% Check if individual cycles need to be removed from ensemble average 
% HAPPY = char(questdlg({'Do you need to remove individual cycles from average?'},'Improve average?',{'YES'},{'NO'},{'NO'}) );
HAPPY = 'NO';
col = {'k','b'};

while(strcmp(HAPPY,'YES'))  
    figure(44)
    hold on
    
    ix_all = zeros(max(diff(index)),nb); %a matrix to sort indices of all cycles
    for k = 1:nb
        ind = [index(k):index(k+1)];
        plot(ind, Wave((index(k)-delay):(index(k+1)-delay)) ,'color',col{mod(k,2)+1}, 'linewidth',1);
        ix_all([1:length(ind)],k) = ind;
    end

    % How many cycles to remove?
    defaultanswer={'1'};
    removeNb = inputdlg('How many cycles need to be removed from average?','Number',1,defaultanswer);
    if (numel(removeNb) ~= 1)
        err = errordlg('Number?','Erreur','modal');
        return
    end
    remNb = str2num(char(removeNb)); %convert to number

    % Click on plot to select cycles to remove
    disp('Click on plot to select cycles to remove')
    figure(44)
    [xcut,ycut] = ginput(remNb);
    for k=1:length(xcut)
        [tmpx, tmpy] = find(ix_all == round(xcut(k)));
        icut(k) = tmpy;
    end
    close(figure(44))
    
    %--Repeat average
    clear UniqueSingleCycle SingleCycle
    nb_ok = setdiff([1:1:nb],icut);
    TruncateLength = min( lcycle(nb_ok) );
    
    % Extract individual cycles
    j=1;
    close(figure(20)); figure(20);
    hold on
    for k = 1:length(nb_ok)
        SingleCycle = Wave((index(nb_ok(k))-delay):(index(nb_ok(k)+1)-delay));
        
        plot(SingleCycle,'color', [0.5 0.5 0.5])
        
        %truncate to the smaller cycle
        UniqueSingleCycle(:,j) = SingleCycle(1:TruncateLength);
        j=j+1;
    end
    
    %New Ensemble average:
    MeanCycle = mean(UniqueSingleCycle,2);
    SDCycle = std(UniqueSingleCycle,0,2); %standard deviation
    
    
    figure(20)
    set(gca, 'fontsize',14)
    hold on;
    title(['Detection of ', detect_mark])
    plot(MeanCycle,'r', 'linewidth',3)
    disp(['   Ensemble average of ',num2str(length(nb_ok)),' cycles (cycles ', num2str(nb_ok),')']);

    HAPPY = char(questdlg({'Do you need to remove individual cycles from new average?'},'Improve new average?',{'YES'},{'NO'},{'NO'}) );
    
end

end

 