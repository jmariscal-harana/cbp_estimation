function filtered = ApplyFiltering(signal_raw, method, param, show)

% Apply filtering to raw signals (signal_raw) with some noise
% method = 1: Savitzky-golay filter,
%               polynomial order = param(1)
%               nb points = param(2)
% method = 2: low-pass filter
%               f cut = param(1)
%               sampling frequency = param(2)
% param are parameters required for the filter
% show different from 0: makes a plot of filtered signal
%
% Marie Willemet, King's College London, UK
% July 2015 - v1

% Example:
% load test data in /example/PPG_finger_temporal.m
% filtered = ApplyFiltering(signal_raw(:,1),1, [2 83], 1);
% filtered = ApplyFiltering(signal_raw(:,2),2, [250 1000], 1);

%-- Filter
switch method
	case 1
		% Savitzky-Golay filter
		filtered = sgolayfilt(signal_raw, param(1),param(2));
		
	case 2
		warning('FFT filter - ensure heart rate is constant, and signal is periodic!')
		Fs = length(signal_raw);
		samplF = param(2);
		
		Y1 = fft(signal_raw)/samplF;
		f = (samplF/1)*linspace(0,1,Fs)';
		
		Mag = 2 * abs(Y1(1:round(length(Y1)/2+1)));
		Phs = unwrap(angle(Y1(1:round(length(Y1)/2+1))));
		
		time = [0:1:length(signal_raw)-1]./1000;
		
		% plot Magnitude and phase of fourier analysis
		if(show)
			figure(100)
			subplot(2,1,1)
			plot(Mag)
			subplot(2,1,2)
			plot(Phs)
		end
		
		% Apply low-pass filter
		lowPass = param(1);
		rect = zeros(size(Y1));
		rect(1:lowPass+1) = 1; %preserve low positive frequencies
		rect(end-lowPass+1:end) = 1; %perserve low negative frequencies
		filtered = ifft(Y1.*samplF.*rect);
		
% Reconstruct Pressure Waveforms - convolution (BUG: offset!)
%         y1 = Mag(1)/2;
%         for i=2:lowPass
%             y1 = y1 + Mag(i) * cos(f(i)*2*pi.*time + Phs(i));
%         end
		
end

if(show)
	figure
	hold on;
	plot(signal_raw,'-k','LineWidth',2)
	plot(filtered,'--b','LineWidth',2)	
end

end


