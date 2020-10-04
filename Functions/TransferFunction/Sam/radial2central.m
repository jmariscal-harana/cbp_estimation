%% Function estimating radial pressure from central pressure through a transfer function
% 
% Samuel Vennin
% King's College London
% 20 April 2018
% 
%%
function [est_radial_pulse] = radial2central(centralP,MBP,DBP,best_coeff_a_shape,best_coeff_b_shape)
    
    %- centralP: input pulse - must be an aortic pressure waveform
    %- MBP: Mean Blood Pressure used for calibration
    %- DBP: Diastolic Blood Pressure used for calibration

    %-- Wipe from memory est_cath_pulse
    clear est_cath_pulse
    
    %- Normalise and resample waveform
    centralP = scaledata(centralP,0,1);
    centralP = resample(centralP,250,length(centralP));
    
    %-- Get the correct shape
    est_radial_pulse = filter(best_coeff_b_shape,best_coeff_a_shape,centralP);
    
    %-- Filtering
    [tf_b,tf_a] = butter(5,10/(100/2));
    est_radial_pulse = filtfilt(tf_b,tf_a,est_radial_pulse);
        
    %-- Makes the start and the end of the pulse at the same pressure
    N = length(est_radial_pulse);
    redress_coef = (est_radial_pulse(end)-est_radial_pulse(1))/N;
    est_radial_pulse = est_radial_pulse - est_radial_pulse(1);
    est_radial_pulse = est_radial_pulse - redress_coef*(0:N-1);
    
    %-- Recalibrate the estimated pulse
    est_radial_pulse = est_radial_pulse / max(est_radial_pulse);
    est_radial_pulse = est_radial_pulse * (MBP - DBP) / mean(est_radial_pulse) + DBP;
    est_radial_pulse = est_radial_pulse';
    %-- Wipe from memory est_cath_pulse
    clear est_cath_pulse
    clear est_cath_pulse_sbp


function dataout = scaledata(datain,minval,maxval)

dataout = datain - min(datain(:));
dataout = (dataout/range(dataout(:)))*(maxval-minval);
dataout = dataout + minval;