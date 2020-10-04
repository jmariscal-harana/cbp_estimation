%% Function to plot the estimated cath pulse and the actual cath pulse 

function [error_mat] = plot_tf_result(nb_files,best_coeff_a_sbp,best_coeff_b_sbp,best_coeff_a_shape,best_coeff_b_shape,array_avg_pulse_length,array_cath_avg_pulse,array_cuff_avg_pulse,array_sbp,array_mbp,array_dbp)


for k = 1:nb_files  
    %-- Get the estimated sbp
    est_cath_pulse = array_cuff_avg_pulse(k,1:array_avg_pulse_length(k));
    est_cath_pulse = est_cath_pulse * (array_mbp(1,k) - array_dbp(1,k)) / mean(est_cath_pulse) + array_dbp(1,k);    
    est_sbp = max(est_cath_pulse);
 
    
    %-- Wipe from memory est_cath_pulse
    clear est_cath_pulse
    
    %-- Get the correct shape
    est_cath_pulse = filter(best_coeff_b_shape,best_coeff_a_shape,array_cuff_avg_pulse(k,1:array_avg_pulse_length(k)));
    
    %-- Filtering
    [tf_b,tf_a] = butter(5,10/(100/2));
    est_cath_pulse = filtfilt(tf_b,tf_a,est_cath_pulse);
        
    %-- Makes the start and the end of the pulse at the same pressure
    N = length(est_cath_pulse);
    redress_coef = (est_cath_pulse(end)-est_cath_pulse(1))/N;
    est_cath_pulse = est_cath_pulse - est_cath_pulse(1);
    est_cath_pulse = est_cath_pulse - redress_coef*(0:N-1);
    est_cath_pulse = est_cath_pulse / max(est_cath_pulse);
    
    %-- Recalibrate the estimated pulse
%      est_cath_pulse = est_cath_pulse * (est_sbp(1,k) - array_dbp(1,k)) + array_dbp(1,k);
    est_cath_pulse = est_cath_pulse * (array_mbp(1,k) - array_dbp(1,k)) / mean(est_cath_pulse) + array_dbp(1,k);
    
    %-- Recalibrate the measured pulses (cath and cuff)
    cath_pulse = array_cath_avg_pulse(k,1:array_avg_pulse_length(k))
    N = length(cath_pulse);
    redress_coef = (cath_pulse(end)-cath_pulse(1))/N;
    cath_pulse = cath_pulse - cath_pulse(1);
    cath_pulse = cath_pulse - redress_coef*(0:N-1);
    cath_pulse = cath_pulse / max(cath_pulse);
    array_cath_avg_pulse(k,1:array_avg_pulse_length(k)) = cath_pulse * (array_mbp(1,k) - array_dbp(1,k)) / mean(cath_pulse) + array_dbp(1,k);
    cuff_pulse = array_cuff_avg_pulse(k,1:array_avg_pulse_length(k));
    N = length(cuff_pulse);
    redress_coef = (cuff_pulse(end)-cuff_pulse(1))/N;
    cuff_pulse = cuff_pulse - cuff_pulse(1);
    cuff_pulse = cuff_pulse - redress_coef*(0:N-1);
    cuff_pulse = cuff_pulse / max(cuff_pulse);
    array_cuff_avg_pulse(k,1:array_avg_pulse_length(k)) = cuff_pulse * (array_mbp(1,k) - array_dbp(1,k)) / mean(cuff_pulse) + array_dbp(1,k);
    
    %-- Create a table of the central systolic blood pressure estimation serror
%     sbp_error(1,k) = est_sbp(1,k) - array_sbp(1,k);
    error_mat(k,1) = array_sbp(1,k);
    error_mat(k,2) = est_sbp;
    error_mat(k,3) = max(est_cath_pulse);
    
    %-- Plot pulses
    figure
    plot(est_cath_pulse,'b'); hold on;
    plot(array_cath_avg_pulse(k,1:array_avg_pulse_length(k)),'k'); hold on;
    plot(array_cuff_avg_pulse(k,1:array_avg_pulse_length(k)),'k:');
%     plot(est_cath_pulse_sbp,'r'); hold on;
    legend('estimated','Peripheral','Central');
    
    %-- Wipe from memory est_cath_pulse
    clear est_cath_pulse
    clear est_cath_pulse_sbp
end