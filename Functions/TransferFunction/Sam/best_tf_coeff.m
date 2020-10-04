%% Function to determine which are the best coefficients

function [best_coeff_a,best_coeff_b,best_rms_error,best_error_table,rms_error] = best_tf_coeff(nb_files,coeff_a,coeff_b,best_coeff_a,best_coeff_b,best_rms_error,best_error_table,array_avg_pulse_length,array_cath_avg_pulse,array_cuff_avg_pulse,array_sbp,array_dbp,array_mbp,error_type)

% nb_files = round(nb_files/2);

if strcmp(error_type,'sbp')
    error_table = zeros(1,nb_files);
    error_table_bis = zeros(1,nb_files);
    for k = 1:nb_files
        %-- Get the estimated pulse
        est_cath_pulse = filter(coeff_b,coeff_a,array_cuff_avg_pulse(k,1:array_avg_pulse_length(1,k)));
        
        %-- Calibrate the estimated pulse
        est_cath_pulse = est_cath_pulse * (array_mbp(1,k) - array_dbp(1,k)) / mean(est_cath_pulse);
        
        %-- Calculate the error between estimated and actual central systolic blood pressure
        error_table(1,k) = max(est_cath_pulse) - (array_sbp(1,k) - array_dbp(1,k));
        error_table_bis(1,k) = max(est_cath_pulse) - (array_sbp(1,k) - array_dbp(1,k));
        
        %-- Wipe from memory est_cath_pulse
        clear est_cath_pulse;
    end
elseif strcmp(error_type,'shape')
    error_table = zeros(1,nb_files);
    error_table_bis = zeros(1,nb_files);
    for k = 1:nb_files
        %-- Get the estimated pulse
        est_cath_pulse = filter(coeff_b,coeff_a,array_cuff_avg_pulse(k,1:array_avg_pulse_length(1,k)));
        
        %-- Filtering
        [tf_b,tf_a] = butter(5,10/(100/2));
        est_cath_pulse = filtfilt(tf_b,tf_a,est_cath_pulse);
        
        %-- Makes the start and the end of the pulse at the same pressure
        N = length(est_cath_pulse);
        redress_coef = (est_cath_pulse(end)-est_cath_pulse(1))/N;
        est_cath_pulse = est_cath_pulse - est_cath_pulse(1);
        est_cath_pulse = est_cath_pulse - redress_coef*(0:N-1);
        est_cath_pulse = est_cath_pulse / max(est_cath_pulse);
        
        %-- Makes the start and the end of the pulse at the same pressure
        clear N
        cath_pulse = array_cath_avg_pulse(k,1:round(array_avg_pulse_length(1,k)));
        N = length(cath_pulse);
        redress_coef = (cath_pulse(end)-cath_pulse(1))/N;
        cath_pulse = cath_pulse - cath_pulse(1);
        cath_pulse = cath_pulse - redress_coef*(0:N-1);
        cath_pulse = cath_pulse / max(cath_pulse);
        
        %-- Calculate the fitting rms error on the pulse bit to fit
        pulse_bit = 1;
        error_table(1,k) = norm(est_cath_pulse(1:round(array_avg_pulse_length(1,k)*pulse_bit)) - cath_pulse(1:round(array_avg_pulse_length(1,k)*pulse_bit)))/sqrt(round(array_avg_pulse_length(1,k)*pulse_bit));
        
        est_cath_pulse = est_cath_pulse * (array_mbp(1,k) - array_dbp(1,k)) / mean(est_cath_pulse);
        
        error_table_bis(1,k) = (max(est_cath_pulse) - min(est_cath_pulse)) - (array_sbp(1,k) - array_dbp(1,k));
        
        %-- Wipe from memory est_cath_pulse
        clear est_cath_pulse cath_pulse;
    end
end

rms_error = norm(error_table(1,:),2) / sqrt(nb_files);

if rms_error < best_rms_error || best_rms_error == -1
    best_coeff_a = [];
    cbest_coeff_b = [];
    best_coeff_a = coeff_a;
    best_coeff_b = coeff_b;
    best_rms_error = rms_error;
    best_error_table(1,:) = error_table_bis(1,:);
end
