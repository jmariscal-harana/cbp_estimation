%% Central to Radial FT
% 
% Samuel Vennin
% King's College London
% 19 April 2018
% 
%% Load data
clear all, close all, clc
load('Data_Rad_Cen')

%% Prepare dataset
for k=1:size(cBP,2)
    Pcentral = cBP{k};
    Pradial = rBP{k};
    %- Extract landmarks that will be used for rescaling waveforms
    cDBP(k) = min(Pcentral);
    cMBP(k) = mean(Pcentral);
    cSBP(k) = max(Pcentral);
    pDBP(k) = min(Pradial);
    pMBP(k) = mean(Pradial);
    pSBP(k) = max(Pradial);
    %- Normalise waveforms
    Pradial_norm = scaledata(Pradial,0,1);
    Pcentral_norm = scaledata(Pcentral,0,1);
%     figure()
%     plot(Pradial,'r');hold on;plot(Pcentral,'b')
    %- Resample waveforms
    Pcentral_norm = double(Pcentral_norm);
    Pcentral = double(Pcentral);
    Pradial_norm = double(Pradial_norm);
    Pradial = double(Pradial);
    central_norm(k,:) = resample(Pcentral_norm,250,length(Pcentral_norm));
    radial_norm(k,:) = resample(Pradial_norm,250,length(Pradial_norm));
    central(k,:) = resample(Pcentral,250,length(Pcentral));
    radial(k,:) = resample(Pradial,250,length(Pradial));
    clear Pcentral Pradial
    k
end

%% Generate TF
best_rms_error_shape = -1;
best_rms_error_sbp = -1;
best_coeff_a_shape = [];
best_coeff_b_shape = [];
best_coeff_a_sbp= [];
best_coeff_b_sbp = [];
best_error_table_shape = [];
best_error_table_sbp = [];
nb_files = size(central,1);
array_avg_pulse_length = size(central,2)*ones(1,nb_files);

for na = 0:0
    for nb = 1:1
    [coeff_a,coeff_b] = calc_tf_coeff(nb_files,na,nb,array_avg_pulse_length,radial,central);
    [best_coeff_a_sbp,best_coeff_b_sbp,best_rms_error_sbp,best_error_table_sbp,rms_error] = best_tf_coeff(nb_files,coeff_a,coeff_b,best_coeff_a_sbp,best_coeff_b_sbp,best_rms_error_sbp,best_error_table_sbp,array_avg_pulse_length,radial,central,cSBP,cDBP,cMBP,'sbp');
    end
    na
end

for na = 0:40
    for nb = 1:40
    [coeff_a,coeff_b] = calc_tf_coeff(nb_files,na,nb,array_avg_pulse_length,radial,central);
    [best_coeff_a_shape,best_coeff_b_shape,best_rms_error_shape,best_error_table_shape,rms_error] = best_tf_coeff(nb_files,coeff_a,coeff_b,best_coeff_a_sbp,best_coeff_b_sbp,best_rms_error_sbp,best_error_table_sbp,array_avg_pulse_length,radial,central,cSBP,cDBP,cMBP,'shape');
    RMS(na+1,nb) = rms_error;
    end
    na
end

minMatrix = min(RMS(:))
[na,nb] = find(RMS==minMatrix)
[coeff_a,coeff_b] = calc_tf_coeff(nb_files,na,nb,array_avg_pulse_length,radial_norm,central_norm);
[best_coeff_a_shape,best_coeff_b_shape,best_rms_error_shape,best_error_table_shape,rms_error] = best_tf_coeff(nb_files,coeff_a,coeff_b,best_coeff_a_sbp,best_coeff_b_sbp,best_rms_error_sbp,best_error_table_sbp,array_avg_pulse_length,radial,central,cSBP,cDBP,cMBP,'shape');
    
%-- Display the number of coefficients used in the TF
disp('');
disp(['nb coeff a for the shape: ',num2str(length(best_coeff_a_shape))]);
disp(['nb coeff b for the shape: ',num2str(length(best_coeff_b_shape))]);
disp('');
disp(['nb coeff a for the sbp: ',num2str(length(best_coeff_a_sbp))]);
disp(['nb coeff b for the sbp: ',num2str(length(best_coeff_b_sbp))]);
disp('');
        
%-- Trace comparison between the estimated and actual trace
 error_mat = plot_tf_result(nb_files,best_coeff_a_sbp,best_coeff_b_sbp,best_coeff_a_shape,best_coeff_b_shape,array_avg_pulse_length,radial_norm,central_norm,cSBP,cMBP,cDBP);

%% Try transfer function
for k=1:size(central,1)
%     [radial_est(k,:),t] = lsim(hz,central_norm(k,:),t);
%     [c imax] = max(radial_est(k,:));
%     radial_est(k,1:imax) = scaledata(radial_est(k,1:imax),0,1);
%     radial_est(k,imax:end) = scaledata(radial_est(k,imax:end),0,1);
%     radial_est(k,:) = radial_est(k,:)*(MBP(k) - DBP(k))/mean(radial_est(k,:)) + DBP(k);
    [radial_est(k,:)] = radial2central(central_norm(k,:),cMBP(k),cDBP(k),best_coeff_a_shape,best_coeff_b_shape);
    radial(k,:) = radial_norm(k,:)*(pMBP(k)-pDBP(k))/mean(radial_norm(k,:)) + pDBP(k);
    central(k,:) = central_norm(k,:)*(cMBP(k)-cDBP(k))/mean(central_norm(k,:)) + cDBP(k);
    PP(k) = max(radial(k,:)) - min(radial(k,:));
    PP_est(k) = max(radial_est(k,:)) - min(radial_est(k,:));
    r(k) = rmse(radial(k,:),radial_est(k,:));
    figure()
    hold all;plot(t,radial(k,:),'linewidth',2);plot(t,radial_est(k,:),'r','linewidth',2);plot(t,central(k,:),'k','linewidth',2);
    legend('Radial','Estimated','Central')
    k
end

mean(r)
[m s] = bland_altman(PP',PP_est',50,200,-30,30)


