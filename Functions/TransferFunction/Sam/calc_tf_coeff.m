%% Function to calculate the coefficients of the transfer function

function [a,b] = calc_tf_coeff(nb_files,na,nb,array_avg_pulse_length,array_cath_avg_pulse,array_cuff_avg_pulse)

a = zeros(1,na+1);
b = zeros(1,nb);

% nb_files = round(nb_files/2);

Y = zeros(sum(array_avg_pulse_length(1,1:nb_files)),1);
A = zeros(sum(array_avg_pulse_length(1,1:nb_files)),na+nb);
B = zeros(na+nb);

for k = 1:nb_files
    for j = 0:array_avg_pulse_length(1,k)-(max(na,nb)+1)
        Y(sum(array_avg_pulse_length(1,1:k))+1 + j - array_avg_pulse_length(1,k)) = array_cath_avg_pulse(k,j + max(na,nb-1)+1);
        %sum(array_avg_pulse_length(1,1:k))+1 + j - array_avg_pulse_length(1,k)
        A(sum(array_avg_pulse_length(1,1:k))+1 + j - array_avg_pulse_length(1,k),:) = [array_cath_avg_pulse(k,j-1 + max(na,nb-1)+1:-1:j-na + max(na,nb-1)+1),array_cuff_avg_pulse(k,j + max(na,nb-1)+1:-1:j-nb+1 + max(na,nb-1)+1)];
    end
end

T = A' * A;
    
for n = 1:na+nb
    for m = 1:na+nb
        B(n,m) = T(n,m);
    end
end
    
D = inv(B);
C = D * A' * Y;

a(1,1) = 1;
a(1,2:end) = -C(1:na,1);
b(1,:) = C(na+1:end,1);
