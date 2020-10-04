% Create_inflow_Nektar.m
% Inputs: t = time vector
%         func = function vector
%         NHarm = number of harmonics
%         OutputFile
%
% No direct outputs

% ====================
% Jorge Aramburu
% KCL, June 2018
% ====================

function Create_inflow_Nektar(t, func, NHarm, OutputFile)
    T = t(end); %period
    n_pts = length(func);
    FT = zeros(n_pts,1);
    H_amp = zeros(NHarm,1);
    H_ph = zeros(NHarm,1);

    for k=1:n_pts
        omega=2*pi*(k-1)/n_pts; %omega
        for j=1:n_pts
            FT(k)=FT(k)+func(j)*exp(-i*omega*(j-1))/n_pts;
        end
    end

    for j = 1:NHarm
        H_ph(j)  = atan(-real(FT(j+1))/imag(FT(j+1)));
        H_amp(j) = 2*real(FT(j+1))/sin(H_ph(j));
    end

    % Save the bcs file!
    fid = fopen(OutputFile,'wt');
    fprintf(fid,'%d %0.5g %0.5g\n',NHarm, T, FT(1));
    for j = 1:NHarm
        fprintf(fid,'%0.5g %0.5g\n',H_amp(j),H_ph(j));
    end

end