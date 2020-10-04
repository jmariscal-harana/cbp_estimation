% Create_inflow.m
% Inputs: T = period
%         F = [a0  0 ...
%              a1 b1 ...
%              a2 b2 ...
%              a3 b3 ...
%              a4 b4 ...
%              aN bN] ...
%         OutputFile
%
% No direct outputs

% ====================
% Jorge Aramburu
% KCL, June 2018
% ====================

function Create_inflow(T, F, OutputFile)
[nfil, ncol] = size(F);
NHarm = nfil - 1;

a0 = F(1,1)/2;

for j = 1:NHarm
   H_amp(j) = sqrt(F(j+1,1)^2 + F(j+1,2)^2);
   H_ph(j)  = atan2(F(j+1,1),F(j+1,2));
end

% Save the bcs file!
fid = fopen(OutputFile,'wt');
fprintf(fid,'%d %0.5g %0.5g\n',NHarm, T, a0);
for j = 1:NHarm
    fprintf(fid,'%0.5g %0.5g\n',H_amp(j),H_ph(j));
end

end

