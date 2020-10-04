function [sum_sq] = c_SumSqrt_sgolay(P,U,dens,K,F)

%Jordi Alastruey, April 2010

%==========================================================================
%	Jorge Mariscal-Harana, King's College London (modifications)
%	v1.0 (12/06/19) - Ensure that P and U are row vectors
%	v1.1 (19/06/19) - Error messages are now using error instead of disp
%
%==========================================================================

%P: Array with pressure points (Pa)
%U: Array with velocity points (m/s)
%dens: blood density (Kg/m3)
%F: Window length of the Savitzky-Golay filter (must be odd)
%K: Order of polynomial fit of the Savitzky-Golay filter (K < F)

% Coefficients of the Savitzky-Golay filter used for the derivatives
[b,g]=sgolay(K,F);

% Data lengths
p=length(P);
u=length(U);

%	Row vectors required for calculations
P = P(:)';
U = U(:)';

if(p>u)
    error('The pressure array is longer than the velocity array');
end
if(p<u)
    error('The velocity array is longer than the pressure array');
end

% First order pressure and velocity derivatives
for n = (F+1)/2:p-(F+1)/2
    dP(n)=g(:,2)'*P(n - (F+1)/2 + 1: n + (F+1)/2 - 1)';
end
for n = (F+1)/2:u-(F+1)/2
    dU(n)=g(:,2)'*U(n - (F+1)/2 + 1: n + (F+1)/2 - 1)';
end

% Sum-of-the-squares
sum_sq=sqrt(sum(dP.^2)/sum(dU.^2))/dens;
return
