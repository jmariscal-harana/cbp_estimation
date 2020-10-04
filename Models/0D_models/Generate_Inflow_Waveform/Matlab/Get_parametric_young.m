function Qtot_young = Get_parametric_young(dt, Qi, Qmax, Qmin, ti,tmax,tmin, ts, td, T)


%% Q2:second order, no horizontal slope = 0 at t=0 - YOUNG wave

% [P14 P13 P12 P11 P10     P22 P21 P20    P34 P33 P32 P31 P30]
LHS_Y = [
    tmax^4      tmax^3      tmax^2      tmax^1  1        0 0 0    0 0 0 0 0;  %Q1(tmax) = Qmax
    0           0           0           0       1        0 0 0    0 0 0 0 0;  %Q1(0) = 0
%    0           0           0           1       0        0 0 0    0 0 0 0 0;  %Q1'(0) = 0
    4*tmax^3    3*tmax^2    2*tmax^1    1       0        0 0 0    0 0 0 0 0;  %Q1'(tmax) = 0
    4*3*ti^2    3*2*ti      2           0       0        0 0 0    0 0 0 0 0;  %Q1''(ti) =0
    ti^4        ti^3        ti^2        ti      1        0 0 0    0 0 0 0 0;  %Q1(ti) = Qi  
    0           0           0           0       0          tmax^2   tmax 1       0 0 0 0 0; %Q2(tmax)=Qmax
    0           0           0           0       0          2*tmax   1    0       0 0 0 0 0; %Q2'(tmax)=0
    0           0           0           0       0          ts^2     ts   1       -ts^4      -ts^3     -ts^2    -ts -1; %Q2(ts)=Q3(ts)
    0           0           0           0       0           2*ts     1    0       -4*ts^3    -3*ts^2   -2*ts    -1  0;  %Q2'(ts)=Q3'(ts)  
    0 0 0 0 0    0 0 0     td^4        td^3        td^2    td      1; %Q3(td) = 0
    0 0 0 0 0    0 0 0     4*td^3      3*td^2      2*td    1       0; %Q3'(td)=0
    0 0 0 0 0    0 0 0     tmin^4      tmin^3      tmin^2  tmin    1; %Q3(tmin) = Qmin
    0 0 0 0 0    0 0 0     4*tmin^3    3*tmin^2    2*tmin  1       0  %Q3'(tmin) = 0
    ];

%RHS
RHS_Y = [
    Qmax;
    0;
%     0; %Q1'(0) = 0
    0;
    0;
    Qi;
    Qmax;
    0;
    0;
    0;
    0;
    0;
    Qmin;
    0
];

%Solve system A*X=B: x = A\B;
X = LHS_Y\RHS_Y;

%% Equation of flow curve solution - YOUNG
%Q1, t1
t1=[0:dt:tmax];
Q1 = X(1).*t1.^4 + X(2).*t1.^3 + X(3).*t1.^2 +X(4).*t1 + X(5);

%Q2,t2
t2=[tmax+dt:dt:ts];
Q2 = X(6).*t2.^2 + X(7).*t2 + X(8);

%Q3,t3
t3 = [ts+dt:dt:td];
Q3 = X(9).*t3.^4 + X(10).*t3.^3 + X(11).*t3.^2 + X(12).*t3 + X(13);

%Q4, t4
t4 = [td+dt:dt:T];
Q4 = t4.*0;

Qtot_young = [Q1, Q2, Q3, Q4]; 

