clc;
clear all;
%%loading On-axis Stiffness matrix
load('stiffness');
syms teta
m=cos(teta);
n=sin(teta);
%%Positive strain transformation matrix
T_strain_positive=[m^2,n^2,m*n;n^2,m^2,-m*n;-2*m*n,2*m*n,(m^2-n^2)];
%%Negative stress transformation matrix
T_stress_negative=[m^2,n^2,-2*m*n;n^2,m^2,2*m*n;m*n,-m*n,(m^2-n^2)];
%%properties of the composite
T300_5208=struct('Ex',181,...
    'Ey',10.3,...
    'Vx',0.28,...
    'Es',7.17);
% Compute Vy dynamically and add it to the structure. It is computed by
% supposing symmetricity.
T300_5208.Vy = T300_5208.Vx * (T300_5208.Ey / T300_5208.Ex);
%%Calculation of on-axis stiffness for T300_5208
Q_T300_5208=subs(Q,fieldnames(T300_5208), struct2cell(T300_5208));
%%Calculation of Positive strain transformation matrix for 45 degree
TP=subs(T_strain_positive,[teta],deg2rad(45));
%%Calculation of Negative stress transformation matrix for 45 degree
TN=subs(T_stress_negative,[teta],deg2rad(45));
%%Calculation of off-axis stiffness for T300_5208
Q_off=TN*Q_T300_5208*TP;
%%Given off-axis strain
offaxisstrain=[0.001; -0.003; 0.004];
%%Calculation of off-axis stress
Sigma_off=Q_off*offaxisstrain;
disp(vpa(Sigma_off, 3));
%%Saving transformation tensors in case we need it later
T_strain_negative=[m^2,n^2,-m*n;n^2,m^2,m*n;2*m*n,-2*m*n,(m^2-n^2)];
T_stress_positive=[m^2,n^2,2*m*n;n^2,m^2,-2*m*n;-m*n,m*n,(m^2-n^2)];
save('Transformations.mat', 'T_strain_positive', 'T_strain_negative',...
'T_stress_positive', 'T_stress_negative');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('compliance')
load('Transformations')
S_T300_5208=subs(S,fieldnames(T300_5208), struct2cell(T300_5208));
%%Calculation of Positive stress transformation matrix for 45 degree
TP=subs(T_stress_positive,[teta],deg2rad(45));
%%Calculation of Negative strain transformation matrix for 45 degree
TN=subs(T_strain_negative,[teta],deg2rad(45));
%%Calculation of off-axis compliance for T300_5208
S_off=TN*S_T300_5208*TP;
%%Given off-axis strain
offaxisstress=[1; 0; 0];
%%Calculation of off-axis stress
strain_off=S_off*offaxisstress;
disp(vpa(strain_off, 3));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Given off-axis strain
torque=[0; 0; 1];
%%Calculation of off-axis strain for positive 45
strain=S_off*torque;
disp(vpa(strain, 3));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Calculation of Positive stress transformation matrix for 45 degree
TP_minus45=subs(T_stress_positive,[teta],deg2rad(-45));
%%Calculation of Negative strain transformation matrix for 45 degree
TN_minus45=subs(T_strain_negative,[teta],deg2rad(-45));
%%Compliance for -45
S_off_minus45=TN_minus45*S_T300_5208*TP_minus45;
%%strains for negative 45 degree
strain_minus45=S_off_minus45*torque;
disp(vpa(strain_minus45, 3));