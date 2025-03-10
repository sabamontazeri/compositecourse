load('Transformations')
load('compliance')
syms teta
%%properties of the composite
T300_5208=struct('Ex',181,...
    'Ey',10.3,...
    'Vx',0.28,...
    'Es',7.17);
% Compute Vy dynamically and add it to the structure. It is computed by
% supposing symmetricity.
T300_5208.Vy = T300_5208.Vx * (T300_5208.Ey / T300_5208.Ex);
%%Calculation of on-axis stiffness for T300_5208
S_T300_5208=subs(S,fieldnames(T300_5208), struct2cell(T300_5208));
%%Calculation of Positive stress transformation matrix for 45 degree
TP=subs(T_stress_positive,[teta],deg2rad(45));
%%Calculation of Negative strain transformation matrix for 45 degree
TN=subs(T_strain_negative,[teta],deg2rad(45));
%%Calculation of off-axis compliance for T300_5208
S_off=TN*S_T300_5208*TP,4;
% Create a structure
Constants_off = struct(...
    'E1', double(1/S_off(1,1)), ...
    'v21',double(-S_off(2,1)/S_off(1,1)), ...
    'v61',double(-S_off(3,1)/S_off(1,1)), ...
    'E2', double(1/S_off(2,2)), ...
    'v12',double(-S_off(1,2)/S_off(2,2)), ...
    'v62',double(S_off(3,2)/S_off(2,2)), ...
    'E6', double(1/S_off(3,3)),...
    'v16',double(S_off(1,3)/S_off(3,3)), ...
    'v26',double(S_off(2,3)/S_off(3,3)) ...
);
save('Constants_off.mat','Constants_off')
%%%%%%%%%%%%%%%%%%%
syms Q11 Q12 Q16 Q21 Q22 Q26 Q61 Q62 Q66 m n
Q=[Q11,Q12,Q16;Q12,Q22,Q26;Q16,Q26,Q66];
T_stress_neg=[m^2,n^2,-2*m*n;n^2,m^2,2*m*n;m*n,-m*n,(m^2-n^2)];
T_strain_pos=[m^2,n^2,m*n;n^2,m^2,-m*n;-2*m*n,2*m*n,(m^2-n^2)];
Q_transformed=simplify(T_stress_neg*Q*T_strain_pos);
%%%%%%%%%%%%%%%%%
syms S11 S12 S16 S21 S22 S26 S66 m n
Q=[Q11,Q12,Q16;Q12,Q22,Q26;Q16,Q26,Q66];
T_stress_pos=[m^2,n^2,2*m*n;n^2,m^2,-2*m*n;-m*n,m*n,(m^2-n^2)];
T_strain_neg=[m^2,n^2,-m*n;n^2,m^2,m*n;2*m*n,-2*m*n,(m^2-n^2)];
S_transformed=simplify(T_strain_neg*Q*T_stress_pos)