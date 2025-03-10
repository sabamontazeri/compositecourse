load('compliance')
load('stiffness')
%%properties of the composite
T300_5208=struct('Ex',181,...
    'Ey',10.3,...
    'Vx',0.28,...
    'Es',7.17);
% Compute Vy dynamically and add it to the structure. It is computed by
% supposing symmetricity.
T300_5208.Vy = T300_5208.Vx * (T300_5208.Ey / T300_5208.Ex);

S_T300_5208=subs(S,fieldnames(T300_5208), struct2cell(T300_5208));
Q_T300_5208=subs(Q,fieldnames(T300_5208), struct2cell(T300_5208));

%%Show S
disp('Numerical Compliance Matrix (S):');
disp(vpa(S_T300_5208, 4));
%%Show Q
disp('Numerical Stiffness Matrix (Q):');
disp(vpa(Q_T300_5208, 4));
%%define strain
epsilon=[0.001;0.003;0.002]
%%constitutive law
sigma=Q_T300_5208*epsilon,
disp(vpa(sigma, 4));