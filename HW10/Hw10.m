%%%Load the stiffness matrix Q_ij, Transformations
clc
clear all
load('stiffness')
load('Transformations')
syms teta
%%%%%Pre-allocation for an array that stores all plies' directions (angles)
all_plies = [];


%%Take the properties of the composite as inputs
%%%%e.g properties of the T300/5208
Ex = input('Enter Ex (GPa): ');
Ey = input('Enter Ey (GPa): ');
Vx = input('Enter Vx: ');
Es = input('Enter Es (GPa): ');

Properties = struct('Ex', Ex, ...
                   'Ey', Ey, ...
                   'Vx', Vx, ...
                   'Es', Es);
% Compute Vy dynamically and add it to the structure. It is computed by
% supposing symmetricity.
Properties.Vy = Properties.Vx * (Properties.Ey / Properties.Ex);
%%Calculation of on-axis stiffness for MaterialData
Q_on=subs(Q,fieldnames(Properties), struct2cell(Properties));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Calculation of invariants
U1=double(1/8*(3*Q_on(1,1)+3*Q_on(2,2)+2*Q_on(1,2)+4*Q_on(3,3)));
U2=double(0.5*(Q_on(1,1)-Q_on(2,2)));
U3=double(1/8*(Q_on(1,1)+Q_on(2,2)-2*Q_on(1,2)-4*Q_on(3,3)));
U4=double(1/8*(Q_on(1,1)+Q_on(2,2)+6*Q_on(1,2)-4*Q_on(3,3)));
U5=double(1/8*(Q_on(1,1)+Q_on(2,2)-2*Q_on(1,2)+4*Q_on(3,3)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Inputs from
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%user
sy = input('Is the composite symmetric? (yes/no): ', 's');
types=input('enter the number of all kinds of the plies:');
h = input(' Thickness of total layer (mm) = ' ) ;
C= input('c='); %Half thickness of the core used for simplification
n=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Entring stress resultants
fprintf ('\n Enter Stress resultants (MN/m) :\n\n' );
N1 = input(' N1 = ' );
N2 = input(' N2 = ' );
N6 = input(' N6 = ' );
fprintf ('\n Enter Moment resultants (kN) :\n\n' );
M1 = input(' M1 = ' );
M2 = input(' M2 = ' );
M6 = input(' M6 = ' );
M = [M1 ; M2 ; M6] ;
%%%Make a matrix with 2 column, each row is dedicated to a specific ply(Its number and its orientation).
plies=zeros(types,2);
%%%Fill the matrix with the plies' details which are taken as inputs.
for i = 1:types
    % Determine the correct suffix
    if i == 1
        suffix = 'st';
    elseif i == 2
        suffix = 'nd';
    elseif i == 3
        suffix = 'rd';
    else
        suffix = 'th';
    end
    
    % Input number of plies
    plies(i,1) = input(sprintf('Enter the number of the %d%s kind of plies: ', i, suffix));
    
    % Input angle and convert to radians
    plies(i,2) = deg2rad(input(sprintf('Enter the angle of the %d%s kind of plies: ', i, suffix)));
  
    
    % Update total count
    n = n + plies(i,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Make a list with angles of plies n times that are given from user.
%%This is essential for using it in loop
for i = 1:types
    reps = plies(i,1);           % Number of repetitions
    val = plies(i,2);            % Value to repeat
    all_plies = [all_plies; repmat(val, reps, 1)];   
end
%%%%%double the plies in inverse order if the user says the composite is
%%%%%symmetrical
if strcmpi(sy, 'YES')
   all_plies = [all_plies; flipud(all_plies)]
end
%%%They become two times if the symmetry condition is applied
if strcmpi(sy, 'YES')
    n=2*n;
else
    n=n;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%%% Calculating the numbers essential
Nc = 2 * C ; % (Thickness of the core / Thickness of each layer )
h0 = h / (n + Nc) ; % (mm) (Thickness of each layer)
zc = C * h0;      %Thickness of core
I_star=(h)^3/12;  
Zc_star=(2*zc)/h;
h_star = (1-(Zc_star)^3)*I_star ; % (mm)^3
z2 = h/2 ;  %%%% The magnitude of height above the layer
z1 = h/2 - h0 ;  %%%% The magnitude of height under the layer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Preallocation of V_i for A,B,D
V1_A=0;
V2_A=0;
V3_A=0;
V4_A=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V1_B=0;
V2_B=0;
V3_B=0;
V4_B=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
V1_D=0;
V2_D=0;
V3_D=0;
V4_D=0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Calculating V_i for A,B,D for every orientation and adding them together
%%for total V_i calculation
for i = 1:length(all_plies)
    %%%%subtract core thickness if it exists
  if z2 == zc 
    z2=-2*zc+z2
    z1=-2*zc+z1
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  V1_A=V1_A+cos(2*all_plies(i))*(z2 - z1);
  V2_A=V2_A+cos(4*all_plies(i))*(z2 - z1);
  V3_A=V3_A+sin(2*all_plies(i))*(z2 - z1);
  V4_A=V4_A+sin(4*all_plies(i))*(z2- z1);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  V1_B=V1_B+0.5*cos(2*all_plies(i))*(z2^2 - z1^2);
  V2_B=V2_B+0.5*cos(4*all_plies(i))*(z2^2 - z1^2);
  V3_B=V3_B+0.5*sin(2*all_plies(i))*(z2^2 - z1^2);
  V4_B=V4_B+0.5*sin(4*all_plies(i))*(z2^2 - z1^2);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  V1_D=V1_D+((1/(3))*cos(2*all_plies(i))*(z2^3 - z1^3));
  V2_D=V2_D+((1/(3))*cos(4*all_plies(i))*(z2^3 - z1^3));
  V3_D=V3_D+((1/(3))*sin(2*all_plies(i))*(z2^3 - z1^3));
  V4_D=V4_D+((1/(3))*sin(4*all_plies(i))*(z2^3 - z1^3));
    %%%%%%%%%%%%%%%%%%building an array for saving the heights at top and
    %%%%%%%%%%%%%%%%%%bottom of layer
  Z(2*i-1)=z2
  Z(2*i)=z1
    %%%%%%%%%%%%%%%%%%subtract a layer thickness
  z2=z2-h0;
  z1=z1-h0;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
  %%%%%%%%%%%%%%%%%%%Q_matrix
  theta=all_plies(i)
  %%Calculation of Positive strain transformation matrix for 45 degree
  TP=subs(T_strain_positive,[teta],theta);
  %%Calculation of Negative stress transformation matrix for 45 degree
  TN=subs(T_stress_negative,[teta],theta);
  % %Calculation of off-axis stiffness for Properties
  Q_off(:,:,i)=TN*Q_on*TP;
  
  %%%%%%%%%%%%%%%%%%%%%Transformation matrix
  Transformation_strain(:,:,i)=subs(T_strain_positive,[teta],theta);
  Transformation_stress(:,:,i)=subs(T_stress_positive,[teta],theta);
  
end
%%%%%% Add zero element if we have core
if zc~=0
   mid = ceil(length(Z)/2);
   Z = [Z(1:mid), 0, Z(mid+1:end)];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Matrix_A=[U1,V1_A,V2_A;U1,-V1_A,V2_A;U4,0,-V2_A;U5,0,-V2_A;
    0 ,0.5*V3_A, V4_A;0 ,0.5*V3_A ,-V4_A];

Matrix_B=[U1,V1_B,V2_B;U1,-V1_B,V2_B;U4,0,-V2_B;U5,0,-V2_B;
    0 ,0.5*V3_B, V4_B;0 ,0.5*V3_B ,-V4_B];

Matrix_D=[U1, V1_D, V2_D;U1,-V1_D ,V2_D;U4, 0 ,-V2_D;
    U5, 0, -V2_D;0 ,0.5*V3_D, V4_D;0 ,0.5*V3_D ,-V4_D];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
U=[1;U2;U3];
U_A=[h;U2;U3];
U_B=[0;U2;U3];
U_D=[h_star;U2;U3];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%Derivation of A matrix.
A=Matrix_A*U_A;
%%%Derivation of B matrix.
B=Matrix_B*U_B;
%%%Derivation of D matrix.
D=Matrix_D*U_D;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Making a square matrix of A matrix.
A_square=[A(1),A(3),A(5);A(3),A(2),A(6);A(5),A(6),A(4)];

%%%Making a square matrix of DA matrix.
B_square=[B(1),B(3),B(5);B(3),B(2),B(6);B(5),B(6),B(4)];

%%%Making a square matrix of D matrix.
D_square=[D(1),D(3),D(5);D(3),D(2),D(6);D(5),D(6),D(4)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Stiffness=[A_square B_square
    B_square D_square];

%%%Inverse square S matrix to calculate d matrix.
C_matrix=inv(Stiffness)*1000
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_matrix=C_matrix(1:3,1:3);
beta_matrix = C_matrix(1:3,4:6) ;
delta_matrix=C_matrix(4:6,4:6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General strain (epsilon_o) & curvature (k)
strain_vector = C_matrix * [N1 ; N2 ; N6 ; M1 ; M2 ; M6] ;
epsilon_o = strain_vector(1:3) / 10^3 ; % General strain
k = strain_vector(4:6) / 10^3 ; % General curvature (mm)^(-1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf ('\n\n---------------------------------------\n' ) ;

Stiffness
fprintf('\n\nA_square (MN/m) = \n\n ') ;
fprintf('%.3f %.3f %.3f\n %.3f %.3f %.3f\n %.3f %.3f %.3f\n\n' ,A_square);
fprintf('B_square (kN) = \n\n ') ;
fprintf('%.3f %.3f %.3f\n %.3f %.3f %.3f\n %.3f %.3f %.3f\n\n' ,B_square);
fprintf('D_square (N.m) = \n\n ') ;
fprintf('%.3f %.3f %.3f\n %.3f %.3f %.3f\n %.3f %.3f %.3f\n\n' ,D_square);
fprintf ('---------------------------------------\n' ) ;

C_matrix
fprintf('\n\nalpha_square (GN/m)^(-1) = \n\n ' ) ;
fprintf('%.3f %.3f %.3f\n %.3f %.3f %.3f\n %.3f %.3f %.3f\n\n' ,alpha_matrix);
fprintf('beta_square (MN)^(-1) = \n\n ' ) ;
fprintf('%.3f %.3f %.3f\n %.3f %.3f %.3f\n %.3f %.3f %.3f\n\n' ,beta_matrix);
fprintf('delta_square (kN.m)^(-1) = \n\n ' ) ;
fprintf('%.3f %.3f %.3f\n %.3f %.3f %.3f\n %.3f %.3f %.3f\n\n' ,delta_matrix);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R = 1 ;
for i = 1 : length(all_plies)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format short e ;
epsilon_off_up(:,i) = epsilon_o + k * Z(2*i) ;
epsilon_off_down(:,i) = epsilon_o + k * Z(2*i-1) ;
epsilon_off_up(:,i) = round(epsilon_off_up(:,i),6) ;
epsilon_off_down(:,i) = round(epsilon_off_down(:,i),6) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilon_on_up(:,i) = Transformation_strain(:,:,i) * epsilon_off_up(:,i) ;
epsilon_on_down(:,i) = Transformation_strain(:,:,i) * epsilon_off_down(:,i) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format short g ;
sigma_off_up(:,i) = Q_off(:,:,i) * epsilon_off_up(:,i) * 1000 ; % (MPa)
sigma_off_down(:,i) = Q_off(:,:,i) * epsilon_off_down(:,i) * 1000 ; % (MPa)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma_on_up(:,i) = Transformation_stress(:,:,i) * sigma_off_up(:,i) ; % (MPa)
sigma_on_down(:,i) = Transformation_stress(:,:,i) * sigma_off_down(:,i) ; % (MPa)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf ('---------------------------------------' ) ;
fprintf ('\n\n ply %d  ' ,i) ;
fprintf ('\n\nFor the %d-radian ply : \n\n' ,deg2rad(all_plies(i))) ;
fprintf ('sigma_off_up (MPa) = \n\n %.3f\n %.3f\n %.3f\n\n' ,sigma_off_up(:,i)) ;
fprintf ('sigma_off_down (MPa) = \n\n %.3f\n %.3f\n %.3f\n\n' ,sigma_off_down(:,i));
fprintf ('sigma_on_up (MPa) = \n\n %.3f\n %.3f\n %.3f\n\n' ,sigma_on_up(:,i)) ;
fprintf ('sigma_on_down (MPa) = \n\n %.3f\n %.3f\n %.3f\n\n' ,sigma_on_down(:,i)) ;
fprintf ('epsilon_off_up = \n\n %.6f\n %.6f\n %.6f\n\n' ,epsilon_off_up(:,i)) ;
fprintf ('epsilon_off_down = \n\n %.6f\n %.6f\n %.6f\n\n' ,epsilon_off_down(:,i)) ;
fprintf ('epsilon_on_up = \n\n %.6f\n %.6f\n %.6f\n\n' ,epsilon_on_up(:,i)) ;
fprintf ('epsilon_on_down = \n\n %.6f\n %.6f\n %.6f\n\n' ,epsilon_on_down(:,i)) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma_off(:,2*i) = sigma_off_up(:,i) ;
sigma_off(:,2*i-1) = sigma_off_down(:,i) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma_on(:,2*i) = sigma_on_up(:,i) ;
sigma_on(:,2*i-1) = sigma_on_down(:,i) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilon_off(:,2*i) = epsilon_off_up(:,i) ;
epsilon_off(:,2*i-1) = epsilon_off_down(:,i) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilon_on(:,2*i) = epsilon_on_up(:,i) ;
epsilon_on(:,2*i-1) = epsilon_on_down(:,i) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Null = zeros (3,1) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%add a zero element in middle if we have core
if zc ~= 0
    mid = ceil(16/2);
    sigma_off = [sigma_off(:, 1:mid), Null, sigma_off(:, mid+1:end)];
    sigma_on = [sigma_on(:, 1:mid), Null, sigma_on(:, mid+1:end)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma_off_1 = round(double(sigma_off(1,:)), 2);
sigma_off_2 = round(double(sigma_off(2,:)), 2);
sigma_off_6 = round(double(sigma_off(3,:)), 2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma_on_x = round(double(sigma_on(1,:)), 2);
sigma_on_y = round(double(sigma_on(2,:)), 2);
sigma_on_s = round(double(sigma_on(3,:)), 2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%creating labels of sigma magnitudes of plots
Z = round(double(Z), 2);
label1 = cellstr(num2str([sigma_off_1(:) Z(:)], '%.2f , %.3f'));
label2 = cellstr(num2str([sigma_off_2(:) Z(:)], '%.2f , %.3f'));
label6 = cellstr(num2str([sigma_off_6(:) Z(:)], '%.2f , %.3f'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
labelx = cellstr(num2str([sigma_on_x(:) Z(:)] , '%.2f , %.3f'));
labely = cellstr(num2str([sigma_on_y(:) Z(:)] , '%.2f , %.3f'));
labels = cellstr(num2str([sigma_on_s(:) Z(:)] , '%.2f , %.3f'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Off-axis stresses
F1 = figure('Name','sigma-1(off-axis)','NumberTitle','off');
plot(sigma_off_1,Z,'-o')
title('Off-axis Stress Distribution')
xlabel('\sigma1 (MPa)')
ylabel('Z (mm)')
grid on
text(sigma_off_1(:),Z,label1,'VerticalAlignment','bottom',...
'HorizontalAlignment','left')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F2 = figure('Name','sigma-2(off-axis)','NumberTitle','off');
plot(sigma_off_2,Z,'-o')
title('Off-axis Stress Distribution')
xlabel('\sigma2 (MPa)')
ylabel('Z (mm)')
grid on
text(sigma_off_2(:),Z,label2,'VerticalAlignment','bottom',...
'HorizontalAlignment','left')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F3 = figure('Name','sigma-6(off-axis)','NumberTitle','off');
plot(sigma_off_6,Z,'-o')
title('Off-axis Stress Distribution')
xlabel('\sigma6 (MPa)')
ylabel('Z (mm)')
grid on
text(sigma_off_6(:),Z,label6,'VerticalAlignment','bottom',...
'HorizontalAlignment','left')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%On-axis stresses
F4 = figure('Name','sigma-x(on-axis)','NumberTitle','off');
plot(sigma_on_x,Z,'-o')
title('On-axis Stress Distribution')
xlabel('\sigmax (MPa)')
ylabel('Z (mm)')
grid on
text(sigma_on_x(:),Z,labelx,'VerticalAlignment','bottom',...
'HorizontalAlignment','left')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F5 = figure('Name','sigma-y(on-axis)','NumberTitle','off');
plot(sigma_on_y,Z,'-o')
title('On-axis Stress Distribution')
xlabel('\sigmay (MPa)')
ylabel('Z (mm)')
grid on
text(sigma_on_y(:),Z,labely,'VerticalAlignment','bottom',...
'HorizontalAlignment','left')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F6 = figure('Name','sigma-s(on-axis)','NumberTitle','off');
plot(sigma_on_s,Z,'-o')
title('On-axis Stress Distribution')
xlabel('\sigmas (MPa)')
ylabel('Z (mm)')
grid on
text(sigma_on_s(:),Z,labels,'VerticalAlignment','bottom',...
'HorizontalAlignment','left')
