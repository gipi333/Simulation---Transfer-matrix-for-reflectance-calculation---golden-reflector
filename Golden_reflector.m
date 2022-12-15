clear all
clc
close all



%------------------------------------------------------
% Définition des constantes et paramètres 
%------------------------------------------------------
theta = 0

epsilon_sub = 1.5
epsilon_inc = 1
nbr_couches = 8+8*120


for j = 1:2:nbr_couches
 epsilon_j(j) = 1.76^2;
 epsilon_j(j+1) = 1.5^2;
end

    
    for k=0:1:120
        for j = (1+8*(k)):2:(8+8*(k))
            d_j(j) = (0.46)*(160+k)*10^(-9);
            d_j(j+1) = (0.54)*(160+k)*10^(-9);
        end
    end






c = 3*10^8;

for lambda = 300 :900  
    
w = 2* pi * (c/lambda) * 10^(9);


%------------------------------------------------------
% Définition de Q_inc pour le cas TE et TM
%------------------------------------------------------

k_z_inc = sqrt( (w/c)^2 * epsilon_inc - (w/c)^2 * epsilon_inc^2 * sin( theta * (180/pi))^2 );

Qinc_TE = [1,1;k_z_inc,-k_z_inc];
Qinc_TE_inv = inv(Qinc_TE);

Qinc_TM = [1,1;(k_z_inc/epsilon_inc),(-k_z_inc/epsilon_inc)];
Qinc_TM_inv = inv(Qinc_TM);


%------------------------------------------------------
% Définition de R_sub pour le cas TE et TM
%------------------------------------------------------

k_z_sub = sqrt( (w/c)^2 * epsilon_sub - (w/c)^2 * epsilon_inc^2 * sin( theta * (180/pi) )^2    );

R_sub_TE = [1,1;k_z_sub,-k_z_sub];
R_sub_TM = [1,1; (k_z_sub/epsilon_sub) ,(-k_z_sub/epsilon_sub)];





%------------------------------------------------------
% Définition de la matrice de transfert
%------------------------------------------------------

%------------------------------------------------------
% Cas TE
%------------------------------------------------------


T_init_TE = eye(2);

for j=1:nbr_couches
    
    k_z_j(j) = sqrt( (w/c)^2 * epsilon_j(j) - (w/c)^2 * epsilon_inc^2 * sin(theta * (180/pi))^2 );
       
    m_TE = exp( 1i * k_z_j(j) * d_j(j) );
    n_TE = exp( -1i * k_z_j(j) * d_j(j) );
    o_TE = k_z_j(j) * exp( 1i * k_z_j(j) * d_j(j) );
    p_TE = -k_z_j(j) * exp( -1i * k_z_j(j) * d_j(j) );
    Q_TE_j = [m_TE,n_TE;o_TE,p_TE];
    
    R_TE_j = [ 1,1; k_z_j(j), - k_z_j(j) ];
    
    T_j_TE = R_TE_j * inv(Q_TE_j);
    
    T_init_TE = T_init_TE * T_j_TE;
    
end 


T_final_TE = Qinc_TE_inv *  T_init_TE * R_sub_TE;

Transmitance_TE(lambda-299) = (abs(1/T_final_TE(1,1)))^2 * (real(k_z_sub)/ k_z_inc) ;
Reflectance_TE(lambda-299)  = (abs(T_final_TE(2,1)/T_final_TE(1,1)))^2;

%------------------------------------------------------
% Cas TM
%------------------------------------------------------


T_init_TM = eye(2);

for j=1:nbr_couches
       
    k_z_j(j) = sqrt( (w/c)^2 * epsilon_j(j) - (w/c)^2 * epsilon_inc^2 * sin(theta * (180/pi))^2 );
    
    m_TM = exp( 1i * k_z_j(j) * d_j(j) );
    n_TM = exp( -1i * k_z_j(j) * d_j(j) );
    o_TM = k_z_j(j) * exp( 1i * k_z_j(j) * d_j(j) );
    p_TM = -k_z_j(j) * exp( -1i * k_z_j(j) * d_j(j) );
    Q_TM_j = [m_TM,n_TM;(o_TM/epsilon_j(j)),p_TM/epsilon_j(j)];
    
    R_TM_j = [ 1,1; k_z_j(j)/epsilon_j(j), - k_z_j(j)/epsilon_j(j) ];
    
    T_j_TM = R_TM_j * inv(Q_TM_j);
    
    T_init_TM = T_init_TM * T_j_TM;
    
end 


T_final_TM = Qinc_TM_inv * T_init_TM * R_sub_TM;

Transmitance_TM(lambda-299) =  (epsilon_inc/epsilon_sub)   * (abs(1/T_final_TM(1,1)))^2 * (real(k_z_sub)/ k_z_inc) ;
Reflectance_TM(lambda-299)  = (abs(T_final_TM(2,1)/T_final_TM(1,1)))^2;


end



for i=300:900
    lambd(i-299)=i;
    Reflectance(i-299) = (Reflectance_TM(i-299) + Reflectance_TE(i-299))/2;
    Transmitance(i-299) = (Transmitance_TM(i-299) + Transmitance_TE(i-299))/2;
end
  



 plot(lambd,real(Reflectance),'r')
















