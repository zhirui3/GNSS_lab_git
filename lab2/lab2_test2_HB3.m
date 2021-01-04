% lab2 test test 
close all 
clear all

%% Define all constants 
c = 299792458; % vacuum speed of light [m/s]
omega_E = 7.2921151467e-5; % Earth rotation rate [rad/s]

    
% KIT
xr_k = 4239065.798;
yr_k = 828841.745;
zr_k = 4678313.623;

%% 2.1 Pseudorange solution 
%% read satellite position and observations 
GARM_C1 = load('GARM_C1.obs');       % C1 observation
GARM_P2 = load('GARM_P2.obs');       % P2 observation
GARM_cor = load('GARM.crd');         % station coordinates 
 
% satellite coordinates 
Sat_x_cor = load('SAT_X.coo');
Sat_y_cor = load('SAT_Y.coo');
Sat_z_cor = load('SAT_Z.coo');
% satellite velocities
Sat_x_vel = load('SAT_X.vel');
Sat_y_vel = load('SAT_Y.vel');
Sat_z_vel = load('SAT_Z.vel');
% satellite clock correction 
Sat_clk = load('SAT_T.clk');
% epochs 
epoch = load('EPOCHS.lst');

% biases 
P1C1 = load('SATP1C1.dcb');
P1P2 = load('SATP1P2.dcb');
%% loop over epochs and apply corrections 
n = length(epoch);

% define other parameters in the loop 
T0 = 2.3;
omega = [0; 0; omega_E];
f0 = 10.23;
f1 = 154*f0;
f2 = 120*f0;
% initialize matrix for storing results 
del_p_mat = zeros(120,4);
sigma_p_mat = zeros(120,4);

pos_x=zeros(n,1);
pos_y=zeros(n,1);
pos_z=zeros(n,1);

% compute coordinate transformation 
% computing longitude and latitude of the receiver 
%lat = atan(zr_G/6371000);
%long = 2*atan(yr_G/(sqrt(xr_G^2 + yr_G^2) + xr_G));
%R = [-sin(lat)*cos(long), -sin(lat)*sin(long), cos(lat); -sin(long), cos(long),0;cos(lat)*cos(long), cos(lat)*sin(long), sin(lat)];

% define max number of iterations 
n_iter = 2;


%switches for different corrections. For all apart from 'clock', keep as
%'y' to include correction or change to anything else to remove (i.e 'n')
% for the clock option choose 'rel' for the relativistic correction,
% 'non_rel' for standard clock correction, or any other letter (i.e 'n')
% for no clock correction at all 

ion='n';
biases='n';
tropo='n';
rot='n';
light='n';
clock='c';

title1="Using";

if ion=='y'
    title1=title1+" Ionospheric Correction,";
end

if biases=='y'
    title1=title1+" Code Biases,";
end

if tropo=='y'
    title1=title1+" Tropospheric Correction,";
end

if rot=='y'
    title1=title1+" Sagnac Correction,";
end

if light=='y'
    title1=title1+" Light Travel Correction,";
end

if clock=='r'
    title1=title1+" and Relativistic Correction";
elseif clock=='c'
    title1=title1+" and Non-Relativistic Clock Correction";
else 
    title1=title1+" and No Clock Correction";
end
    
title2=convertStringsToChars(title1);


for i = 1:n 
    % initialize del p and sigma p vectors 
    del_p = zeros(1,4);
    sigma_p = zeros(1,4);

    % define receiver clock correction 
    del_t_r = 0;
    rece_corr = 0;
    % initilize iteration 
    iter = 1; 
    
    % priori coordinates [m] 
    % Garmisch  
    xr_G = 4235957.17962; 
    yr_G = 834342.53713;
    zr_G = 4681540.90061;
    
    while  (1)
        
        %fprintf('Iteration %d \n',iter);
        
        % update initial condition
        
        %rece_corr = del_t_r + del_p(4);
       
        
 
        % difference between raw satellite coordinates and receiver 
        diff_sat_x = Sat_x_cor(i,:) - xr_G ;
        diff_sat_y = Sat_y_cor(i,:) - yr_G ;
        diff_sat_z = Sat_z_cor(i,:) - zr_G ;
        
        p_raw = sqrt(diff_sat_x.^2 + diff_sat_y.^2 + diff_sat_z.^2);
        
        % 2. light travel correction 
        del_t = p_raw./c;
        
        
        % calculating satellite radial velocity 
        sate_vel_rad = zeros(1,7);
        
        %sate_vel_radx = diff_sat_x.*Sat_x_vel(i,:);
        %sate_vel_rady = diff_sat_y.*Sat_y_vel(i,:);
        %sate_vel_radz = diff_sat_z.*Sat_z_vel(i,:);
        
        dot1 = [diff_sat_x; diff_sat_y; diff_sat_z];
        dot2 = [Sat_x_vel(i,:); Sat_y_vel(i,:); Sat_z_vel(i,:)]; 
        for j = 1:7
            sate_vel_rad(1,j) = dot(dot1(:,j),dot2(:,j))/p_raw(1,j);
        end
        
        
        % 3. Sagac correction 
        del_OmegaE = del_t*omega_E;

        % correct satellite coordinates
        Sat_x_corr = Sat_x_cor(i,:).*cos(del_OmegaE) + Sat_y_cor(i,:).*sin(del_OmegaE);
        Sat_y_corr = -1.*Sat_x_cor(i,:).*sin(del_OmegaE) + Sat_y_cor(i,:).*cos(del_OmegaE);
        Sat_z_corr = Sat_z_cor(i,:);
        
        % recompute difference 
        diff_sat_xx = Sat_x_corr - xr_G;
        diff_sat_yy = Sat_y_corr - yr_G;
        diff_sat_zz = Sat_z_corr - zr_G;

        % recompute the geometric distance and previous corrections before
        p_corr0 = sqrt(diff_sat_xx.^2 + diff_sat_yy.^2 + diff_sat_zz.^2);
        
        % 4. relativisic correction 
        %define rotation vector 
        cross_sate_omega = zeros(3,7);
        
        % satellite velocity matrices 
        % 3 x 7 (3 velocity vectors and 7 satellites)
        sate_vel_mat = [Sat_x_vel(i,:); Sat_y_vel(i,:); Sat_z_vel(i,:)];
    
        % satellite position matrices ( not too sure if you are using corrected
        % or raw satellite coordinates)
        r_s_mat = [Sat_x_cor(i,:); Sat_y_cor(i,:); Sat_z_cor(i,:)];
        
        % compute the cross product between the earth rotation vector and
        % satellite position vector 
        for j = 1:7
            cross_sate_omega(:,j) = cross(omega,r_s_mat(:,j));
        end
        
        % compute r prime 
        r_s_prime_mat = sate_vel_mat + cross_sate_omega;
    
        % compute dot product between r_s and r_prime
        r_dot = dot(-2*r_s_mat,r_s_prime_mat);
    
        % compute Delta t 
        Del_t = r_dot/(c^2);
        
        % now compute relativistic correction 
        sate_rel_corr = Sat_clk(i,:) + Del_t;
        
        % 5. Tropospheric correction 
        
        e_G=[xr_G; yr_G; zr_G]./sqrt(xr_G^2+yr_G^2+zr_G^2);
        
        if rot =='y'
            e_Sat = [diff_sat_xx; diff_sat_yy; diff_sat_zz]./ sqrt(diff_sat_xx.^2 + diff_sat_yy.^2 + diff_sat_zz.^2);
        else 
            e_Sat = [diff_sat_x; diff_sat_y; diff_sat_z]./ sqrt(diff_sat_x.^2 + diff_sat_y.^2 + diff_sat_z.^2);
        end
        
        cosz = zeros(1,7);
        for j = 1:7
            cosz(1,j) = dot(e_G, e_Sat(:,j));
        end
        
        T = T0./cosz;
        
        % 6. Ionospheric delay 
        PC = (f1^2*GARM_C1(i,:) - f2^2*GARM_P2(i,:))./(f1^2 - f2^2);
        
        % 7. bias 
        bias = -f2^2/(f1^2 - f2^2)*P1P2 - P1C1;
        PC_bias = - f1^2/(f1^2 - f2^2)*P1C1;
        
        %use/don't use sagnac correction
        
        if rot=='y'
            p = p_corr0+c*rece_corr;
        else 
            p=p_raw+c*rece_corr;
        end
        
        
        %biases
            
        if biases =='y' && ion == 'y'
            p=p + PC_bias;
        elseif biases == 'y' && ion == 'n'
            p = p + bias;
        elseif biases == 'n'
           
        end
            
        %troposphere correction
        
        if tropo=='y'
            p=p+T;
        end 
        
        %relativistic correction, clock correction, or no clock correction
        
        if clock=='r'
            p=p-c.*sate_rel_corr;
            
        elseif clock=='c'
            p=p-c.*Sat_clk(i,:);
            
        else 
            %do nothing
        end
        
        %light travel correction
        
        if light=='y'
            p=p-sate_vel_rad.*del_t;
        end
        
          
        %ionospheric correction
        
        if ion=='y'
            lhs = PC- p;
        else 
            lhs = GARM_C1(i,:) - p;
        end
        
        
          
        % compute unit vectors 
    
        if rot=='y'
            e_x = diff_sat_xx./p_corr0;
            e_y = diff_sat_yy./p_corr0;
            e_z = diff_sat_zz./p_corr0;
        else 
            e_x = diff_sat_x./p_raw;
            e_y = diff_sat_y./p_raw;
            e_z = diff_sat_z./p_raw;
        end
            

        % initialize grand design matrix 
        A = [-1*e_x', -1*e_y', -1*e_z', c*ones(7,1)];
        del_p = (A'*A)\(A'*lhs');
        
        % compute cofactor matrix 
        Q = inv(A'*A); 
        
        % computer residual 
        r_p = lhs' - A*del_p;
        
        % compute posteriroir RMS of the weight
        m0_p = sqrt((r_p'*r_p)/3);
        
        for j = 1:4 
            sigma_p(1,j) = m0_p*sqrt(Q(j,j));
        end
        
            rece_corr = del_t_r + del_p(4);
            
            xr_G = xr_G + del_p(1);
        
       
            yr_G = yr_G + del_p(2);
        
        
        
            zr_G = zr_G + del_p(3);
          
            fprintf('Iteration: %d \n',iter);
            fprintf('X: %f \n',xr_G);
            fprintf('Y: %f \n',yr_G);
            fprintf('Z: %f \n',zr_G);
            fprintf('t: %f \n',rece_corr);
         
            fprintf('Correction x:%f\n', del_p(1)); 
            fprintf('Error x:%f\n', sigma_p(1)); 
            fprintf('Correction y:%f\n', del_p(2)); 
            fprintf('Error y:%f\n', sigma_p(2)); 
            fprintf('Correction z:%f\n', del_p(3)); 
            fprintf('Error z:%f\n', sigma_p(3)); 
          
            iter = iter + 1;  
            
            if(((abs(del_p(1)) <= sigma_p(1)/1000000 && abs(del_p(2)) <= sigma_p(2)/1000000 && abs(del_p(3)/1000000 <= sigma_p(3))) && abs(del_p(4)/1000000 <= sigma_p(4)) || iter>=n_iter))
                break; 
            end
        
    end
    
    %long = atan(zr_G/6371000);
    %lat = 2*atan(yr_G/(sqrt(xr_G^2 + yr_G^2) + xr_G));
    %R = [-sin(lat)*cos(long), -sin(lat)*sin(long), cos(lat); -sin(long), cos(long),0;cos(lat)*cos(long), cos(lat)*sin(long), sin(lat)];
    
    long = 11.083379;
    lat = 47.303331;
    R = [-sind(lat)*cosd(long), -sind(lat)*sind(long), cosd(lat); -sind(long), cosd(long),0; cosd(lat)*cosd(long), cosd(lat)*sind(long), sind(lat)];
    % Crop out the time part, and do transformation for station
    % coordinates
    A_new = A(:,1:3);
    Q_new = inv(A_new'*A_new);
    Q_new_trans = R*Q_new*R';
    
    % Transforming station coordiantes
    del_p_sta = R*del_p(1:3);
    
    
    for j = 1:3
        sigma_p(1,j) = m0_p*sqrt(Q_new_trans(j,j));
    end
    % Storing results 
    del_p_mat(i,:) = [del_p_sta', del_p(4)];
    %del_p_mat(i,:) = del_p;
    sigma_p_mat(i,:) = sigma_p;
    
    pos_x(i)=xr_G;
    pos_y(i)=yr_G;
    pos_z(i)=zr_G;
    
end 


% plot 
figure (1) 
errorbar(epoch./3600 ,del_p_mat(:,1),sigma_p_mat(:,1),'.');
hold on 
errorbar(epoch./3600 ,del_p_mat(:,2),sigma_p_mat(:,2),'.');
errorbar(epoch./3600 ,del_p_mat(:,3),sigma_p_mat(:,3),'.');
xlabel('Epoch [hours]')
ylabel('Offset [m]')
legend('N','E','U')
title(title2,'FontSize',10)

figure(2) 
errorbar(epoch./3600 ,del_p_mat(:,4),sigma_p_mat(:,4),'.');
xlabel('Epoch [hours]')
ylabel('Receiver clock correction [s]')
title(title2,'FontSize',10)

figure(3) 
errorbar(epoch./3600 ,pos_x,sigma_p_mat(:,1),'.');
hold on 
errorbar(epoch./3600 ,pos_y,sigma_p_mat(:,2),'.');
errorbar(epoch./3600 ,pos_z,sigma_p_mat(:,3),'.');
xlabel('Epoch [hours]')
ylabel('Position [m]')
legend('N','E','U')
title(title2,'FontSize',10)