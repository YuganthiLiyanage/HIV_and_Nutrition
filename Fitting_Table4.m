clear all
close all
clc

global  tforward initial_cond ViralData CD4Data ProteinLevel dt




ViralData = [3.5273  5.0143  5.9952  6.2903  5.9324  5.4763...
    5.0858  4.7602  4.5493  4.4706  4.3752  4.4809  4.4172  4.1996  3.8672]';

CD4Data = [897.8974  592.6728  631.8200  727.0780  646.4628  551.5711]'; %cells/\mu l

ProteinLevel = [136.1  88  101.3  106.4  101.5  76.4  86.9...
    91  95.6  91.7  87.5  87  85.1  94.7 116.4]'; %gr/L

initial_cond = [2613 0 1048 69];

tVLdata = [1.9, 5.8, 9.7, 13.8, 17.6, 20.8, 24.7, 27.7,...
   31.7, 40.8, 48.8, 63.3, 94.1, 174.6, 257.4]';
tCD4data = [2.0, 17.8, 32.1, 49.0, 94.0, 259.3]';
tProteindata = [2.2, 6.3, 10.2, 14.0, 18.2, 21.3,...
    25.0, 28.3, 32.2, 41.2, 49.2, 68.0, 93.5, 178.5, 254.1 ]';

dt = 0.1;
tforward = 0:dt:300;

lb = zeros(1,11);



%best fitted values to be used for MCS

k = [82.6351921930605,0.0986268162102787,1.91085674848266e-05,...
    0.907108191942541,10975.2832810487,1.18132063474118,...
    1.68828831603924,7.72709312288070e-11,1.13913660951503e-08,...
    0.0129209172690775,1.23811964165194];

[k,fval] =  fminsearchbnd(@err_in_data,k,lb,[],optimset('Display','iter',...
    'MaxFunEvals', 4000, 'MaxIter', 4000));



[t_r, y_r] = ode23s(@(t,y)Model_HIV_WithinHost(y,k),tforward,initial_cond);

figure(1)
plot(tforward, log10(y_r(:,3)),'-b','LineWidth',2)
hold on 
plot(tVLdata, ViralData, 'ro')
title('Viral Load')

figure(2)
plot(tforward, (y_r(:,3)),'-b','LineWidth',2)
hold on 
plot(tVLdata, 10.^ViralData, 'ro')
title('Viral Load')
figure(3)

plot(tforward, log10(y_r(:,1)),'-b','LineWidth',2)
hold on 
plot(tCD4data, log10(CD4Data), 'ro')
title('CD4 cells')
figure(4)
plot(tforward, (y_r(:,1)),'-b','LineWidth',2)
hold on 
plot(tCD4data, (CD4Data), 'ro')
title('CD4 cells')

figure(5)

plot(tforward, log10(y_r(:,4)),'-b','LineWidth',2)
hold on 
plot(tProteindata, log10(ProteinLevel), 'ro')
title('Total Protein')

fprintf('r = %g\n',  k(1));   
fprintf('d = %g\n', k(2));
fprintf('beta = %g\n', k(3));
fprintf('delta = %g\n',  k(4));
fprintf('pi = %g\n', k(5));
fprintf('c = %g\n',  k(6));
fprintf('c1 = %g\n',  k(7));
fprintf('c2 = %g\n',  k(8));
fprintf('gamma = %g\n',  k(9));
fprintf('mu = %g\n', k(10));
fprintf('lambda = %g\n', k(11));




 function error_in_data = err_in_data(k)
 
 global  tforward initial_cond ViralData CD4Data ProteinLevel 
 
 %initial_cond = [k(12), k(13), k(14), k(15)];
 [~,y] = ode23s(@(t,y)Model_HIV_WithinHost(y,k),tforward,initial_cond);
 
 t_v_measure = [1.9, 5.8, 9.7, 13.8, 17.6, 20.8, 24.7,...
     27.7, 31.7, 40.8, 48.8, 63.3, 94.1, 174.6, 257.4].*10 + 1;
 t_cd4_measure = [2.0, 17.8, 32.1, 49.0, 94.0, 259.3].*10 + 1;
 t_protein_measure = [2.2, 6.3, 10.2, 14.0, 18.2, 21.3,...
     25.0, 28.3, 32.2, 41.2, 49.2, 68.0, 93.5, 178.5, 254.1 ].*10 + 1;
 
 Model_Viral = log10(y(t_v_measure(:),3));
 Model_CD4 = log10(y(t_cd4_measure(:),1));
 Model_Protein = log10(y(t_protein_measure(:),4));
 
 error_in_data = sum((Model_Viral - ViralData).^2) +...
                 sum((Model_CD4 - log10(CD4Data)).^2) +...
                 sum((Model_Protein - log10(ProteinLevel)).^2);           

 end

 function dy = Model_HIV_WithinHost(y,k)

dy = zeros(4,1);

%params = [r d beta  delta pi c c_1 c_2  gamma mu lambda]
r = k(1);
d = k(2);
beta = k(3);
delta = k(4);
pi = k(5);
c = k(6);
c_1 = k(7);
c_2 = k(8);
gamma = k(9);
mu = k(10);
lambda = k(11);



T = y(1);
T_i = y(2);
V = y(3);
P = y(4);



dy(1) = r - beta* V.*T./(1 + c_1*P) - d*T ;
dy(2) = beta* V.*T./(1 + c_1*P)  - delta*T_i;
dy(3) = pi*T_i - c*V - c_2*P.*V;
dy(4)= lambda + gamma*P.*V - mu*P;

 
end