clear all
close all
clc

global  tforward initial_cond t_v_measure t_cd4_measure t_p_measure

   
numiter = 1000; 

     
true_params = [80.3611,0.0930935,2.06809e-05,0.931085,11726.4,1.21228,...
               1.94174,8.93301e-09,0.0130875,1.24097];


        
 
X = zeros(length(true_params),numiter); 


t_v_measure =  [1.9, 5.8, 9.7, 13.8, 17.6, 20.8, 24.7, 27.7,...
   31.7, 40.8, 48.8, 63.3, 94.1, 174.6, 257.4].*10+1;
t_cd4_measure = [2.0, 17.8, 32.1, 49.0, 94.0, 259.3].*10+1;
t_p_measure = [2.2, 6.3, 10.2, 14.0, 18.2, 21.3,...
    25.0, 28.3, 32.2, 41.2, 49.2, 68.0, 93.5, 178.5, 254.1 ].*10+1;

dt = 0.1;

tforward = 0:dt:300;


initial_cond = [2613 0 1048 69];


[~, y_trp] = ode23s(@(t,y)Model_HIV_WithinHost(y,true_params),tforward,initial_cond);

 noiselevel = [0.01, 0.05, 0.1, 0.2];

total_ARE =  zeros(length(noiselevel), length(true_params));

%params = [r d beta  delta pi  c c_1 c_2 gamma mu lambda]
 total_ARE_Table = {'r', 'd', 'beta',  'delta', 'pi', 'c', 'c_1', ...
     'gamma', 'mu','lambda'};


for noisei = 1:4
    
rng default
noiselev = noiselevel(noisei)

    for i= 1:numiter
            i

    Model_Viral = log10(y_trp(t_v_measure(:),3));
    Model_CD4 =  log10(y_trp(t_cd4_measure(:),1));
    Model_Protein = log10(y_trp(t_p_measure(:),4));


%  ViralData = (noiselev*Model_Viral.*randn(length(t_v_measure),1))+ Model_Viral;
%  CD4Data = (noiselev*Model_CD4.*randn(length(t_cd4_measure),1))+ Model_CD4;
%  PData   = (noiselev*Model_Protein.*randn(length(t_p_measure),1))+ Model_Protein;

  ViralData = (noiselev*randn(length(t_v_measure),1))+ Model_Viral;
  CD4Data = (noiselev*randn(length(t_cd4_measure),1))+ Model_CD4;
  PData   = (noiselev*randn(length(t_p_measure),1))+ Model_Protein; 
 
            k = true_params;
            lb = zeros(1,10);
             
             k = fminsearchbnd(@(k)err_in_data(k,ViralData, CD4Data,PData),...
                 k,lb,[],optimset('MaxFunEvals', 9000,'MaxIter',9000));
             
             X(:,i) = k';
             
     end
        
        arescore = zeros(1,length(true_params));

    for i = 1:length(true_params)
        arescore(i) = 100*sum(abs(true_params(i) - X(i,:))/abs(true_params(i)))/numiter;
    end
    
    total_ARE(noisei,:) = round(arescore,1);
    total_ARE_Table(noisei+1,:) = num2cell(total_ARE(noisei,:));

end

function error_in_data = err_in_data(k, ViralData, CD4Data,PData) 
 
 global  tforward initial_cond 
 
 [~,y] = ode23s(@(t,y)Model_HIV_WithinHost(y,k),tforward,initial_cond);
 
 t_v_measure =  [1.9, 5.8, 9.7, 13.8, 17.6, 20.8, 24.7, 27.7,...
   31.7, 40.8, 48.8, 63.3, 94.1, 174.6, 257.4].*10+1;
t_cd4_measure = [2.0, 17.8, 32.1, 49.0, 94.0, 259.3].*10+1;
t_p_measure = [2.2, 6.3, 10.2, 14.0, 18.2, 21.3,...
    25.0, 28.3, 32.2, 41.2, 49.2, 68.0, 93.5, 178.5, 254.1 ].*10+1;
 
 Model_Viral_f = log10(y(t_v_measure(:),3));
 Model_CD4_f = log10(y(t_cd4_measure(:),1));
 Model_Protein_f =log10(y(t_p_measure(:),4));
 
 error_in_data = sum((Model_Viral_f - ViralData).^2) +...
                 sum((Model_CD4_f- CD4Data).^2)+...
                 sum((Model_Protein_f - PData).^2);           

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
gamma = k(8);
mu = k(9);
lambda = k(10);


T = y(1);
T_i = y(2);
V = y(3);
P = y(4);


dy(1) = r - beta* V.*T/(1 + c_1*P) - d*T ;
dy(2) = beta* V.*T/(1+c_1*P)  - delta*T_i;
dy(3) = pi*T_i - c*V;
dy(4)= lambda + gamma*P.*V - mu*P;
 
end