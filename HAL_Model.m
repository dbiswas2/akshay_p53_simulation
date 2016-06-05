
clear
close all
clc

global k_dph1 K_dph1 k_ph1 K_ph1 k1 K1 p0 delta_0 Vr p1 delta_1 k_Sm k_Spm K_Spm p2 delta_2 k_tm k_S p5 

global delta_5 k_Sw k_Spw K_Spw p6 delta_6 k_tw k_dph2 K_dph2 k_ph2 K_ph2 E ATM_TOT 


k_dph1 = 78;
K_dph1 = 25;
k_ph1 = 3;
K_ph1 = 0.1;
k1 = 10;
K1 = 1.01;
p0 = 0.083;
delta_0 = 0.2;
Vr = 10;
p1 = 0.04;
delta_1 = 0.16;
k_Sm = 0.005;
k_Spm = 1;
K_Spm = 0.1;
p2 = 0.083;
delta_2 = 0.0001;
k_tm = 1;
k_S = 0.015;
p5 = 0.04;
delta_5 = 0.2;
k_Sw =0.003;
k_Spw =1;
K_Spw = 0.1;
p6 = 0.083;
delta_6 = 0.001;
k_tw = 1;
k_dph2 = 96;
K_dph2 = 26;
k_ph2 = 1;
K_ph2 = 0.1;
E = 0.1;
ATM_TOT = 1.3;

initial = zeros(1,14);
[t,x] = ode45(@HAL_ODE,[0 400],zeros(1,14)');

p53_n = x(:,1);
Mdm2_n = x(:,2);
Mdm2_mRNA_n = x(:,3);
p53p_n = x(:,4);
ATMp_n = x(:,5);
Wip1_n = x(:,6);
Wip1_mRNA_n = x(:,7);

p53_c = x(:,8);
Mdm2_c = x(:,9);
Mdm2_mRNA_c = x(:,10);
p53p_c = x(:,11);
ATMp_c = x(:,12);
Wip1_c = x(:,13);
Wip1_mRNA_c = x(:,14);

figure('color','white')
plot(t,p53_n + p53p_n,'LineWidth', 1)
hold on
plot(t,Mdm2_n,'LineWidth',1)
legend('Total p53 = p53 + p53^*','Mdm2','Location','northeast')
xlabel('Time (dimesionless)', 'Fontsize', 12, 'FontWeight','bold')
ylabel('Concentration (dimesionless)', 'Fontsize', 12, 'FontWeight','bold')
title('In Nucleus', 'Fontsize', 14, 'FontWeight','bold')


figure('color','white')
plot(t,p53_c + p53p_c,'LineWidth', 1)
hold on
plot(t,Mdm2_c,'LineWidth',1)
legend('Total p53 = p53 + p53^*','Mdm2','Location','northeast')
xlabel('Time (dimesionless)', 'Fontsize', 12, 'FontWeight','bold')
ylabel('Concentration (dimesionless)', 'Fontsize', 12, 'FontWeight','bold')
title('In Cytoplasm', 'Fontsize', 14, 'FontWeight','bold')
