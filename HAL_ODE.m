
function dxdt = HAL_ODE(t,x)

global k_dph1 K_dph1 k_ph1 K_ph1 k1 K1 p0 delta_0 Vr p1 delta_1 k_Sm k_Spm K_Spm p2 delta_2 k_tm k_S p5 

global delta_5 k_Sw k_Spw K_Spw p6 delta_6 k_tw k_dph2 K_dph2 k_ph2 K_ph2 E ATM_TOT 
u0 = x(1);
u1 = x(2);
u2 = x(3);
u3 = x(4);
u4 = x(5);
u5 = x(6);
u6 = x(7);

v0 = x(8);
v1 = x(9);
v2 = x(10);
v3 = x(11);
v4 = x(12);
v5 = x(13);
v6 = x(14);


du0 = k_dph1*u5*u3/(K_dph1 + u3) - k1*u1*u0/(K1 + u0) - k_ph1*u4*u0/(K_ph1 + u0) - p0*Vr*(u0 - v0);
du1 = -p1*Vr*(u1 - v1) -delta_1*u1;
du2 = k_Sm + k_Spm*u3^4/(K_Spm^4 + u3^4) - p2*Vr*u2 - delta_2*u2;
du3 = k_ph1*u4*u0/(K_ph1 + u0) - k_dph1*u5*u3/(K_dph1 + u3);
du4 = k_ph2*E*(ATM_TOT - u4)/(K_ph2 + 0.5*(ATM_TOT - u4)) - 2*k_dph2*u5*u4^2/(K_dph2 + u4^2);
du5 = p5*Vr*v5 - delta_5*u5;
du6 = k_Sw + k_Spw*u3^4/(K_Spw^4 + u3^4) - p6*Vr*u6 - delta_6*u6;

dv0 = k_S - k1*v1*v0/(K1 + v0) - p0*(v0-u0) - delta_0*v0;
dv1 = k_tm*v2 - p1*(v1 - u1) - delta_1*v1;
dv2 = p2*u2 - k_tm*v2 - delta_2*v2;
dv3 = 0;
dv4 = 0;
dv5 = k_tw*v6 - p5*v5 - delta_5*v5;
dv6 = p6*u6 - k_tw*v6 - delta_6*v6;

dxdt = [du0 du1 du2 du3 du4 du5 du6 dv0 dv1 dv2 dv3 dv4 dv5 dv6]';