function dy = brain_network(t, y)
% BRAIN_NETWORK 
% A system of 17 differential equations modeling the dynamics of interconnected brain regions 
%
% Author:   AmirMohammad Azarmehri


dy = zeros(17, 1);

% indices
CA3_e = 1;
CA3_i = 2;
CA1_e = 3;
CA1_i = 4;
CA3_c = 5;
CA1_c = 6;
CX_e  = 7;
CX_i  = 8;
CX_c  = 9;
CXp_e = 10;
CXp_i = 11;
CXp_c = 12;
TH_t  = 13;
u_t   = 14;
TH_r  = 15;
u_r   = 16;
RE_e  = 17;

% Time constant
tau_e   = 0.02;
tau_i   = 0.01;
tau_c   = 0.50;
tau_t   = 0.04;
tau_r   = 0.02;
tau_u_r = 0.15;
tau_u_t = 0.15;
tau_re  = 0.04;

% Strength of the connections
c_s  = 10;
g_c  = 3;
J_ee = @(c_n) 1./(1+exp((c_n-c_s)/g_c));


% Eq. 1 ------------------------------------------------------------------1

NP_CA3_ee = 1000*0.33;
NP_CA3_ie = 100*0.47;

J_CA3_ee = 1.6*(0.125/0.33) * J_ee(y(CA3_c));
J_CA3_ie = 5.5*0.27/0.47;

r0         = 0.1;
r1         = 70;
Vs_CA3_e   = 26;
g_CA3_e    = 6;
r_CA3_e    = @(V_CA3_e) r0+(r1./(1+exp(-(V_CA3_e-Vs_CA3_e)/g_CA3_e)));

Vs_CA3_i   = 26;
g_CA3_i    = 3;
r_CA3_i    = @(V_CA3_i) r0+(r1./(1+exp(-(V_CA3_i-Vs_CA3_i)/g_CA3_i)));

% Eq. 1
dy(CA3_e) = -y(CA3_e)/tau_e ...
    + NP_CA3_ee*J_CA3_ee*r_CA3_e(y(CA3_e)) ...
    - NP_CA3_ie*J_CA3_ie*r_CA3_i(y(CA3_i));

% Eq. 2 ------------------------------------------------------------------2

NP_CA3_ei = 1000*0.34;
NP_CA3_ii = 100*0.44;

J_CA3_ei = 0.1*0.06/0.34;
J_CA3_ii = 0;

% Eq. 2
dy(CA3_i) = -y(CA3_i)/tau_i ...
    + NP_CA3_ei*J_CA3_ei*r_CA3_e(y(CA3_e)) ...
    - NP_CA3_ii*J_CA3_ii*r_CA3_i(y(CA3_i));

% Eq. 3 ------------------------------------------------------------------3

NP_CA1_ie = 100*0.24;
NP_CA3CA1_ee = 1000*0.21;
NP_RECA1_ee = 1000*0.02;

J_CA1_ie = 1.55*0.7/0.24*2;
J_CA3CA1_ee = 1.5*(0.083/0.21) * J_ee(y(CA1_c));

J_RECA1_ee =  14.*6* J_ee(y(CA1_c));

Vs_CA1_i   = 20;
g_CA1_i    = 1;
r_CA1_i    = @(V_CA1_i) r0+(r1./(1+exp(-(V_CA1_i-Vs_CA1_i)/g_CA1_i)));

Vs_RE_e    = 26;
g_RE_e     = 3;
r_RE_e     = @(V_RE_e) r0+(r1./(1+exp(-(V_RE_e-Vs_RE_e)/g_RE_e)));

% Eq. 3
dy(CA1_e) = -y(CA1_e)/tau_e ...
    - NP_CA1_ie*J_CA1_ie*r_CA1_i(y(CA1_i)) ...
    + NP_CA3CA1_ee*J_CA3CA1_ee*r_CA3_e(y(CA3_e)) ...
    + NP_RECA1_ee*J_RECA1_ee*r_RE_e(y(RE_e));

% Eq. 4 ------------------------------------------------------------------4

NP_CA1_ei = 1000*0.34;
NP_CA3CA1_ei = 1000*0.21;
NP_CA1_ii    = 100*0.26;
NP_RECA1_ei  = 1000*0.02;

J_CA1_ei = 5*0.069/0.34;
J_CA3CA1_ei = 0.06*0.005/0.21;
J_CA1_ii    = 0.2/0.26;
J_RECA1_ei  = 0.01*0.6;

Vs_CA1_e = 19;
g_CA1_e  = 2;
r_CA1_e  = @(V_CA1_e) r0+(r1./(1+exp(-(V_CA1_e-Vs_CA1_e)/g_CA1_e)));

% Eq. 4
dy(CA1_i) = -y(CA1_i)/tau_i ...
    + NP_CA1_ei*J_CA1_ei*r_CA1_e(y(CA1_e)) ...
    + NP_CA3CA1_ei*J_CA3CA1_ei*r_CA3_e(y(CA3_e)) ...
    - NP_CA1_ii*J_CA1_ii*r_CA1_i(y(CA1_i)) ...
    + NP_RECA1_ei*J_RECA1_ei*r_RE_e(y(RE_e));

% Eq. 5 ------------------------------------------------------------------5

Delta_c_CA3  = 0.011*0.05/0.33;

% Eq. 5
dy(CA3_c) = -y(CA3_c)/tau_c ...
    + NP_CA3_ee*Delta_c_CA3*r_CA3_e(y(CA3_e));

% Eq. 6 ------------------------------------------------------------------6

Delta_c_CA1 = 0.027*0.045/0.21;

% Eq. 6
dy(CA1_c) = -y(CA1_c)/tau_c ...
    + NP_CA3CA1_ee*Delta_c_CA1*r_CA3_e(y(CA3_e)) ...
    + NP_RECA1_ee*Delta_c_CA1*r_RE_e(y(RE_e));

% Eq. 7 ------------------------------------------------------------------7

NP_CX_ee     = 0.8*10000*0.2;
NP_CX_ie     = 0.2*10000*0.2;
NP_CX2CX1_ee = 0.8*10000*0.02;
NP_RECX_ee   = 1000*0.5*0.02;
NP_CA1CX_ee  = 1000*0.02;

J_CX_ee     = 0.375*J_ee(y(CX_c))*0.905;
J_CX_ie     = 1.5;
J_CX2CX1_ee = 1.3*1.7*0.19*0.59 * J_ee(y(CX_c));

J_CA1CX_ei  = 0.2;
Delta_c_CXp = 0.0048*1;

J_RECX_ee   = 0.01* J_ee(y(CX_c));
J_CA1CX_ee  = 0.05* J_ee(y(CX_c));

r_CX_e    = @(V_CX_e) 0.1+(70./(1+exp(-(V_CX_e-30)/5)));
r_CX_i    = @(V_CX_i) 0.1+(70./(1+exp(-(V_CX_i-30)/2)));

% Eq. 7
dy(CX_e) = -y(CX_e)/tau_e ...
    + NP_CX_ee*J_CX_ee*r_CX_e(y(CX_e)) ...
    - NP_CX_ie*J_CX_ie*r_CX_i(y(CX_i)) ...
    + NP_CX2CX1_ee*J_CX2CX1_ee*r_CX_e(y(CXp_e)) ...
    + NP_RECX_ee*J_RECX_ee*r_RE_e(y(RE_e)) ...
    + NP_CA1CX_ee*J_CA1CX_ee*r_CA1_e(y(CA1_e));

% Eq. 8 ------------------------------------------------------------------8

NP_CX_ii = 0.2*10000*0.2;
NP_CX_ei = 0.8*10000*0.2;
NP_CX2CX1_ei = 0.8*10000*0.02;
NP_RECX_ei   = 1000*0.5*0.02;
NP_CA1CX_ei  = 1000*0.02;

J_CX_ii = 0.285*0.9;
J_CX_ei = 0.43;
J_CX2CX1_ei = 0.8;
J_RECX_ei   = 0.5*0.05;

% Eq. 8
dy(CX_i) = -y(CX_i)/tau_i ...
    - NP_CX_ii*J_CX_ii*r_CX_i(y(CX_i)) ...
    + NP_CX_ei*J_CX_ei*r_CX_e(y(CX_e)) ...
    + NP_CX2CX1_ei*J_CX2CX1_ei*r_CX_e(y(CXp_e)) ...
    + NP_RECX_ei*J_RECX_ei*r_RE_e(y(RE_e)) ...
    + NP_CA1CX_ei*J_CA1CX_ei*r_CA1_e(y(CA1_e));

% Eq. 9 ------------------------------------------------------------------9

Delta_c_CX = 0.0025*1.1;

% Eq. 9
dy(CX_c) = -y(CX_c)/tau_c ...
    + NP_CX_ee*Delta_c_CX*r_CX_e(y(CX_e)) ...
    + NP_CX2CX1_ee*Delta_c_CX*r_CX_e(y(CXp_e)) ...
    + NP_RECX_ee*Delta_c_CX*r_RE_e(y(RE_e)) ...
    + NP_CA1CX_ee*Delta_c_CX*r_CA1_e(y(CA1_e));

% Eq. 10 ----------------------------------------------------------------10

NP_CX_ie = 0.2*10000*0.2;
NP_CX_ee = 0.8*10000*0.2;
NP_CX1CX2_ee = 0.8*10000*0.02;
NP_THCX_te   = 1000*0.5*0.02;

J_CXp_ie = 1.4;
J_CXp_ee = 1.6*0.3075*J_ee(y(CXp_c))*0.82;
J_CX1CX2_ee =0.5*1.379*J_ee(y(CXp_c));
J_THCX_te   = 1.5*3 * J_ee(y(CXp_c));

R_T_t = 100;
L_t   = 0.04;
V_T_t = 26;
g_T_t = 6;
R_B_t = 250;
V_B_t = 15;
g_B_t = 0.9;
r_TH_t  = @(V_TH_t, u_t) ...
    R_T_t*exp(L_t*u_t)./(1+exp(-(V_TH_t-V_T_t)/g_T_t)) + ...
    R_B_t*(1-exp(L_t*u_t))./(1+exp(-(V_TH_t-V_B_t)/g_B_t));

% Eq. 10
dy(CXp_e) = -y(CXp_e)/tau_e ...
    - NP_CX_ie*J_CXp_ie*r_CX_i(y(CXp_i)) ...
    + NP_CX_ee*J_CXp_ee*r_CX_e(y(CXp_e)) ...
    + NP_CX1CX2_ee*J_CX1CX2_ee*r_CX_e(y(CX_e)) ...
    + NP_THCX_te*J_THCX_te*r_TH_t(y(TH_t), y(u_t));

% Eq. 11 ----------------------------------------------------------------11

NP_CX1CX2_ei = 0.8*10000*0.02;
NP_THCX_ti = 1000*0.5*0.02;

J_CXp_ei = 0.485;
J_CXp_ii = 0.7*0.72;
J_CX1CX2_ei = 0.61;
J_THCX_ti = 0.1*6;

% Eq. 11
dy(CXp_i) = -y(CXp_i)/tau_i ...
    + NP_CX_ei*J_CXp_ei*r_CX_e(y(CXp_e)) ...
    - NP_CX_ii*J_CXp_ii*r_CX_i(y(CXp_i)) ...
    + NP_CX1CX2_ei*J_CX1CX2_ei*r_CX_e(y(CX_e)) ...
    + NP_THCX_ti*J_THCX_ti*r_TH_t(y(TH_t), y(u_t));

% Eq. 12 ----------------------------------------------------------------12

% Eq. 12
dy(CXp_c) = -y(CXp_c)/tau_c ...
    + NP_CX_ee*Delta_c_CXp*r_CX_e(y(CXp_e)) ...
    + NP_CX1CX2_ee*Delta_c_CXp*r_CX_e(y(CX_e)) ...
    + NP_THCX_te*Delta_c_CXp*r_TH_t(y(TH_t), y(u_t));

% Eq. 13 ----------------------------------------------------------------13

NP_TH_rt = 1000*0.5*0.01;
NP_CXTH_et = 10000*0.8*0.2;

J_TH_rt = 1.65*2.2*4.6;
J_CXTH_et = 0.2*2;

f_max_t = 250;
f_th_t  = 50;
q_t     = 0.005;
f_u_t = @(u_t) f_max_t./(1+exp((u_t+f_th_t)/q_t));
A_t = 0.2;

R_T_r = 100;
L_r   = 0.04;
V_T_r = 26;
g_T_r = 6;
R_B_r = 250;
V_B_r = 15;
g_B_r = 0.9;
r_TH_r  = @(V_TH_r, u_r) ...
    R_T_r*exp(L_r*u_r)./(1+exp(-(V_TH_r-V_T_r)/g_T_r)) + ...
    R_B_r*(1-exp(L_r*u_r))./(1+exp(-(V_TH_r-V_B_r)/g_B_r));

% Eq. 13
dy(TH_t) = 0.9*( ...
    -y(TH_t)/tau_t ...
    + f_u_t(u_t)./A_t ...
    - NP_TH_rt*J_TH_rt*r_TH_r(y(TH_r), y(u_r)) ...
    + NP_CXTH_et*J_CXTH_et*r_CX_e(y(CXp_e)) ...
    );

% Eq. 14 ----------------------------------------------------------------14

if y(TH_t) >= -0.1, b_t = 0; else, b_t = -200; end % mA

% Eq. 14
dy(u_t) = 0.9*((b_t-y(u_t))/tau_u_t);

% Eq. 15 ----------------------------------------------------------------15

NP_TH_rr = 1000*0.5*0.001;
NP_TH_tr = 1000*0.5*0.02;
NP_CXTH_er = 10000*0.8*0.2;

J_TH_rr = 1.*1.5;
J_TH_tr = 1.2*1.5*0.9;
J_CXTH_er = 0.43*0.5;

f_max_r = 250;
f_th_r  = 50;
q_r     = 0.005;
f_u_r = @(u_r) f_max_r./(1+exp((u_r+f_th_r)/q_r));
A_r = 0.2;



% Eq. 15
dy(TH_r) = 0.9*( ...
    -y(TH_r)/tau_r ...
    + f_u_r(u_r)/A_r ...
    - NP_TH_rr*J_TH_rr*r_TH_r(y(TH_r), y(u_r)) ...
    + NP_TH_tr*J_TH_tr*r_TH_t(y(TH_t), y(u_t)) ...
    + NP_CXTH_er*J_CXTH_er*r_CX_e(y(CXp_e)) ...
    );

% Eq. 16 ----------------------------------------------------------------16

if y(TH_r) >= 0, b_r = 0; else, b_r = -200; end % mA

% Eq. 16
dy(u_r) = 0.9*( ...
    (b_r-y(u_r))/tau_u_r ...
    );

% Eq. 17 ----------------------------------------------------------------17

NP_THRE_rre  = 1000*0.5*0.02;
NP_THRE_tre  = 1000*0.5*0.02;
NP_CXRE_ere  = 0.8*10000*0.2;
NP_CA1RE_ere = 100*0.1;

J_THRE_rre   = 1.4;
J_THRE_tre   = 2.2*0;
J_CXRE_ere   = 0.35;
J_CA1RE_ere  = 2;

% Eq. 17
dy(RE_e) = 0.9*( ...
    -y(RE_e)/tau_re ...
    - NP_THRE_rre*J_THRE_rre*r_TH_r(y(TH_r), y(u_r)) ...
    + NP_THRE_tre*J_THRE_tre*r_TH_t(y(TH_t), y(u_t)) ...
    + NP_CXRE_ere*J_CXRE_ere*r_CX_e(y(CX_e)) ...
    + NP_CA1RE_ere*J_CA1RE_ere*r_CA1_e(y(CA1_e)) ...
    );

end
