clear; clc; close all;

s = tf('s');

%% ============================================================
%  PART 0: Nominal two-zone thermal model
%  Selected SISO loop: u = Q1  -->  y = T2
%  ============================================================

C1  = 3.0e6;
C2  = 2.5e6;
R1o = 0.010;
R2o = 0.012;
R12 = 0.020;

A = [-(1/(C1*R1o) + 1/(C1*R12)),   1/(C1*R12);
      1/(C2*R12),                 -(1/(C2*R2o) + 1/(C2*R12))];

Bu = [1/C1; 0];

% Disturbances: wT = [Tout; qint1; qint2; Q2]
Bw = [1/(C1*R1o), 1/C1,      0,    0;
      1/(C2*R2o), 0,      1/C2, 1/C2];

C = [0 1];
D = 0;

G  = minreal(ss(A, Bu,     C, 0));      % Q1 -> T2
GT = minreal(ss(A, Bw(:,1), C, 0));     % Tout -> T2

fprintf('\n================ SISO nominal plant ================\n');
disp('G(s) = Q1 -> T2');
G

fprintf('\nGT(s) = Tout -> T2\n');
GT

%% ============================================================
%  PART 1: Final SISO controller from Project Update 1
%  K(s) = Kp * KPI * Klead * Kroll
%  ============================================================

wc = 3e-3;

magG_wc = abs(squeeze(freqresp(G, wc)));
Kp = 1 / magG_wc;

phi_max_deg = 60;
phi_max = deg2rad(phi_max_deg);

alpha = (1 - sin(phi_max)) / (1 + sin(phi_max));
Tl    = 1 / (wc * sqrt(alpha));
Tf    = 1 / (10*wc);

% Keep the same Ti you used before
Ti = 3333.333333;

KPI   = (Ti*s + 1) / (Ti*s);
Klead = (Tl*s + 1) / (alpha*Tl*s + 1);
Kroll = 1 / (Tf*s + 1);

K = minreal(Kp * KPI * Klead * Kroll);

fprintf('\n================ Final controller ================\n');
fprintf('Kp    = %.6e\n', Kp);
fprintf('alpha = %.6f\n', alpha);
fprintf('Tl    = %.6f s\n', Tl);
fprintf('Tf    = %.6f s\n', Tf);
fprintf('Ti    = %.6f s\n', Ti);
disp('K(s) = ');
K

%% ============================================================
%  Classical closed-loop functions
%  ============================================================

L  = minreal(G*K);
S  = minreal(feedback(1, L));      % 1/(1+GK)
T  = minreal(feedback(L, 1));      % GK/(1+GK)
KS = minreal(feedback(K, G));      % K/(1+GK)

fprintf('\n================ Classical loop quantities ================\n');
disp('L(s) = G(s)K(s)'); L
disp('S(s) = 1/(1+G(s)K(s))'); S
disp('T(s) = G(s)K(s)/(1+G(s)K(s))'); T
disp('KS(s) = K(s)/(1+G(s)K(s))'); KS

[GM, PM, Wcg, Wcp] = margin(L);
fprintf('\nMargins:\n');
fprintf('GM   = %.4f (absolute), %.4f dB\n', GM, 20*log10(GM));
fprintf('PM   = %.4f deg\n', PM);
fprintf('Wcg  = %.6e rad/s\n', Wcg);
fprintf('Wcp  = %.6e rad/s\n', Wcp);

%% ============================================================
%  PART 1: Four-block formulation of the SISO control problem
%
%  Exogenous inputs:
%     v = [r; d; n]
%  where
%     r = reference
%     d = outdoor disturbance Tout
%     n = measurement noise
%
%  Plant equations:
%     y = G*u + GT*d
%     e = r - y - n = r - G*u - GT*d - n
%
%  Choose:
%     z  = [e; u; y]
%     yk = e
%
%  Then:
%     [z; yk] = [P11 P12; P21 P22] [v; u]
%  ============================================================

P11 = [ 1,   -GT,  -1;
        0,     0,   0;
        0,    GT,   0 ];

P12 = [ -G;
         1;
         G ];

P21 = [ 1,   -GT,  -1 ];
P22 = -G;

P = minreal(ss([P11 P12;
                P21 P22]));

fprintf('\n================ Four-block generalized plant ================\n');
disp('P11 ='); minreal(P11)
disp('P12 ='); minreal(P12)
disp('P21 ='); minreal(P21)
disp('P22 ='); minreal(P22)

% Closed-loop map from v=[r;d;n] to z=[e;u;y]
Tzv = minreal(lft(P, K));

fprintf('\nClosed-loop generalized transfer Tzv:\n');
Tzv

%% ============================================================
%  PART 2: H-infinity evaluation of the generalized closed loop
%  ============================================================

fprintf('\n================ Generalized closed-loop norm ================\n');
Hinf_Tzv = norm(Tzv, inf);
fprintf('||Tzv||_Hinf = %.8e\n', Hinf_Tzv);

%% ============================================================
%  PART 2: Individual closed-loop channels
%  No shaping filters, only raw channels
%  ============================================================

% Reference to error
Ter_r = minreal(S);

% Disturbance to error
Ter_d = minreal(-S*GT);

% Noise to output
Ty_n = minreal(-T);

% Reference to control
Tu_r = minreal(KS);

% Disturbance to output
Ty_d = minreal(S*GT);

fprintf('\n================ Individual channel H-infinity norms ================\n');
fprintf('||S||_Hinf       = %.8e\n', norm(Ter_r, inf));
fprintf('||S*GT||_Hinf    = %.8e\n', norm(Ter_d, inf));
fprintf('||T||_Hinf       = %.8e\n', norm(Ty_n, inf));
fprintf('||KS||_Hinf      = %.8e\n', norm(Tu_r, inf));

%% ============================================================
%  PART 2: H2 norms only for proper channels
%
%  IMPORTANT:
%  S has direct feedthrough 1, so raw H2 of S is not used.
%  T is strictly proper here, so H2(T) is valid.
%  KS is strictly proper here, so H2(KS) is valid.
%  S*GT is strictly proper, so H2(S*GT) is valid.
%  ============================================================

fprintf('\n================ Proper-channel H2 norms ================\n');

try
    H2_T = norm(T, 2);
    fprintf('||T||_H2        = %.8e\n', H2_T);
catch
    fprintf('||T||_H2        = not available\n');
end

try
    H2_KS = norm(KS, 2);
    fprintf('||KS||_H2       = %.8e\n', H2_KS);
catch
    fprintf('||KS||_H2       = not available\n');
end

try
    H2_SGT = norm(minreal(S*GT), 2);
    fprintf('||S*GT||_H2     = %.8e\n', H2_SGT);
catch
    fprintf('||S*GT||_H2     = not available\n');
end

%% ============================================================
%  Optional plots
%  ============================================================

figure;
margin(L); grid on;
title('Final loop L(s) = G(s)K(s)');

figure;
nyquist(L); grid on;
title('Nyquist of final loop');

figure;
bodemag(S, T, KS, {1e-6, 1e-1}); grid on;
legend('S','T','KS','Location','best');
title('Sensitivity functions');

%% ============================================================
%  PART 3: MIMO system characterization from Assignment 3
%  poles / zeros / minimal realization / balanced realization / HSV
%  ============================================================

Ebat = 13.5 * 3.6e6;

A3 = [A,           [0;0];
      0, 0, 0];

B2_3 = [1/C1, 0,     0;
        0,    1/C2,  0;
        0,    0,   -1/Ebat];

C2_3 = eye(3);
D22_3 = zeros(3,3);

Gyu = minreal(ss(A3, B2_3, C2_3, D22_3));

fprintf('\n================ MIMO characterization ================\n');

disp('Poles of full MIMO u->y plant:');
disp(pole(Gyu));

disp('Transmission zeros of full MIMO u->y plant:');
disp(tzero(Gyu));

Gyu_min = minreal(Gyu, 1e-12);

fprintf('Original order of Gyu     = %d\n', order(Gyu));
fprintf('Order after minreal(Gyu)  = %d\n', order(Gyu_min));

if order(Gyu_min) == order(Gyu)
    fprintf('Result: no removable pole-zero cancellations detected.\n');
else
    fprintf('Result: minreal reduced the model order.\n');
end

%% ============================================================
%  Balanced realization of the stable thermal subsystem only
%  ============================================================

AT = A;
BT = [1/C1, 0;
      0,   1/C2];
CT = eye(2);
DT = zeros(2,2);

sysT = minreal(ss(AT, BT, CT, DT));

Wc = gram(sysT, 'c');
Wo = gram(sysT, 'o');

fprintf('\nThermal controllability Gramian Wc:\n');
disp(Wc);

fprintf('Thermal observability Gramian Wo:\n');
disp(Wo);

eigWc = eig(Wc);
eigWo = eig(Wo);

fprintf('\nEigenvalues of Wc:\n');
disp(sort(eigWc, 'descend'));

fprintf('Eigenvalues of Wo:\n');
disp(sort(eigWo, 'descend'));

[sysb, hsv] = balreal(sysT);

fprintf('\nHankel singular values of thermal subsystem:\n');
disp(hsv);

disp('Balanced realization sysb = ');
sysb

fprintf('\nDone.\n');
