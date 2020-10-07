%% Physical constants
g = 9.81;                       % gravity acceleration [m/sec^2]
%% Physical parameters
m = 0.023;						% wheel weight [kg]
R = 0.027;						% wheel radius [m]
Jw = m * R^2 / 2;				% wheel inertia moment [kgm^2]
M =  0.556;                     % body weight [kg]
W = 0.105;						% body width [m]
D = 0.1;						% body depth [m]
h = 0.21;       				% body height [m]
L = 0.12;						% distance of the center of mass from the wheel axle [m]
l = 0.12;                      % distance entre les centre des roues 
Jpsi = M * L^2 / 3;				% body pitch inertia moment [kgm^2]
Jphi = M * (W^2 + D^2) / 12;	% body yaw inertia moment [kgm^2]
fm = 0.0022;					% friction coefficient between body & DC motor
fw = 0;           				% friction coefficient between wheel & floor
%% Motors parameters
Jm = 1e-5;						% DC motor inertia moment [kgm^2]
Rm = 6.69;						% DC motor resistance [Om]
Kb = 0.468;						% DC motor back EMF constant [Vsec/rad]
Kt = 0.317;						% DC motor torque constant [Nm/A]
K_PWM = 8.087;                  % Volts to PWM value coefficient [1/V]
Ref = [5;5];

%% Constantes definition
alpha =  Kt / Rm;
beta =  Kt * Kb / Rm + fm;
E_11 = (2 * m + M) * R^2 + 2 * Jw + 2 * Jm;
E_12 = M * L * R - 2  * Jm;
E_22 = M * L^2 + Jpsi + 2 * Jm;


%% Part 1: Model

E=[E_11 E_12;E_12 E_22]; %definition de la matrice E
F = 2*beta*[1 -1;-1 1];  %definition de la matrice F
H = 2*alpha*[1;-1];      %definition de la matrice H
G = [0 0;0 -M*g*L];      %definition de la matrice G
A1 = [zeros(2,2) eye(2,2);-inv(E)*G -inv(E)*F]; % definition de la matrice A1
B1 = [0;0;inv(E)*H];  % definition de la matrice B1
C1 = eye(4);          % definition de la matrice C1 pour avoir tous les etats en sortie
D1 = zeros(4,1);       % definition de la matrice D1

s1 = ss(A1, B1, C1, D1); %cr�ation du mod�le d'�tat du syst�me avec la commmande ss
s1.StateName = {'theta', 'psi', 'theta_dot', 'psi_dot'};
s1.InputName = {'U'};

% ici on a un seul �tat d'�quilibre  celui o� toutes les variables d'�tats valent z�ro
%(X = 0) car l'�quation X_dot = A1*X+B1*U = 0 (avec U=0) => X = 0.

% cet �tat d'�quilibre n'est pas stable car:
% 1) si on positionne le segway debout, et qu'on le maintient ainsi nous
% sommes dans une configuration o� X = 0 c'est-�-dire que nous sommes � l'�tat d'�quilibre
% ensuite si on le lache on constate que le segway tombe(diverge de cet �tat) donc cette
% �quilibre n'est pas stable.
pol = pzmap(s1)
% 2) � l'aide la fonction pzmap on peut visualer que ce syst�me poss�de un
% p�le � partie r�el positive donc il est instable.
%les poles renvoy�s par pzmap sont 0,-375.7201,6.8216 et -6.5067


%% For discrete control and simulation
Ts = 0.0033;                     % Control system sample time
Psi0 = deg2rad(10);             % Initial value to disturb the system
s1z = c2d(s1,Ts,'zoh');         % �chantillonnage du mod�le avec prise en compte de l'effet de blocage du CNA

% v�rification de l'�galit�
Az = expm(s1.a * Ts)            % calcul de exp(A1*Ts)
Phi = s1z.a;                    % extraction de Phi la matrice syst�me du mod�le �chantillonn�
% lorsque l'on excute le code on constate au vu des r�sultats que l'�galit�
% est bien v�rifi�e.

%dans cette relation lorsque sp parours l'axe des imaginaires purs, zp
%d�crit bien le cercle de centre 0 et rayon 1, le cercle unit� esr donc
%bien la limite de stabilit� en discret.

% commandabilit�
co = ctrb(s1z.a, s1z.b); % calcul de la matrice de commandabilit� du syst�me �chantillonn�
rank(co);                % calcul de son rang
% le rang de la matrice de commmandabilit� du syst�me �chantillonn� est 4
% donc ce syst�me est commandable.
% On v�rifie la commandabilit� pour nous assurer que l'on peut effectivement
% appliquer une commande par retour d'�tat sur ce syst�me.

%% Stabilisation par retour d'�tat 
% 1) dans le cas donn� dans cette question o� lambda1 > 0 et lambda2 < 0
% on aura la variable x1k qui va diverger vers l'infini, et la variable x2k
% qui va converger vers 0.

% 2) les p�les dominant en temps discret sont les p�les qui sont proches du
% cercle unit� en �tant � l'int�rieure de celui ci. car en continue les
% p�les dominants sont les p�les les proches de l'axe des imaginaires tout
% en �tant � gauche de celui-ci et si on utilise la relation de passage du
% continue vers le discret les poles dominants seront donc les plus pr�s du
% cercle unit� tout en �tant � l'int�rieur.

% les constantes de temps associ�es � ses p�les(les plus rapides) sont: pour le pole 0.2225
% on a 0.0027s et pour le p�le 0.9743 on a 0.1537s
tau1 = -Ts / log(0.2225)  % constante de temps associ�e au pole stable 0.2225
tau2 = -Ts / log(0.9743)  % constante de temps associ�e au pole stable 0.9743
% formule d'Ackerman pour le calcul du gain de retour d'�tat
K = acker(s1z.a, s1z.b, [0.2225 0.9743 0.987 0.987])

%0.985
%% impl�mentation du retour d'�tat
% round(qui retourne l'entier le plus proche) mod�lise la r�solution des capteurs
% en effet les capteurs li�s au segway ne d�livre que des mesures enti�res.

%% etude du filtre
gamma = 0.999;
Fi = tf([1-gamma 0],[1 -gamma], Ts)
bode(Fi) % diagramme de bode du filtre
fc = bandwidth(Fi)/(2*pi) % calcul la bande passante � -3dB en hz
% c'est un filtre passe bas d'apr�s son digramme de bode
% on a fr�quence de coupure  0.0397 hz
% on peut dire que ce filtre est ad�quat pour isoler le biais(qui un
% signal constant donc un signal basse fr�quence)  car sa fr�quence de coupure est assez faible pour
% isoler la composante continue d'un signal qui lui est mis en entr�e.



%%%%%%%%%%%LQR
A_BAR = [A1, zeros(4, 1);1 ,0,0,0,0];
B_BAR = [[B1/2 B1/2]; 0, 0];
QQ = [
1000, 0, 0, 0, 0
0, 1e03, 0, 0, 0
0, 0, 1e04, 0, 0
0, 0, 0, 1e05, 0
0, 0, 0, 0, 100
];
%QQ = eye(5)*5;
RR = 10000*eye(2);
[KK,S,e] = lqrd(A_BAR, B_BAR, QQ, RR,Ts);
k_f = KK(1, 1:4) % feedback gain
k_i = KK(1, 5) % integral gain
% suppress velocity gain because it fluctuates NXTway-GS
k_f(3) = k_f(3) * 0.85;
%% Controller la vitesse


C = [1 0 0 0];% definition de C pour que theta soit la sortie

N = inv(C*inv(eye(4)+s1z.b*K-s1z.a)*s1z.b); % calcul du pr�filtre

%au vu des r�sultats de l'exp�rience pour l'�vitement d'obstacle, envoyer
% une vitesse oppos�e au moteur des que le robot d�tecte un obstacle, n'est
% pas un choix judicieux car ce changement brusque de consigne de vitesse
% n'est pas bon pour la stabilit� du syst�me.

%pour palier � cette difficult� nous avons pensez � rajouter un
%pr�filtre(qui ici est un filtre passe bas) dans le but d'adoucir le
%changement brusque de consigne qui excite les modes non mo                 d�lis�s du
%syst�me.


%%Follow line
L = [1 -2 4];
 xg=[0 1];

%%%%%%%%%
waypoints=[0 1 1 0;1 1 0 0];