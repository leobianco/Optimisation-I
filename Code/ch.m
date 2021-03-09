%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          PROJET OPT201                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Etudiants : LAMBERT Guillaume et MARTINS BIANCO Leonardo

%%%%%%%%%%%%%% VARIABLE PILOTANT LE COMPORTEMENT DU SIMULATEUR %%%%%%%%%%%%%%%%

indic = 1;

% indic = 1 : chs dessine la cha�ne
% indic = 2 : chs calcule e et c
% indic = 4 : chs calcule e, c, g, a
% indic = 5 : chs calcule hl
% indic = 6 : chs affiche l'�tude des erreurs du gradient de l'�nergie e

%%%%%%%%%%%%%%%%%%%%%%%% STRUCTURE OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.tol = [10**(-8),10**(-8)]; % On arr�te l'algorithme quand |grad| < tol
options.maxit = 500; % Nombre maximal d'it�rations.
options.rl = 0; % 0 pour recherche lin�aire, 1 pour l'algorithme de Newton.
options.verb = 1; % 1 pour afficher les r�sultats, 2 pour d�tails de la RL.

%%%%%%%%%%%%%%%%%%%%%%%%%%%% TESTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Pour utiliser un cas-test, il faut juste effacer %{ au-dessus de son code.
%
%%%%%%%%%%%%%%%%%%%%%%% CAS TEST 2a %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{

% Coordonn�es du dernier point
global A; A = 1;
global B; B = -1;

% Tailles des barres
global L; L = [0.7; 0.5; 0.3; 0.2; 0.5];
global Nb; Nb = length(L);

% Tableau des coordonn�es des points int�rieurs et sa taille
xy = [0.2; 0.4; 0.6; 0.8; -1.0; -1.5; -1.5; -1.3];
global N; N = length(xy);

% Nombre de points int�rieurs de la cha�ne
global Nn; 
Nn = N/2;

% Multiplicateurs de Lagrange (lm) d�finis arbitrairement
lm = [0.5077,0.4223,0.5190,0.6156,0.8774]';

[e, c, g, a, hl, indic] = chs(1, xy, lm); %�tat initial

[x, lm, info] = sqp('chs', xy, lm, options);

[e, c, g, a, hl, indic] = chs(1, x, lm); %�tat final

%}

%%%%%%%%%%%%%%%%%%%%%%% CAS TEST 2b %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{

% Coordonn�es du dernier point
global A; A = 1;
global B; B = -1;

% Tailles des barres
global L; L = [0.7; 0.5; 0.3; 0.2; 0.5];
global Nb; Nb = length(L);

% Tableau des coordonn�es des points int�rieurs et sa taille
xy = [0.2; 0.4; 0.6; 0.8; 1.0; 1.5; 1.5; 1.3];
global N; N = length(xy);

% Nombre de points int�rieurs de la cha�ne
global Nn; 
Nn = N/2;

% Multiplicateurs de Lagrange (lm) d�finis arbitrairement
lm = [0.5077,0.4223,0.5190,0.6156,0.8774]';

[e, c, g, a, hl, indic] = chs(1, xy, lm); %�tat initial

[x, lm, info] = sqp('chs', xy, lm, options);

[e, c, g, a, hl, indic] = chs(1, x, lm); %�tat final

%}

%%%%%%%%%%%%%%%%%%%%%%% CAS TEST 2c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{

% Coordonn�es du dernier point
global A; A = 1;
global B; B = -1;

% Tailles des barres
global L; L = [0.7; 0.5; 0.3; 0.2; 0.5];
global Nb; Nb = length(L);

% Tableau des coordonn�es des points int�rieurs et sa taille
xy = [0.2; 0.4; 0.6; 0.8; -1.0; -1.5; 1.5; -1.3];
global N; N = length(xy);

% Nombre de points int�rieurs de la cha�ne
global Nn; 
Nn = N/2;

% Multiplicateurs de Lagrange (lm) d�finis arbitrairement
lm = [0.5077,0.4223,0.5190,0.6156,0.8774]';

[e, c, g, a, hl, indic] = chs(1, xy, lm); %�tat initial

[x, lm, info] = sqp('chs', xy, lm, options);

[e, c, g, a, hl, indic] = chs(1, x, lm); %�tat final

%}

%%%%%%%%%%%%%%%%%%%%%%% CAS TEST 2d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{

% Coordonn�es du dernier point
global A; A = 1;
global B; B = -1;

% Tailles des barres
global L; L = [0.7; 0.5; 0.3; 0.2; 0.5];
global Nb; Nb = length(L);

% Tableau des coordonn�es des points int�rieurs et sa taille
xy = [0.2; 0.4; 0.6; 0.8; 1.0; -1.2; 1.5; -1.3];
global N; N = length(xy);

% Nombre de points int�rieurs de la cha�ne
global Nn; 
Nn = N/2;

% Multiplicateurs de Lagrange (lm) d�finis arbitrairement
lm = [0.5077,0.4223,0.5190,0.6156,0.8774]';

[e, c, g, a, hl, indic] = chs(1, xy, lm); %�tat initial

[x, lm, info] = sqp('chs', xy, lm, options);

[e, c, g, a, hl, indic] = chs(1, x, lm); %�tat final

%}

%%%%%%%%%%%%%%%%%%%%%%% CAS TEST 3a %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{

% Coordonn�es du dernier point
global A; A = 1;
global B; B = 0;

% Tailles des barres
global L; L = [0.6; 0.6];
global Nb; Nb = length(L);

% Tableau des coordonn�es des points int�rieurs et sa taille
xy = [0.5;0.4];
global N; N = length(xy);

% Nombre de points int�rieurs de la cha�ne
global Nn; 
Nn = N/2;

% Multiplicateurs de Lagrange (lm) d�finis arbitrairement
lm = [0.5077,0.4223]';

[e, c, g, a, hl, indic] = chs(1, xy, lm); %�tat initial

[x, lm, info] = sqp('chs', xy, lm, options);

[e, c, g, a, hl, indic] = chs(1, x, lm); %�tat final


%}

%%%%%%%%%%%%%%%%%%%%%%% CAS TEST 3b %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{

% Coordonn�es du dernier point
global A; A = 1;
global B; B = 0;

% Tailles des barres
global L; L = [2; 1];
global Nb; Nb = length(L);

% Tableau des coordonn�es des points int�rieurs et sa taille
xy = [0.5; 0.3];
global N; N = length(xy);

% Nombre de points int�rieurs de la cha�ne
global Nn; 
Nn = N/2;

% Multiplicateurs de Lagrange (lm) d�finis arbitrairement
lm = [0.5077,0.4223]';

[e, c, g, a, hl, indic] = chs(1, xy, lm); %�tat initial

[x, lm, info] = sqp('chs', xy, lm, options);

[e, c, g, a, hl, indic] = chs(1, x, lm); %�tat final

%}

%%%%%%%%%%%%%%%%%%%%%%% CAS TEST 3c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{

% Coordonn�es du dernier point
global A; A = 0;
global B; B = -1;

% Tailles des barres
global L; L = [2; 1];
global Nb; Nb = length(L);

% Tableau des coordonn�es des points int�rieurs et sa taille
xy = [0.3; 0.3];
global N; N = length(xy);

% Nombre de points int�rieurs de la cha�ne
global Nn; 
Nn = N/2;

% Multiplicateurs de Lagrange (lm) d�finis arbitrairement
lm = [0.5077,0.4223]';

[e, c, g, a, hl, indic] = chs(1, xy, lm); %�tat initial

[x, lm, info] = sqp('chs', xy, lm, options);

[e, c, g, a, hl, indic] = chs(1, x, lm); %�tat final

%}