function [x_k, lm, info] = sqp (simul, x, lm, options)
  
  %%%%%%%%%%%%%%%%%%%%%%%% VERIFICATIONS PARAMETRES %%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  status = verif_param_sqp(x, lm, options.tol(1),options.tol(2),options.maxit, options.verb, options.rl);
    
  % Si les param�tres ne sont pas corrects, on arr�te la fonction
  if status == 1
    x_k=0; lm=0; info.niter=0; info.status = 1;
 
  elseif status == 0
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  
  x_k = x;
  n = length(x_k);
  m = length(lm);
  
  %%%%%%%%%%%%%%%%%%%%%%%%% CONDITIONS INITIALES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  [e_, c_, g_, a_, hl, indic_sortie] = feval(simul,5,x,lm);
  [e, c, g, a, hl_, indic_sortie] = feval(simul,4,x,lm);
  
  cond1 = max(abs(g + a.'*lm)); #condition 1 sur le gradient du lagrangien
  cond2 = max(abs(c)); #condition 2 sur l'activation des contraintes
  nb_iter = 0; #compteur du nombre d'it�rations effectu�es
  appels_sim = 0; #nombre d'appels au simulateur (cas options.verb=2)
  
  %il y a deux mani�res de stopper la boucle :
  % - seuil d'optimalit� atteint
  % - nombre maximal d'it�rations atteint
 
  % Affichage des donn�es
 if options.verb == 1 && options.rl == 0
  data = []; % Donn�s � afficher
  fprintf("----------------------------------------------------------------------------\n")
  fprintf("iter     |gl|         |ce|        |x|      |lm|       alpha        phi\n")
 endif
 
 while nb_iter < options.maxit && (cond1 > options.tol(1) || cond2 > options.tol(2))
   
   appels_sim_debut = appels_sim; %on sauvegarde le nombre d'appel au simulateur 
                                  %au d�but de l'it�ration k
   
   %%%%%%%%%%%%%%%%%%%%%%%%% SYST�ME NEWTON %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   #Calcul de l'it�r� x_k+1 :
   A_1 = [hl,a.'];
   A_2 = [a,zeros(m,m)];
   A = [A_1;A_2];    
   
   % La solution doit �tre calcul�e selon le cas options.rl == 1 ou == 0
      
   %%%%%%%%%%%%%%%%%%%%%%%%% NEWTON SANS RL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   if options.rl == 1
    
    b = [g;c];
    sol = A\b; 
     
    x_k = x_k - sol(1:n);
    
    #Calcul de l'it�r� lambda_k+1 :
    lm = - sol(n+1:n+m);
    
    %Mise � jour des conditions d'arr�t
    nb_iter += 1; %incr�mentation du nombre d'it�rations
    [e, c, g, a, hl_, indic_sortie] = feval(simul,4,x_k,lm);
    [e_, c_, g_, a_, hl, indic_sortie] = feval(simul,5,x_k,lm);
    cond1 = max(abs(g + a.'*lm)); #condition 1 sur le gradient du lagrangien
    cond2 = max(abs(c)); #condition 2 sur l'activation des contraintes
    
   %%%%%%%%%%%%%%%%%%%%%%%%% NEWTON AVEC RL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   elseif options.rl == 0     # On utilise la r�gle d'Armijo
   
    b = [g + a.' * lm; c];
    sol = A\b; 
    
    omega = 10**(-4); % Param�tre omega
    
    # On rappele que b = F(z) dans l'�nonc�
    phi_1 = (1/2) .* norm(b)**2; # Phi initial : phi(z_k)
    
    # Boucle pour tester les possibles pas
    
    % Pour options.verb = 2
    data_verb_2 = [];
    
    for i = 0:options.maxit   
     
     # Calcul du nouveau point
     alpha_k = (1/2) ** i;
     x_pas = x_k - alpha_k .* sol(1:n);
     lm_pas = lm - alpha_k .* sol(n+1:n+m);
     
     # On doit calculer phi dans le point d'arriv�.
     
     # Donn�es dans le nouveau point
     [e, c, g, a, hl_, indic_sortie] = feval(simul,4,x_pas,lm_pas);
     appels_sim += 1;
     
     b_pas = [g + a.' * lm_pas;c];   # F(z_k + alpha_k .* p_k)
     phi_pas = (1/2) .* norm(b_pas)**2;    # Phi dans le nouveau point
     
     % Stockage des donn�es pour options.verb = 2
     data_verb_2 = [data_verb_2; alpha_k (phi_pas - phi_1) ((phi_pas - phi_1) ./ alpha_k)];
     
     if phi_pas <= phi_1 .* (1 - 2 * omega * alpha_k) 
      break # Et la bonne valeur d'alpha est stock�e dans la variable alpha_k
     endif
  
    endfor
    
    % Une fois qu'on a le bon alpha, la nouvelle it�ration est x_k et lm_k :
    x_k = x_pas;
    lm = lm_pas;
    
    %Mise � jour des conditions d'arr�t
    [e_, c_, g_, a_, hl, indic_sortie] = feval(simul,5,x_k,lm);
    appels_sim += 1;
    cond1 = max(abs(g + a.'*lm)); #condition 1 sur le gradient du lagrangien
    cond2 = max(abs(c)); #condition 2 sur l'activation des contraintes    
    
    nb_iter += 1;
    
    % Stockage des donn�s pour affichage dans le cas options.verb = 1
    if options.verb == 1 && options.rl == 0
     data = [data; nb_iter max(abs(g + a.' * lm)) max(abs(c)) max(abs(x_pas)) max(abs(lm)) alpha_k phi_1];
    endif
    
    % Affichage du cas options.verb = 2
    if options.verb == 2 && options.rl == 0
    fprintf("---------------------------------------------------------------------------\n")
    formatSpec = "iter %i ,   simul %i ,       phi %.5e ,       pente %.5e \n \n";
    fprintf(formatSpec, [nb_iter appels_sim_debut phi_1 (-1) .* 2.* phi_1])
  
    formatSpec = "Recherche Lineaire d'Armijo : |d| = %.2e\n \n";
    fprintf(formatSpec, [max(abs(sol))])

    fprintf("  alpha         phip-phi        DF(phi) \n")
    formatSpec = "%.4e    %.5e    %.5e\n";
    
    for k = 1:size(data_verb_2)(1)
     fprintf(formatSpec, data_verb_2(k,:))
    endfor
  
    formatSpec = "\n|gl| = %.3e , |ce| = %.3e \n";
    fprintf(formatSpec, [max(abs(g + a.' * lm)) max(abs(c))])
    endif
    
   endif
   
 endwhile
 
 % Affichage des param�tres dans le cas options.verb = 1
 if options.verb == 1 && options.rl == 0
  formatSpec = '%2.i    %.4e   %.4e   %.1e   %.1e   %.3e   %.5e \n';
  for k = 1:nb_iter
   fprintf(formatSpec, data(k,:))
  endfor
  fprintf("---------------------------------------------------------------------------\n")
 endif
 
 %Stockage du nombre d'it�rations effectu�es (info.niter)
 info.niter = nb_iter; 
  
 %D�termination de la terminaison de l'algorithme (info.status = 0 ou 2)
 info.status = fun_terminaison(nb_iter,options.maxit);
 
 %Affichage des conditions d'optimalit�
 printf("Les conditions d'optimalites sont %e (gradient) et %e (contraintes).\n",[cond1,cond2])
 
 endif % If du status
  
 %Affichage du nombre d'it�rations de l'algorithme et du comportement de l'optimiseur
 fprintf("Le nombre maximal d'iterations autorisees est %d et le nombre d'iterations effectuees est %d \n",[options.maxit,info.niter]);
 fprintf("Le comportement de l'optimiseur est %d \n",[info.status]);
 
  [e_, c_, g_, a_, hl, indic_sortie] = feval(simul,5,x_k,lm);
  [e, c, g, a, hl_, indic_sortie] = feval(simul,4,x_k,lm);
  # V�rification de la hessienne r�duite
  # On trouve le noyau de la jacobienne des contraintes
  N = null(a);
  # On calcule la Hessienne r�duite:
  Hess_red = N.' * hl * N;
  # On calcule le spectre de la Hessienne r�duite:
  printf('Le spectre de la hessienne reduite est \n')
  
  full(eig(Hess_red))
  

endfunction