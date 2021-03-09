function [e, c, g, a, hl, indic_sortie] = chs(indic, xy, lm)
  
  %%%%%%%%%%%%%%%%%%%%%%%% VARIABLES GLOBALES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  global A; global B; global Nb; global N; global Nn; global L;
  
  %%%%%%%%%%%%%%%%%%%%%%%% VERIFICATIONS PARAMETRES %%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  indic_sortie = verif_param(indic, xy, lm);
  
  % Si les param�tres ne sont pas corrects, on arr�te la fonction
  if indic_sortie == 1
    e=0;c=0;g=0;a=0;hl=0;
    return
  endif
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NOTATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Coordonn�es x et y de tous les points de la cha�ne
  % Il y a un d�calage +1 pour ces tableaux : x(j+1) = x_j et y(j+1) = y_j
  x = [0;xy(1:Nn);A]; 
  y = [0;xy(Nn + 1: 2 * Nn);B];
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  switch indic
    case 1  
      %repr�sentation de la cha�ne
      fun_trace(x,y);
      
      %variables par d�faut
      e=0;c=0;g=0;a=0;hl=0;
      
    case 2
      e = fun_e(x,y);
      c = fun_c(x,y);
      
      %variables par d�faut
      g=0;a=0;hl=0;
      
    case 4
      e = fun_e(x,y);
      c = fun_c(x,y);
      g = fun_g(x,y);
      a = fun_a(x,y);
      
      %variable par d�faut
      hl=0;
      
    case 5
      hl = fun_hl(x,y,lm);
      
      %variables par d�faut
      e=0;c=0;g=0;a=0;
      
    case 6
      test_derivees(xy,x,y);
      
      %variables par d�faut
      e=0;c=0;g=0;a=0;hl=0;
      
   endswitch

endfunction
