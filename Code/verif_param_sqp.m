function status = verif_param_sqp(x, lm, tol1, tol2, maxit, verb, rl)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % verif_param :
  % Renvoie la valeur de l'indice de sortie :
  %   - 0 si les param�tres sont corrects
  %   - 1 sinon
  %
  % SYNOPSIS status = verif_param_sqp(x, lm, tol1, tol2, maxit)
  %          
  % INPUT * x : coordonn�es des points
  %         lm : multiplicateurs de Lagrange
  %         tol1, tol2 : tol�rances d'optimalit�
  %         maxit : nombre maximal d'it�rations autoris�es
  %         verb : type de tableau dans le cas de Newton avec RL
  %         rl : choix de recherche lin�aire ou pas
  %
  % OUTPUT - val: scalaire
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Variables globales
  global Nn; global Nb;
  
  %V�rification des param�tres :
  % - x est un vecteur-colonne de taille [2*Nn,1]
  % - lm est un vecteur-colonne de taille [Nb,1]
  % - tol1 et tol2 sont dans ]0,1[
  % - maxit est un entier strictement positif
  
  if !(verb == 1 || verb == 2) || !(rl == 0 || rl == 1) || !(tol1 > 0 && tol1 < 1) || !(tol2 > 0 && tol2 < 1) || !(maxit > 0 && mod(maxit,1)==0) || (!isvector(x) || !isequal(size(x),[2*Nn,1])) || (!isvector(lm) || !isequal(size(lm),[Nb,1]))
    % Les param�tres ne sont pas corrects
    status = 1;
  else
    % Les param�tres sont corrects
    status = 0;
  endif
  
endfunction
