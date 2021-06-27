function [t] = diagonal_dominante(A);
  [n,n]= size (A);
   for (i=1:n)
     soma =0.0;
     for (j=1:n)
       soma += abs(A(i,j));
     endfor
     soma = soma-abs(A(i,i));
     if (soma >= abs(A(i,i)))
       t='FALSE';
       return;
     endif
   endfor
   t='TRUE';

endfunction
