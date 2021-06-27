function [x,er,iter]=sor(A,b,tol,nmaxiter,w)
	   tic;
     [n,n]=size(A);
     iter = 1;
     er(1) = 1.0;  
     x0 = zeros(n,1);
     x = x0;
     while (er(iter) > tol )&&(iter < nmaxiter)

            for i=1:n
                soma = 0.0;
                for j = 1:(i-1)
                   soma = soma + A(i,j)*x(j);
                endfor
                for j = (i+1):n
		                soma = soma + A(i,j)*x0(j);
                endfor
	              x(i) = w*(b(i) - soma)/A(i,i) + (1-w)*x0(i);
            endfor
            iter = iter + 1;
	          er(iter) = norm(x-x0,inf)/norm(x,inf);            
            x0 = x;
            x	    
     endwhile;
     normx = norm(x,inf);
     printf("SOR - w=%f\n",w);
     printf("Convergencia obtida apos %d iteracoes\n",iter);
     printf("Norma do erro relativo =%f\n",er(iter));
	   toc;
endfunction;
            
     

           
