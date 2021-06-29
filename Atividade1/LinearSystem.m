classdef LinearSystem
    properties
        A = [];
        b = [];
    endproperties

    methods
        function self = LinearSystem(matrixName)
            sparseMatrix = SparseMatrix(matrixName);
            self.A = sparseMatrix.A;
            self.b = self.A*ones(sparseMatrix.numberOfRows);
        endfunction

        function solve(self)
            %TODO make a switch from methods, gauss, jacobi, sor nargin, gausssedel, all
            % solve method must call show method with log chart
        endfunction

        function gaussianElimination(self)
            %x = A\b
        endfunction

        function [x,er,iter] = jacobi(A,b,tol,nmaxiter)
            tic;
            [n,n]=size(A);
            iter = 1;
            er(1) = 1.0;  
            x0 = zeros(n,1);
            x = x0;
            while ((er(iter) > tol )&&(iter < nmaxiter))

                    for i=1:n
                        soma = 0.0;
                        for j = 1:(i-1)
                        soma = soma + A(i,j)*x0(j);
                        endfor
                        for j = (i+1):n
                                soma = soma + A(i,j)*x0(j);
                        endfor
                        x(i) = (b(i) - soma)/A(i,i);
                    endfor
                    iter = iter + 1;
                    er(iter) = norm(x-x0,inf)/norm(x,inf);            
                    x0 = x;	    
            endwhile;
            normx = norm(x,inf);
            printf("Jacobi - Convergencia obtida apos %d iteracoes\n",iter);
            printf("Norma do erro relativo =%f\n",er(iter));
            toc;
        endfunction

        function [x,er,iter] = sor(A,b,tol,nmaxiter,w)
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

        function show(self)
            
        endfunction
    endmethods
endclassdef