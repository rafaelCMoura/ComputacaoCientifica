classdef LinearSystem
    properties
        A = [];
        b = [];
        x_gauss = [];
        x_jacobi = [];
        x_seidel = [];
        x_sor = [];
        w_sor = 1;
    endproperties

    methods
        function self = LinearSystem(matrixName, w_sor, tol, maxiter)
            sparseMatrix = SparseMatrix(matrixName);
            self.A = sparseMatrix.compressedMatrix;
            self.b = self.A*ones(sparseMatrix.numberOfRows, 1);

            % Gauss Elimination:
            self.x_gauss = self.A\self.b;
            
            % Jacobi
            [x_jacobi, er, iter] = self.jacobi(self.A, self.b, tol, maxiter, '.-b');
            self.x_jacobi = x_jacobi;

            % Gauss Seidel (SOR w=1)
            [x_seidel, er, iter] = self.sor(self.A, self.b, tol, maxiter, 1, '.-r');
            self.x_seidel = x_seidel;

            % SOR OverRelaxation (1<w<2) UnderRelaxation (0<w<1) Gauss Seidel (w=1)
            self.w_sor = w_sor;
            [x_sor, er, iter] = self.sor(self.A, self.b, tol, maxiter, self.w_sor, '.-k');
            self.x_sor = x_sor;
        endfunction

        function [x,er,iter] = jacobi(self,A,b,tol,nmaxiter, color)
            tic;
            [n,n]=size(A);
            iter = 1;
            er(1) = 1.0;  
            x0 = zeros(n,1);
            x = x0;
            while ((er(iter) > tol )&&(iter < nmaxiter))
                    printf("Jacobi: iteracao %d com erro = %f \n", iter, er(iter));
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
                    hold on;
                    set (gca, "yminorgrid", "on");
                    semilogy(iter, er(iter), color);
                    xlabel("iter - Número de Iteraçōes");
                    ylabel("log(er) - Norma do erro relativo");
                    title("Comparativo entre os métodos iterativos");
                    pause(0.1);
            endwhile;
            normx = norm(x,inf);
            printf("Jacobi - Convergencia obtida apos %d iteracoes\n",iter);
            printf("Norma do erro relativo =%f\n",er(iter));
            toc;
        endfunction

        function [x,er,iter] = sor(self,A,b,tol,nmaxiter,w, color)
	        tic;
            [n,n]=size(A);
            iter = 1;
            er(1) = 1.0;  
            x0 = zeros(n,1);
            x = x0;
            while (er(iter) > tol )&&(iter < nmaxiter)
                    printf("SOR w = %f: iteracao %d com erro = %f \n", w, iter, er(iter));
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
                    hold on;
                    set (gca, "yminorgrid", "on");
                    semilogy(iter, er(iter), color);
                    xlabel("iter - Número de Iteraçōes");
                    ylabel("log(er) - Norma do erro relativo");
                    title("Comparativo entre os métodos iterativos");
                    pause(0.1);
            endwhile;
            normx = norm(x,inf);
            printf("SOR - w =%f\n",w);
            printf("Convergencia obtida apos %d iteracoes\n",iter);
            printf("Norma do erro relativo =%f\n",er(iter));
            toc;
        endfunction;

    endmethods
endclassdef