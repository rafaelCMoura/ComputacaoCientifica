classdef LinearSystem
    properties
        matrixName = "";
        A = [];
        b = [];
        x_gauss = [];
        x_jacobi = [];
        x_seidel = [];
        x_sor = [];
        w_sor = 1;
    endproperties

    methods
        function self = LinearSystem(matrixName, w_sor, tol, maxiter, follow)
            self.matrixName = matrixName;

            sparseMatrix = SparseMatrix(matrixName);
            self.A = sparseMatrix.compressedMatrix;
            self.b = self.A*ones(sparseMatrix.numberOfRows, 1);

            % Gauss Elimination:
            self.x_gauss = self.A\self.b;
            
            figure(1);
            % Jacobi
            [x_jacobi, er, iter] = self.jacobi(self.A, self.b, tol, maxiter, '-b', follow);
            self.x_jacobi = x_jacobi;

            % Gauss Seidel (SOR w=1)
            [x_seidel, er, iter] = self.sor(self.A, self.b, tol, maxiter, 1, '-r', follow, 1);
            self.x_seidel = x_seidel;

            % SOR OverRelaxation (1<w<2) UnderRelaxation (0<w<1) Gauss Seidel (w=1)
            self.w_sor = w_sor;
            [x_sor, er, iter] = self.sor(self.A, self.b, tol, maxiter, self.w_sor, '-k', follow, 1);
            self.x_sor = x_sor;

            figure(2);
            w = 0.2:0.2:1.8;
            erro = zeros(length(w),1);
            for i=1:length(w)
                [x_sor, er, iter] = self.sor(self.A, self.b, tol, maxiter, w(i), '-k', 0, 0);
                erro(i) = er(iter);
            end
            grid on;
            plot(w, erro, ".-b", "markersize", 10);
            xlabel("fator w do SOR");
            ylabel("Erro na ultima iteracao");
            title("Comparativo entre o SOR para diferentes valores de w");

        endfunction

        function [x,er,iter] = jacobi(self,A,b,tol,nmaxiter, color, follow)
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
                    if(follow==1)
                        printf("Jacobi: iteracao %d com erro = %f \n", iter, er(iter));
                        hold on;
                        set (gca, "yminorgrid", "on");
                        semilogy(iter, er(iter), color, "markersize", 10);
                        xlabel("iter - Numero de Iteracoes");
                        ylabel("log(er) - Norma do erro relativo");
                        title(strcat("Comparativo entre os metodos iterativos: Matriz -", self.matrixName));
                        pause(0.05);
                    end
            endwhile;
            normx = norm(x,inf);
            printf("Jacobi - Convergencia obtida apos %d iteracoes\n",iter);
            printf("Norma do erro relativo =%f\n",er(iter));
            toc;

            hold on;
            set (gca, "yminorgrid", "on");
            semilogy(1:iter, er, color, "markersize", 10);
            xlabel("iter - Numero de Iteracoes");
            ylabel("log(er) - Norma do erro relativo");
            title(strcat("Comparativo entre os metodos iterativos: Matriz -", self.matrixName));
            pause(0.05);
        endfunction

        function [x,er,iter] = sor(self,A,b,tol,nmaxiter,w, color, follow, final_plot)
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
                    if(follow==1)
                        printf("SOR w = %f: iteracao %d com erro = %f \n", w, iter, er(iter));
                        hold on;
                        set (gca, "yminorgrid", "on");
                        semilogy(iter, er(iter), color, "markersize", 10);
                        xlabel("iter - Numero de Iteracoes");
                        ylabel("log(er) - Norma do erro relativo");
                        title(strcat("Comparativo entre os metodos iterativos: Matriz -", self.matrixName));
                        pause(0.05);
                    end
            endwhile;
            normx = norm(x,inf);
            printf("SOR - w =%f\n",w);
            printf("Convergencia obtida apos %d iteracoes\n",iter);
            printf("Norma do erro relativo =%f\n",er(iter));
            toc;

            if (final_plot==1)
                hold on;
                set (gca, "yminorgrid", "on");
                semilogy(1:iter, er, color, "markersize", 10);
                xlabel("iter - Numero de Iteracoes");
                ylabel("log(er) - Norma do erro relativo");
                title(strcat("Comparativo entre os metodos iterativos: Matriz -", self.matrixName));
                pause(0.05);
            end
        endfunction;

    endmethods
endclassdef