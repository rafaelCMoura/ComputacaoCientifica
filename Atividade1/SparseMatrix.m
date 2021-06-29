classdef SparseMatrix
    properties
        fileName = ""
        compressedMatrix = [];
        L = [];
        U = [];
        P = [];
        Q = [];
        numberOfNonZeros = 0;
        rateOfFill = 0;
        numberOfRows = 0;
        isDiagonallyDominant = -1;
        conditionNumber = 0;
    endproperties

    methods
        function self = SparseMatrix(matrixName)
            % Load
            self.fileName = strcat(matrixName, ".mat");
            load(strcat("matrices/", self.fileName)); % Creates a variable called Problem
            self.compressedMatrix = Problem.A;
            
            % LUPQ decomposition (How about LUPQR?)
            [self.L, self.U, self.P, self.Q] = lu(self.compressedMatrix);

            % Diagonally Dominant
            self.isDiagonallyDominant = self.diagonallyDominant();

            % Number of nonzeros
            self.numberOfNonZeros = nnz(self.compressedMatrix);

            % Number of rows
            self.numberOfRows = rows(self.compressedMatrix);

            % Rate of Fill?
            self.rateOfFill = 100 - 100*self.numberOfNonZeros/(nnz(self.L) + nnz(self.U));

            % Condition Number
            self.conditionNumber = cond(self.compressedMatrix);
        endfunction

        function isDiagonallyDominant = diagonallyDominant(self);
            n = self.numberOfRows;
            for (i=1:n)
                summ = 0.0;
                for (j=1:n)
                    summ += abs(self.compressedMatrix(i,j));
                endfor
                summ = summ - abs(self.compressedMatrix(i,i));
                if (summ >= abs(self.compressedMatrix(i,i)))
                    isDiagonallyDominant = 0;
                    return;
                endif
            endfor
            isDiagonallyDominant = 1;
        endfunction

        function show(self)
            subplot(2,2,1);
            spy(self.compressedMatrix, '.b');
            title(self.fileName);

            subplot(2,2,2);
            spy(self.L, '.b');
            title("L");

            subplot(2,2,3);
            spy(self.U, '.b');
            title("U");

            subplot(2,2,4);
            spy(self.P, '.b');
            title("P");
        endfunction
    endmethods
endclassdef