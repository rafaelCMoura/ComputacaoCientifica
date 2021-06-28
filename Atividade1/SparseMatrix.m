classdef SparseMatrix
    properties (Access = public)
        fileName = ""
        compressedMatrix = [];
        numberOfNonZeros = 0;
        rateOfFill = 0;
        numberOfRows = 0;
        L = [];
        U = [];
        P = [];
        Q = [];
    endproperties

    methods (Access = public)
        function self = SparseMatrix(matrixName)
            % Load
            self.fileName = strcat(matrixName, ".mat");
            load(strcat("matrices/", self.fileName)); % Creates a variable called Problem
            self.compressedMatrix = Problem.A;
            
            % LUPQ decomposition (How about LUPQR?)
            [self.L, self.U, self.P, self.Q] = lu(self.compressedMatrix);

            % Number of nonzeros
            self.numberOfNonZeros = nnz(self.compressedMatrix);

            % Rate of Fill?
            self.rateOfFill = 100 - 100*self.numberOfNonZeros/(nnz(self.L) + nnz(self.U));
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