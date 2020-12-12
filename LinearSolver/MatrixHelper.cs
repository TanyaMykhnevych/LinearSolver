using System;

namespace LinearSolver
{
    public class MatrixHelper
    {
        public static double[][] MatrixCreate(int rows, int cols)
        {
            double[][] result = new double[rows][];
            for (int i = 0; i < rows; ++i)
                result[i] = new double[cols];
            return result;
        }

        public static double[][] MatrixDuplicate(double[][] matrix)
        {
            double[][] result = MatrixCreate(matrix.Length, matrix[0].Length);
            for (int i = 0; i < matrix.Length; ++i)
                for (int j = 0; j < matrix[i].Length; ++j)
                    result[i][j] = matrix[i][j];
            return result;
        }

        public static double[] MatrixProduct(double[][] matrixA, double[] matrixB)
        {
            int aRows = matrixA.Length;
            int aCols = matrixA[0].Length;
            int bRows = matrixB.Length;

            if (aCols != bRows)
                throw new Exception("Non-conformable matrices in MatrixProduct");

            double[] result = new double[aRows];

            for (int i = 0; i < aRows; ++i) // each row of A
                for (int k = 0; k < aCols; ++k) // could use k less-than bRows
                    result[i] += matrixA[i][k] * matrixB[k];

            return result;
        }

        public static double[][] SparseMatrix(double[][] matrix)
        {
            int blockSize = matrix.Length / 5;
            double[][] result = MatrixCreate(matrix.Length, matrix[0].Length);
            for (int i = 0; i < matrix.Length; i++)
            {
                int finish = i / blockSize * blockSize + blockSize;
                for (int j = i / blockSize * blockSize; j < finish; j++)
                {
                    result[i][j] = matrix[i][j];
                }
            }

            return result;
        }
    }
}
