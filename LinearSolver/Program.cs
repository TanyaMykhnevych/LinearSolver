using MatrixInverse;
using System;
using System.Drawing;
using System.IO;

namespace LinearSolver
{
    class Program
    {
        static void Main(string[] args)
        {
            Bitmap bpm = new Bitmap(Image.FromFile(Path.Combine(@"..\..\..\Images\picture.jpg")));

            var matrix = GetGreenValuesMatrix(bpm, 3);
            var b = GetGreenValuesColumn(bpm, 3);

            MatrixInverseProgram.Run(matrix, b);

            var inversedMatrix = MatrixInverse(matrix);

            var res = MatrixProduct(inversedMatrix, b);

            Console.WriteLine("Hello World!");
        }

        private static double[][] GetGreenValuesMatrix(Bitmap bmp, int n)
        {
            double[][] result = MatrixHelper.MatrixCreate(n, n);

            for (int x = 0; x < n; x++)
            {
                result[x] = new double[n];
                for (int y = 0; y < n; y++)
                    result[x][y] = bmp.GetPixel(x, y).G;
            }

            return result;
        }

        private static double[] GetGreenValuesColumn(Bitmap bmp, int n)
        {
            double[] result = new double[n];

            for (int y = 0; y < n; y++)
                result[y] = bmp.GetPixel(bmp.Width - 1, y).G;

            return result;
        }

        static double[][] MatrixInverse(double[][] matrix)
        {
            int n = matrix.Length;
            double[][] result = MatrixHelper.MatrixDuplicate(matrix);

            int[] perm;
            int toggle;
            double[][] lum = MatrixDecompose(matrix, out perm, out toggle);
            if (lum == null)
                throw new Exception("Unable to compute inverse");

            double[] b = new double[n];
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    if (i == perm[j])
                        b[j] = 1.0;
                    else
                        b[j] = 0.0;
                }

                double[] x = HelperSolve(lum, b);

                for (int j = 0; j < n; ++j)
                    result[j][i] = x[j];
            }
            return result;
        }

        static double[] HelperSolve(double[][] luMatrix, double[] b)
        {
            // before calling this helper, permute b using the perm array
            // from MatrixDecompose that generated luMatrix
            int n = luMatrix.Length;
            double[] x = new double[n];
            b.CopyTo(x, 0);

            for (int i = 1; i < n; ++i)
            {
                double sum = x[i];
                for (int j = 0; j < i; ++j)
                    sum -= luMatrix[i][j] * x[j];
                x[i] = sum;
            }

            x[n - 1] /= luMatrix[n - 1][n - 1];
            for (int i = n - 2; i >= 0; --i)
            {
                double sum = x[i];
                for (int j = i + 1; j < n; ++j)
                    sum -= luMatrix[i][j] * x[j];
                x[i] = sum / luMatrix[i][i];
            }

            return x;
        }

        // https://stackoverflow.com/questions/46836908/double-inversion-c-sharp
        static double[][] MatrixDecompose(double[][] matrix, out int[] perm, out int toggle)
        {
            // Doolittle LUP decomposition with partial pivoting.
            // rerturns: result is L (with 1s on diagonal) and U;
            // perm holds row permutations; toggle is +1 or -1 (even or odd)
            int rows = matrix.Length;
            int cols = matrix[0].Length; // assume square
            if (rows != cols)
                throw new Exception("Attempt to decompose a non-square m");

            int n = rows; // convenience

            double[][] result = MatrixHelper.MatrixDuplicate(matrix);

            perm = new int[n]; // set up row permutation result
            for (int i = 0; i < n; ++i) { perm[i] = i; }

            toggle = 1; // toggle tracks row swaps.
                        // +1 -greater-than even, -1 -greater-than odd. used by MatrixDeterminant

            for (int j = 0; j < n - 1; ++j) // each column
            {
                double colMax = Math.Abs(result[j][j]); // find largest val in col
                int pRow = j;

                // reader Matt V needed this:
                for (int i = j + 1; i < n; ++i)
                {
                    if (Math.Abs(result[i][j]) > colMax)
                    {
                        colMax = Math.Abs(result[i][j]);
                        pRow = i;
                    }
                }
                // Not sure if this approach is needed always, or not.

                if (pRow != j) // if largest value not on pivot, swap rows
                {
                    double[] rowPtr = result[pRow];
                    result[pRow] = result[j];
                    result[j] = rowPtr;

                    int tmp = perm[pRow]; // and swap perm info
                    perm[pRow] = perm[j];
                    perm[j] = tmp;

                    toggle = -toggle; // adjust the row-swap toggle
                }
                if (result[j][j] == 0.0)
                {
                    // find a good row to swap
                    int goodRow = -1;
                    for (int row = j + 1; row < n; ++row)
                    {
                        if (result[row][j] != 0.0)
                            goodRow = row;
                    }

                    if (goodRow == -1)
                        throw new Exception("Cannot use Doolittle's method");

                    // swap rows so 0.0 no longer on diagonal
                    double[] rowPtr = result[goodRow];
                    result[goodRow] = result[j];
                    result[j] = rowPtr;

                    int tmp = perm[goodRow]; // and swap perm info
                    perm[goodRow] = perm[j];
                    perm[j] = tmp;

                    toggle = -toggle; // adjust the row-swap toggle
                }
                for (int i = j + 1; i < n; ++i)
                {
                    result[i][j] /= result[j][j];
                    for (int k = j + 1; k < n; ++k)
                    {
                        result[i][k] -= result[i][j] * result[j][k];
                    }
                }

            }

            return result;
        }

        static double[] MatrixProduct(double[][] matrixA, double[] matrixB)
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
    }
}
