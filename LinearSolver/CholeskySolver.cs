using System;
using System.Linq;
using System.Threading.Tasks;

namespace LinearSolver
{
    //    double[][] matrix = new double[3][];
    //    matrix[0] = new double[3] { 81, -45, 45 };
    //    matrix[1] = new double[3] { -45, 50, -15 };
    //    matrix[2] = new double[3] { 45, -15, 38 };

    //    double[] b = new double[3] { 531, -460, 193 };

    // Lab 2
    public static class CholeskySolver
    {
        public static double[] Solve(double[][] matrix, double[] b, bool parallel = false)
        {
            var L = parallel ? CholeskyDecompositionParallel(matrix) : CholeskyDecomposition(matrix);

            var U = Transpose(L);

            var y = ForwardSubstitution(L, b);
            var x = BackwardSubstitution(U, y);

            return x;
        }

        private static double[][] CholeskyDecomposition(double[][] matrix)
        {
            double[][] lower = MatrixHelper.MatrixCreate(matrix.Length, matrix.Length);

            // Decomposing a matrix
            // into Lower Triangular
            for (int i = 0; i < matrix.Length; i++)
            {
                for (int j = 0; j <= i; j++)
                {
                    double sum = 0;

                    // Evaluating L(i, j)
                    // using L(j, j)
                    for (int k = 0; k < j; k++)
                        sum += (lower[i][k] * lower[j][k]);

                    // summation for diagnols
                    if (j == i)
                    {
                        lower[j][j] = Math.Sqrt(matrix[j][j] - sum);
                    }
                    else
                    {
                        lower[i][j] = (matrix[i][j] - sum) / lower[j][j];
                    }
                }
            }

            return lower;
        }

        private static double[][] CholeskyDecompositionParallel(double[][] matrix)
        {
            double[][] lower = MatrixHelper.MatrixCreate(matrix.Length, matrix.Length);

            // Decomposing a matrix
            // into Lower Triangular
            for (int i = 0; i < matrix.Length; i++)
            {
                for (int j = 0; j <= i; j++)
                {
                    // Evaluating L(i, j)
                    // using L(j, j)
                    double sum = 0;
                    Parallel.For(0, j, k =>
                    {
                        sum += (lower[i][k] * lower[j][k]);
                    });

                    // summation for diagnols
                    if (j == i)
                    {
                        lower[j][j] = Math.Sqrt(matrix[j][j] - sum);
                    }
                    else
                    {
                        lower[i][j] = (matrix[i][j] - sum) / lower[j][j];
                    }
                }
            }

            return lower;
        }

        private static double[][] Transpose(double[][] decomposed)
        {
            double[][] lower = MatrixHelper.MatrixCreate(decomposed.Length, decomposed.Length);

            for (int i = 0; i < decomposed.Length; i++)
            {
                for (int j = 0; j < decomposed.Length; j++)
                {
                    lower[j][i] = decomposed[i][j];
                }
            }

            return lower;
        }

        private static double[] ForwardSubstitution(double[][] matrix, double[] b)
        {
            double[] y = new double[b.Length];

            for (int r = 0; r < matrix.Length; ++r)
            {
                if (matrix[r][r] != 0)
                {
                    y[r] = (b[r] - Multiply(matrix[r], y)) / matrix[r][r];
                }
            }

            return y;
        }

        private static double[] BackwardSubstitution(double[][] matrix, double[] b)
        {
            double[] x = new double[b.Length];

            for (int r = matrix.Length - 1; r >= 0; --r)
            {
                if (matrix[r][r] != 0)
                {
                    x[r] = (b[r] - Multiply(matrix[r], x)) / matrix[r][r];
                }
            }
            return x;
        }

        private static double Multiply(double[] a, double[] b)
        {
            double productSum = 0;
            for (int i = 0; i < a.Length; i++)
            {
                productSum += a[i] * b[i];
            }

            return productSum;
        }
    }
}
