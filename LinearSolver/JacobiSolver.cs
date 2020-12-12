using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace LinearSolver
{
    public static class JacobiSolver
    {
        public static double[] Solve(double[][] matrix, double[] b)
        {
            return JacobiMethod(matrix, b, 500);
        }

        private static double[] JacobiMethod(double[][] inputMatrix, double[] b, int iterations)
        {
            double[] solvedVector = new double[b.Count()];
            for (int p = 0; p < iterations; p++)
            {
                for (int i = 0; i < inputMatrix.Length; i++)
                {
                    double sigma = 0;
                    for (int j = 0; j < inputMatrix[0].Length; j++)
                    {
                        if (j != i)
                            sigma += inputMatrix[i][j] * solvedVector[j];
                    }
                    solvedVector[i] = (b[i] - sigma) / inputMatrix[i][i];
                }
            }
            return solvedVector;
        }
    }
}
