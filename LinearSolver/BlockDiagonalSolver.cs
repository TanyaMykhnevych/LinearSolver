using System;

namespace LinearSolver
{
    public static class BlockDiagonalSolver
    {
        public const int BlockSize = 25;

        public static double[] Solve(double[][] matrix, double[] b, bool parallel = false)
        {
            double[][] bdMatrix = ConvertToBlockDiagonal(matrix);

            Triangulate(bdMatrix, b, matrix.Length, BlockSize);
            double[] answers = FindAnswers(bdMatrix, b, BlockSize);

            return answers;
        }

        private static double[][] ConvertToBlockDiagonal(double[][] matrix)
        {
            double[][] result = MatrixHelper.MatrixCreate(matrix.Length, matrix[0].Length);
            for (int i = 0; i < matrix.Length; i++)
            {
                int finish = i / BlockSize * BlockSize + BlockSize;
                for (int j = i / BlockSize * BlockSize; j < finish; j++)
                {
                    result[i][j] = matrix[i][j];
                }
            }

            return result;
        }

        private static void Triangulate(double[][] dMatrix, double[] dResults, int size, int kSize)
        {
            for (int i = 0; i < size - 1; i++)
            {
                int rowBound = Math.Min(size, i + 1 + kSize / 2);

                for (int j = i + 1; j < rowBound; j++)
                {
                    double multiplier = -dMatrix[j][i] / dMatrix[i][i];
                    int columnBound = Math.Min(size, i + kSize);

                    for (int k = i; k < columnBound; k++)
                    {
                        dMatrix[j][k] += dMatrix[i][k] * multiplier;
                    }

                    dResults[j] += dResults[i] * multiplier;
                }
            }
        }

        private static double[] FindAnswers(double[][] triangulated, double[] results, int kSize)
        {
            int size = triangulated.GetLength(0);
            double[] answers = new double[size];

            for (int i = size - 1; i >= 0; i--)
            {
                double currentResult = results[i];
                int columnBound = Math.Min(size, i + 1 + kSize / 2);

                for (int j = i + 1; j < columnBound; j++)
                {
                    currentResult -= answers[j] * triangulated[i][j];
                }

                answers[i] = currentResult / triangulated[i][i];
            }

            return answers;
        }
    }
}
