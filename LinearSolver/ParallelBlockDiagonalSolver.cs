using System;
using System.Threading.Tasks;

namespace LinearSolver
{
    public static class ParallelBlockDiagonalSolver
    {
        public const int BlockSize = 25;

        public static double[] Solve(double[][] matrix, double[] b, bool parallel = false)
        {
            double[][] bdMatrix = ConvertToBlockDiagonal(matrix, BlockSize);

            int h = BlockSize / 2;
            int m = BlockSize / 2;


            ParallelTriangulation(bdMatrix, b, matrix.Length, BlockSize, m, h);

            for (int i = 0; i < m; i++)
            {
                for (int j = m + BlockSize - 1; j > i + 1 + h; j--)
                {
                    double multiplier = -bdMatrix[j][i] / bdMatrix[i][i];
                    int columnBound = Math.Min(matrix.Length, i + BlockSize);

                    for (int k = i; k < columnBound; k++)
                    {
                        bdMatrix[j][k] += bdMatrix[i][k] * multiplier;
                    }

                    b[j] += b[i] * multiplier;
                }
            }

            double[] answers = SolveInParallel(bdMatrix, b, matrix.Length, BlockSize, m);

            return answers;
        }

        private static double[][] ConvertToBlockDiagonal(double[][] matrix, double k)
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

        private static void ParallelTriangulation(double[][] dMatrix, double[] b, int size, int k, int m, int h)
        {
            Task task1 = Task.Run(() =>
            {
                for (int i = 0; i < m; i++)
                {
                    for (int j = i + 1; j < i + 1 + h; j++)
                    {
                        double multiplier = -dMatrix[j][i] / dMatrix[i][i];
                        int columnBound = Math.Min(size, i + k);

                        for (int k = i; k < columnBound; k++)
                        {
                            dMatrix[j][k] += dMatrix[i][k] * multiplier;
                        }

                        b[j] += b[i] * multiplier;
                    }
                }
            });

            Task task2 = Task.Run(() =>
            {
                for (int i = size - 1; i > m + k; i--)
                {
                    for (int j = i + 1; j < i + 1 - h; j++)
                    {
                        double multiplier = -dMatrix[j][i] / dMatrix[i][i];
                        int columnBound = Math.Max(0, i - k);

                        for (int k = i; k >= columnBound; k--)
                        {
                            dMatrix[j][k] += dMatrix[i][k] * multiplier;
                        }

                        b[j] += b[i] * multiplier;
                    }
                }
            });

            Task.WaitAll(task1, task2);
        }

        private static double[] SolveInParallel(double[][] dMatrix, double[] dResults, int size, int kSize, int m)
        {
            double[] answers = new double[size];

            Task task1 = Task.Run(() =>
            {
                for (int i = m - 1; i >= 0; i--)
                {
                    double currentResult = dResults[i];
                    int columnBound = Math.Min(size, i + 1 + kSize / 2);

                    for (int j = i + 1; j < columnBound; j++)
                    {
                        currentResult -= answers[j] * dMatrix[i][j];
                    }

                    answers[i] = currentResult / dMatrix[i][i];
                }
            });

            Task task2 = Task.Run(() =>
            {
                for (int i = m; i < size; i++)
                {
                    double currentResult = dResults[i];
                    int columnBound = Math.Min(size, i + 1 + kSize / 2);

                    for (int j = i + 1; j < columnBound; j++)
                    {
                        currentResult -= answers[j] * dMatrix[i][j];
                    }

                    answers[i] = currentResult / dMatrix[i][i];
                }
            });

            Task.WaitAll(task1, task2);

            return answers;
        }
    }
}
