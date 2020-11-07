using System;
using System.Diagnostics;
using System.Drawing;
using System.IO;

namespace LinearSolver
{
    class Program
    {
        static void Main(string[] args)
        {
            Bitmap bpm = new Bitmap(Image.FromFile(Path.Combine(@"..\..\..\Images\picture.jpg")));

            Solve(bpm, 10);
            Solve(bpm, 100);
            Solve(bpm, 500);

            Console.ReadLine();
        }

        private static void Solve(Bitmap bpm, int n)
        {
            Console.WriteLine($"Started solving for N = {n}");

            double[][] matrix = GetGreenValuesMatrix(bpm, n);
            double[] b = GetGreenValuesColumn(bpm, n);

            Stopwatch watch = Stopwatch.StartNew();
            double[] resSeq = SequentialMatrixSolver.Solve(matrix, b);
            watch.Stop();
            Console.WriteLine($"Sequential solving time = {watch.ElapsedMilliseconds} ms");

            watch.Reset();
            watch.Start();
            double[] resPar = ParallelMatrixSolver.Solve(matrix, b);
            watch.Stop();
            Console.WriteLine($"Parallel solving time = {watch.ElapsedMilliseconds} ms");

            Console.WriteLine();
        }

        private static double[][] GetGreenValuesMatrix(Bitmap bmp, int n)
        {
            double[][] result = MatrixHelper.MatrixCreate(n, n);

            for (int x = 0; x < n; x++)
            {
                result[x] = new double[n];
                for (int y = 0; y < n; y++) result[x][y] = bmp.GetPixel(x, y).G;
            }

            return result;
        }

        private static double[] GetGreenValuesColumn(Bitmap bmp, int n)
        {
            double[] result = new double[n];

            for (int y = 0; y < n; y++) result[y] = bmp.GetPixel(bmp.Width - 1, y).G;

            return result;
        }
    }
}
