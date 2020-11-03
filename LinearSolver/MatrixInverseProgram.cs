using System;

namespace MatrixInverse
{
    class MatrixInverseProgram
    {
        public static void Run(double[][] m, double[] b)
        {
            Console.WriteLine("\nBegin matrix inverse demo\n");

            Console.WriteLine("Original matrix m is ");
            MatShow(m, 4, 8);

            double d = MatDeterminant(m);
            if (Math.Abs(d) < 1.0e-5)
                Console.WriteLine("\nMatrix has no inverse");
            else
                Console.WriteLine("\nDet(m) = " + d.ToString("F4"));

            double[][] inv = MatInverse(m);
            Console.WriteLine("\nInverse matrix inv is ");
            MatShow(inv, 4, 8);

            Console.WriteLine("\nSolving m*x = b");
            double[] x = MatrixProduct(inv, b);  // (1, 0, 2, 1)

            Console.WriteLine("\nSolution x = ");
            VecShow(x, 1, 8);

            Console.WriteLine("\nEnd demo");
        } // Main

        static double[][] MatInverse(double[][] m)
        {
            // assumes determinant is not 0
            // that is, the matrix does have an inverse
            int n = m.Length;
            double[][] result = MatCreate(n, n); // make a copy
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    result[i][j] = m[i][j];

            double[][] lum; // combined lower & upper
            int[] perm;  // out parameter
            MatDecompose(m, out lum, out perm);  // ignore return

            double[] b = new double[n];
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                    if (i == perm[j])
                        b[j] = 1.0;
                    else
                        b[j] = 0.0;

                double[] x = Reduce(lum, b); // 
                for (int j = 0; j < n; ++j)
                    result[j][i] = x[j];
            }
            return result;
        }

        static int MatDecompose(double[][] m, out double[][] lum, out int[] perm)
        {
            // Crout's LU decomposition for matrix determinant and inverse
            // stores combined lower & upper in lum[][]
            // stores row permuations into perm[]
            // returns +1 or -1 according to even or odd number of row permutations
            // lower gets dummy 1.0s on diagonal (0.0s above)
            // upper gets lum values on diagonal (0.0s below)

            int toggle = +1; // even (+1) or odd (-1) row permutatuions
            int n = m.Length;

            // make a copy of m[][] into result lu[][]
            lum = MatCreate(n, n);
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < n; ++j)
                    lum[i][j] = m[i][j];

            // make perm[]
            perm = new int[n];
            for (int i = 0; i < n; ++i)
                perm[i] = i;

            for (int j = 0; j < n - 1; ++j) // process by column. note n-1 
            {
                double max = Math.Abs(lum[j][j]);
                int piv = j;

                for (int i = j + 1; i < n; ++i) // find pivot index
                {
                    double xij = Math.Abs(lum[i][j]);
                    if (xij > max)
                    {
                        max = xij;
                        piv = i;
                    }
                } // i

                if (piv != j)
                {
                    double[] tmp = lum[piv]; // swap rows j, piv
                    lum[piv] = lum[j];
                    lum[j] = tmp;

                    int t = perm[piv]; // swap perm elements
                    perm[piv] = perm[j];
                    perm[j] = t;

                    toggle = -toggle;
                }

                double xjj = lum[j][j];
                if (xjj != 0.0)
                {
                    for (int i = j + 1; i < n; ++i)
                    {
                        double xij = lum[i][j] / xjj;
                        lum[i][j] = xij;
                        for (int k = j + 1; k < n; ++k)
                            lum[i][k] -= xij * lum[j][k];
                    }
                }

            } // j

            return toggle;  // for determinant
        } // MatDecompose

        static double[] Reduce(double[][] luMatrix, double[] b) // helper
        {
            int n = luMatrix.Length;
            double[] x = new double[n];
            //b.CopyTo(x, 0);
            for (int i = 0; i < n; ++i)
                x[i] = b[i];

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
        } // Reduce

        static double MatDeterminant(double[][] m)
        {
            double[][] lum;
            int[] perm;
            double result = MatDecompose(m, out lum, out perm);  // impl. cast
            for (int i = 0; i < lum.Length; ++i)
                result *= lum[i][i];
            return result;
        }

        static double[][] MatCreate(int rows, int cols)
        {
            double[][] result = new double[rows][];
            for (int i = 0; i < rows; ++i)
                result[i] = new double[cols];
            return result;
        }

        static double[][] MatProduct(double[][] matA,
          double[][] matB)
        {
            int aRows = matA.Length;
            int aCols = matA[0].Length;
            int bRows = matB.Length;
            int bCols = matB[0].Length;
            if (aCols != bRows)
                throw new Exception("Non-conformable matrices");

            double[][] result = MatCreate(aRows, bCols);

            for (int i = 0; i < aRows; ++i) // each row of A
                for (int j = 0; j < bCols; ++j) // each col of B
                    for (int k = 0; k < aCols; ++k) // could use bRows
                        result[i][j] += matA[i][k] * matB[k][j];

            return result;
        }

        static void MatShow(double[][] m, int dec, int wid)
        {
            for (int i = 0; i < m.Length; ++i)
            {
                for (int j = 0; j < m[0].Length; ++j)
                {
                    double v = m[i][j];
                    if (Math.Abs(v) < 1.0e-5) v = 0.0;  // avoid "-0.00"
                    Console.Write(v.ToString("F" + dec).PadLeft(wid));
                }
                Console.WriteLine("");
            }
        }

        static void VecShow(int[] vec, int wid)
        {
            for (int i = 0; i < vec.Length; ++i)
                Console.Write(vec[i].ToString().PadLeft(wid));
            Console.WriteLine("");
        }

        static void VecShow(double[] vec, int dec, int wid)
        {
            for (int i = 0; i < vec.Length; ++i)
            {
                double x = vec[i];
                if (Math.Abs(x) < 1.0e-5) x = 0.0;  // avoid "-0.00"
                Console.Write(x.ToString("F" + dec).PadLeft(wid));
            }
            Console.WriteLine("");
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

    } // Program
} //ns
