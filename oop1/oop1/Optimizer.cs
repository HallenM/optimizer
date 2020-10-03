using System;
using System.Collections.Generic;

namespace optimizer
{
    class CurrentProgress : IProgress<(double[] current, double residual, int progresslen, int progressval)>
    {
        public void Report((double[] current, double residual, int progresslen, int progressval) value)
        {
            if (value.progressval % 1000 == 0 || value.progressval <= value.progresslen)
            {
                Console.Write("Residual: \t" + value.residual.ToString() + "\n");
                Console.Write("Curr time: " + value.progressval.ToString() + "\tMax time: " + value.progresslen.ToString() + "\n");
                Console.Write("Progress: " + (Convert.ToInt32(Convert.ToDouble(value.progressval) / Convert.ToDouble(value.progresslen) * 100.0)).ToString() + "\n");
                Console.Write("\nValue of parameters:\n");
                for (int i = 0; i < value.current.Length; i++)
                    Console.Write(value.current[i].ToString() + "\n");
                Console.Write("=============================================================\n");
                Console.Write("\nCalculation...\n\n");
            }
        }
    }

    class SimpleRandomSearch : interfaces.Optimizer
    {
        private double eps;
        public interfaces.FunctionalWithDiff F;
        private DateTime time;

        // a, b числа определяющие интервал
        double[] Randomize(double[] a, double[] b, int len)
        {
            double[] p = new double[len];
            Random rand = new Random();
            for (int i = 0; i < len; i++)
            {
                // Рандомное равномено распределённое число
                double x = rand.NextDouble();
                p[i] = a[i] + (x * (b[i] - a[i]));
            }
            return p;
        }

        interfaces.FunctionalWithDiff interfaces.Optimizer.functional
        {
            set
            {
                F = value;
            }
        }

        double[] interfaces.Optimizer.Optimize(double[] initial,
            IProgress<(double[] current, double residual, int progresslen, int progressval)> progress)
        {
            int userTime = time.Hour * 60 * 60 * 1000 + time.Minute * 60 * 1000 + time.Second * 1000 + time.Millisecond;
            int timesGone = 0;
            (double[] current, double residual, int progresslen, int progressval) realProgress = (new double[initial.Length], 0.0, 0, 0);
            realProgress.progresslen = userTime;
            System.Diagnostics.Stopwatch clock = new System.Diagnostics.Stopwatch();

            // Случайный поиск
            // Выбранная вероятность
            double P = 0.2;
            uint maxiter;
            uint m = 60000;
            double f1 = 0, f2 = 0;
            bool firstLaunch = true;
            bool isTime = true;

            // Диапазоны изменения параметров
            double[] minValue = new double[initial.Length];
            double[] maxValue = new double[initial.Length];
            double[] x = new double[initial.Length];
            double[] xNext = new double[initial.Length];

            for (int i = 0; i < initial.Length; i++)
            {
                minValue[i] = F.Range[i].min;
                maxValue[i] = F.Range[i].max;
                // Присвавивание начальной равномерно распределённой точки
                x[i] = xNext[i] = initial[i];
            }

            Console.Write("\nCalculation...\n\n");
            clock.Start();
            while (isTime && timesGone < userTime)
            {
                if (firstLaunch)
                {
                    double V = 1, Veps = 1;
                    for (int i = 0; i < minValue.Length; i++)
                    {
                        V *= maxValue[i] - minValue[i];
                        double epsParam = eps * (maxValue[i] - minValue[i]);
                        Veps *= /*eps; //*/epsParam;
                    }
                    double Peps = Veps / V;
                    maxiter = (uint)(Math.Log(1 - P) / Math.Log(1 - Peps)) / 10000;
                }
                else maxiter = m;

                for (int i = 1; i <= maxiter; i++)
                {
                    xNext = Randomize(minValue, maxValue, x.Length);
                    f1 = F.Value(x);
                    f2 = F.Value(xNext);
                    if (f1 > f2)
                    {
                        for (int j = 0; j < minValue.Length; j++)
                            x[j] = xNext[j];

                        realProgress.residual = f2;
                        timesGone = Convert.ToInt32(clock.ElapsedMilliseconds);
                        x.CopyTo(realProgress.current, 0);
                        realProgress.progressval = timesGone;
                        progress.Report(realProgress);
                    }
                }
                firstLaunch = !firstLaunch;
                timesGone = Convert.ToInt32(clock.ElapsedMilliseconds);
            }
            clock.Stop();
            Console.Write("Calculation Finished!\n");
            return x;
        }

        double interfaces.Optimizer.Eps
        {
            set 
            {
                eps = value;
            }
        }

        DateTime interfaces.Optimizer.MaxTime
        {
            set
            {
                time = value;
            }
        }
    }

    // Покоординатный спуск
    class Gauss : interfaces.Optimizer
    {
        private double eps;
        private interfaces.FunctionalWithDiff F;
        private DateTime time;

        double RandValue(double a, double b)
        {
            Random rand = new Random();
            // Рандомное равномено распределённое число
            double x = rand.NextDouble();
            return a + (x * (b - a));
        }

        interfaces.FunctionalWithDiff interfaces.Optimizer.functional
        {
            set
            {
                F = value;
            }
        }

        double Dichotomy(double[] xk, double ai, double bi, int i)
        {
            double[] x1 = new double[xk.Length];
            double[] x2 = new double[xk.Length];
            double delta = 1e-7;

            for (int j = 0; j < xk.Length; j++) 
                x1[j] = x2[j] = xk[j];

            double f1 = 0, f2 = 0;

            while (Math.Abs(bi - ai) >= eps)
            {
                x1[i] = (ai + bi - delta) / 2;
                x2[i] = (ai + bi + delta) / 2;

                f1 = F.Value(x1);
                f2 = F.Value(x2);

                if (f1 < f2)
                {
                    // a1 = a0;
                    bi = x2[i];
                }
                else
                {
                    ai = x1[i];
                    // b1 = b0;
                }
            }
            if (f1 < f2) return x1[i];
            else return x2[i];
        }

        double[] interfaces.Optimizer.Optimize(double[] initial,
            IProgress<(double[] current, double residual, int progresslen, int progressval)> progress)
        {
            double[] result = new double[initial.Length];

            int userTime = time.Hour * 60 * 60 * 1000 + time.Minute * 60 * 1000 + time.Second * 1000 + time.Millisecond;
            int timesGone = 0;
            (double[] current, double residual, int progresslen, int progressval) realProgress = (new double[initial.Length], 0.0, 0, 0);
            realProgress.progresslen = userTime;
            System.Diagnostics.Stopwatch clock = new System.Diagnostics.Stopwatch();

            // Покоординатный спуск
            int i = 0;
            uint maxiter = 10000;
            bool isTime = true;
            // Диапазоны изменения параметров
            double[] minValue = new double[initial.Length];
            double[] maxValue = new double[initial.Length];
            double[] xk = new double[initial.Length];
            double[] xkNext = new double[initial.Length];

            for (int j = 0; j < minValue.Length; j++)
            {
                minValue[j] = F.Range[j].min;
                maxValue[j] = F.Range[j].max;
                xk[j] = initial[j];
                xkNext[j] = initial[j] + 10;
            }

            double f1 = F.Value(xk), f2 = F.Value(xkNext);
            int iter = 0;
            Console.Write("\nCalculation...\n\n");

            clock.Start();
            while (isTime && timesGone < userTime)
            {
                //isTime = false;
                int curTime = Convert.ToInt32(clock.ElapsedMilliseconds);
                if (curTime % 1000 == 0)
                    Console.Write("Time: " + curTime + "\n");
                while (Math.Sqrt((f2 - f1) * (f2 - f1)) > eps || iter < maxiter)
                {
                    if (i == xk.Length) { i = 0; iter++; }
                    if (f2 == f1) break;
                    for (; i < xk.Length; i++)
                    {
                        xkNext[i] = Dichotomy(xk, minValue[i], maxValue[i], i);
                        
                        f1 = F.Value(xk);
                        f2 = F.Value(xkNext);
                        if (f1 > f2)
                        {
                            for (int j = 0; j < minValue.Length; j++)
                                xk[j] = xkNext[j];

                            /*Console.WriteLine(xk[0] + "\t" + xk[1] + "\n");
                            Console.WriteLine(xkNext[0] + "\t" + xkNext[1] + "\t" + iter + "\n");
                            Console.WriteLine("Value F " + f1 + "\t" + f2 + "\n");*/

                            timesGone = Convert.ToInt32(clock.ElapsedMilliseconds);
                            realProgress.residual = f2;
                            xk.CopyTo(realProgress.current, 0);
                            realProgress.progressval = timesGone;
                            progress.Report(realProgress);
                        }
                    }
                }
                timesGone = Convert.ToInt32(clock.ElapsedMilliseconds);
            }
            clock.Stop();
            Console.Write("Calculation Finished!\n");
            f1 = F.Value(xk); ; f2 = F.Value(xkNext);
            if (f1 > f2) return xkNext;
            else return xk;
        }

        double interfaces.Optimizer.Eps
        {
            set
            {
                eps = value;
            }
        }

        DateTime interfaces.Optimizer.MaxTime
        {
            set
            {
                time = value;
            }
        }

    }

    class MSG : interfaces.Optimizer
    {
        private double eps;
        private interfaces.FunctionalWithDiff F;
        private DateTime time;

        double Norm(double[] x)
        {
            double res = 0;
            for (int i = 0; i < x.Length; i++)
            {
                res += x[i] * x[i];
            }
            return Math.Sqrt(res);
        }

        double CalcOmega(double[] xk, double[] xkPrev)
        {
            double[] gradFk = new double[xk.Length], gradFk_1 = new double[xkPrev.Length];
            for (int i = 0; i < gradFk.Length; i++)
            {
                gradFk[i] = F.DfDp(i, xk);
                gradFk_1[i] = F.DfDp(i, xkPrev);
            }
            double omega = Norm(gradFk) * Norm(gradFk);

            return omega / (Norm(gradFk_1) * Norm(gradFk_1));
        }

        double DichotomyMSG(double[] xk, double[] Sk, double ai, double bi)
        {
            double x1 = 0, x2 = 0;
            double[] tmp1 = new double[xk.Length];
            double[] tmp2 = new double[xk.Length];
            double delta = 1e-7;
            
            double f1 = 0, f2 = 0;

            while (Math.Abs(bi - ai) >= eps)
            {
                x1 = (ai + bi - delta) / 2;
                x2 = (ai + bi + delta) / 2;
                for (int i = 0; i < xk.Length; i++)
                {
                    tmp1[i] = xk[i] + x1 * Sk[i];
                    tmp2[i] = xk[i] + x2 * Sk[i];
                }
                f1 = F.Value(tmp1);
                f2 = F.Value(tmp2);
                if (f1 < f2)
                {
                    // a1 = a0;
                    bi = x2;
                }
                else
                {
                    ai = x1;
                    // b1 = b0;
                }
            }
            if (f1 < f2) return x1;
            else return x2;
        }

        interfaces.FunctionalWithDiff interfaces.Optimizer.functional
        {
            set
            {
                F = value;
            }
        }
        double[] interfaces.Optimizer.Optimize(double[] initial,
            IProgress<(double[] current, double residual, int progresslen, int progressval)> progress)
        {
            int userTime = time.Hour * 60 * 60 * 1000 + time.Minute * 60 * 1000 + time.Second * 1000 + time.Millisecond;
            int timesGone = 0;
            (double[] current, double residual, int progresslen, int progressval) realProgress = (new double[initial.Length], 0.0, 0, 0);
            realProgress.progresslen = userTime;
            System.Diagnostics.Stopwatch clock = new System.Diagnostics.Stopwatch();

            // МСГ метод
            bool isTime = true;
            int n = initial.Length;
            int maxiter = 1000; //10000
            double omega = 0;
            double lambda = 1;
            double[] gradF = new double[n];
            double[] xk = new double[n], xkNext = new double[n];
            double[] S0 = new double[n], Sk = new double[n];

            double[] tmp = new double[n];
            double norm = 1;
            double fk;

            // Диапазоны изменения параметров
            double[] minValue = new double[n], maxValue = new double[n];
            for (int i = 0; i < n; i++)
            {
                minValue[i] = F.Range[i].min;
                maxValue[i] = F.Range[i].max;
                gradF[i] = F.DfDp(i, initial);
                Sk[i] = -gradF[i];
            }

            xk = initial;
            Console.WriteLine(xk[0] + "\t" + xk[1] + "\t it is xk first\n");
            norm = Norm(Sk);
            Console.Write("\nCalculation...\n\n");

            clock.Start();
            while (isTime && timesGone < userTime)
            {
                for (int i = 0; i < maxiter && norm > eps; i++)
                {
                    lambda = DichotomyMSG(xk, Sk, -1000, 1000);
                    for (int j = 0; j < n; j++)
                        xkNext[j] = xk[j] + lambda * Sk[j];

                    fk = F.Value(xkNext);
                    omega = CalcOmega(xkNext, xk);

                    for (int j = 0; j < n; j++)
                    {
                        gradF[j] = F.DfDp(j, xkNext);
                        Sk[j] = Sk[j] * omega - gradF[j];
                    }
                    // Swap(xkNext, xk);
                    // tmp = xk;
                    xk = xkNext;
                    // xkNext = tmp;
                    norm = Norm(Sk);
                    // Swap(Sk, S0);
                    // tmp = S0;
                    // S0 = Sk;
                    // Sk = tmp;

                    realProgress.residual = fk;
                    timesGone = Convert.ToInt32(clock.ElapsedMilliseconds);
                    xk.CopyTo(realProgress.current, 0);
                    realProgress.progressval = timesGone;
                    progress.Report(realProgress);
                }
                timesGone = Convert.ToInt32(clock.ElapsedMilliseconds);
                if (norm < eps) isTime = false;
            }
            clock.Stop();
            Console.Write("Calculation Finished!\n");
            return xkNext;
        }
        double interfaces.Optimizer.Eps
        {
            set 
            {
                eps = value;
            }
        }
        DateTime interfaces.Optimizer.MaxTime
        {
            set
            {
                time = value;
            }
        }
    }
}