using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
namespace Functionals
{
    class EasyFunction : interfaces.FunctionalWithDiff
    {
        // Параметры
        private double x;
        private double y;
        private double z;
        private double t;

        double interfaces.Functional.Value(double[] parameters)
        {
            x = parameters[0];
            y = parameters[1];
            //z = parameters[2];
            //t = parameters[3];
            return x * x + 3 * y * y + 3 * x * y; // Заменить sinx
            //return x * x + 2 * y * y + 3 * z * z + 4 * t * t;
        }
        (double min, double max)[] interfaces.Functional.Range
        {
            get
            {
                (double min, double max)[] res = new (double min, double max)[2];
                for (int i = 0; i < 2; i++)
                    res[i] = (-10, 15);
                return res;
            }
        }
        double interfaces.FunctionalWithDiff.DfDp(int k, double[] parameters)
        {
            x = parameters[0];
            y = parameters[1];
            //z = parameters[2];
            //t = parameters[3];

            //if (k == 0) return 2 * x; //dx
            //else if (k == 1) return 4 * y; //dy
            //else if (k == 2) return 6 * z; //dz
            //else return 8 * t; //dt

            if (k == 0) return 2 * x + 3 * y; //dx
            else return 6 * y + 3 * x; //dy

        }
    }

    class Spline : basis.Basis, interfaces.FunctionalWithDiff
    {
        private double[] F;
        private int n;
        // количество шагов для одного к.э.
        private int m = 10;
        // узлы сетки
        private double[] xk;
        // значения функции в узлах сетки
        private double[] fk;
        // значения производных в узлах сетки
        private double[] dfk;
        private double minF = double.MaxValue, maxF = double.MinValue;

        public Spline()
        {
            using (StreamReader file = new StreamReader("xy.txt"))
            {
                string[] words = file.ReadToEnd().Split(new char[] { ' ', '\n', '\r' }, StringSplitOptions.RemoveEmptyEntries);
                n = Convert.ToInt32(words[0]);
                xk = new double[n];
                fk = new double[n];
                F = new double[n];
                dfk = new double[n];
                int id = 1;
                for (; id <= n; id++)
                    xk[id - 1] = Convert.ToDouble(words[id]);
                for (int i = 0; i < n; i++)
                {
                    F[i] = Convert.ToDouble(words[id]);
                    if (F[i] < minF) minF = F[i];
                    if (F[i] > maxF) maxF = F[i];
                    id++;
                }
            }
        }
        // Вычисление производных fi
        private void CalcDerivatives()
        {
            double hi, hi1;
            double temp1, temp2, temp3;
            for (int i = 1; i < n - 1; i++)
            {
                hi = xk[i + 1] - xk[i];
                hi1 = xk[i] - xk[i - 1];
                temp1 = -hi / (hi1 * (hi1 + hi)) * fk[i - 1];
                temp2 = (hi - hi1) / (hi1 * hi) * fk[i];
                temp3 = hi1 / (hi * (hi1 + hi)) * fk[i + 1];
                dfk[i] = temp1 + temp2 + temp3;
            }

            hi1 = xk[1] - xk[0]; //h0
            hi = xk[2] - xk[1]; //h1

            temp1 = -(2 * hi1 + hi) / (hi1 * (hi1 + hi)) * fk[0];
            temp2 = (hi1 + hi) / (hi1 * hi) * fk[1];
            temp3 = -hi1 / (hi * (hi1 + hi)) * fk[2];
            dfk[0] = temp1 + temp2 + temp3;

            hi = xk[n - 1] - xk[n - 2]; //hn-2
            hi1 = xk[n - 2] - xk[n - 3]; //hn-3
            temp1 = hi / (hi1 * (hi1 + hi)) * fk[n - 3];
            temp2 = -(hi1 + hi) / (hi1 * hi) * fk[n - 2];
            temp3 = (2 * hi + hi1) / (hi * (hi1 + hi)) * fk[n - 1];
            dfk[n - 1] = temp1 + temp2 + temp3;
        }
        double interfaces.Functional.Value(double[] parameters)
        {
            for (int i = 0; i < parameters.Length; i++)
                fk[i] = parameters[i];

            double h, ksi;
            double[] xk_new, fk_new;
            int id;

            xk_new = new double[m * (n - 1) + 1];
            fk_new = new double[m * (n - 1) + 1];

            CalcDerivatives();
            id = 0;
            for (int k = 0; k < n - 1; k++)
            {
                h = (xk[k + 1] - xk[k]) / m;
                for (int i = 0; i < m; i++)
                {
                    xk_new[id] = xk[k] + i * h;
                    ksi = (xk_new[id] - xk[k]) / (xk[k + 1] - xk[k]);
                    fk_new[id] = fk[k] * psi_[0](ksi, h) + dfk[k] * psi_[1](ksi, h) +
                        fk[k + 1] * psi_[2](ksi, h) + dfk[k + 1] * psi_[3](ksi, h);
                    id++;
                }
            }
            // для последнего узла сетки
            fk_new[id] = fk[n - 1];
            xk_new[id] = xk[n - 1];

            double sum = 0;
            for (int i = 0; i < n; i++)
                sum += Math.Pow(F[i] - fk_new[i * m], 2);
            sum = Math.Sqrt(sum);
            return sum;
        }
        (double min, double max)[] interfaces.Functional.Range
        {
            get
            {
                (double min, double max)[] res = new(double min, double max)[n];
                for (int i = 0; i < n; i++)
                    res[i] = (minF, maxF);
                return res;
            }
        }
        double interfaces.FunctionalWithDiff.DfDp(int j, double[] parameters)
        {
            for (int i = 0; i < parameters.Length; i++)
                fk[i] = parameters[i];

            double h, ksi;
            double[] xk_new, fk_new;
            int id;

            xk_new = new double[m * (n - 1) + 1];
            fk_new = new double[m * (n - 1) + 1];

            CalcDerivatives();
            id = 0;
            for (int k = 0; k < n - 1; k++)
            {
                h = (xk[k + 1] - xk[k]) / m;
                for (int i = 0; i < m; i++)
                {
                    xk_new[id] = xk[k] + i * h;
                    ksi = (xk_new[id] - xk[k]) / (xk[k + 1] - xk[k]);
                    fk_new[id] = fk[k] * psi_[0](ksi, h) + dfk[k] * psi_[1](ksi, h) +
                        fk[k + 1] * psi_[2](ksi, h) + dfk[k + 1] * psi_[3](ksi, h);
                    id++;
                }
            }
            // для последнего узла сетки
            fk_new[id] = fk[n - 1];
            xk_new[id] = xk[n - 1];

            double sum = 0;
            for (int i = 0; i < n; i++)
                sum += Math.Pow(F[i] - fk_new[i * m], 2);
            sum = 1 / Math.Sqrt(sum) * (fk_new[j * m] - F[j]);
            return sum;
        }
    }
    delegate double Func(double x1, double x2);
    class LeastSquares : interfaces.FunctionalWithDiff
    {
        // Число измерений
        private int n = 1000;
        // Число параметров
        private int m = 5;
        // Факторы
        private double[] x1;
        private double[] x2;
        // Ошибки
        private double[] eps;
        // Отклики
        private double[] y;
        private double[] u;
        private double[] f;
        private double[] tetta;
        private double[] minTetta, maxTetta;
        private Func[] F = new Func[5];
        private double[,] X;

        public LeastSquares()
        {
            // Регрессоры
            F[0] = (x1, x2) => 1;
            F[1] = (x1, x2) => x1;
            F[2] = (x1, x2) => x2;
            F[3] = (x1, x2) => x1 * x1;
            F[4] = (x1, x2) => x2 * x2 * x2;
            Random rand = new Random();
            x1 = new double[n];
            x2 = new double[n];
            u = new double[n];
            y = new double[n];
            X = new double[n, m];
            minTetta = new double[m];
            maxTetta = new double[m];
            tetta = new double[m];
            eps = new double[n];
            for (int i = 0; i < n; i++)
            {
                x1[i] = rand.NextDouble() * 2 - 1;
                x2[i] = rand.NextDouble() * 2 - 1;
                // Матрица регрессоров
                for (int j = 0; j < m; j++)
                    X[i, j] = F[j](x1[i], x2[i]);
            }
            for (int i = 0; i < m; i++)
            {
                minTetta[i] = -1;
                maxTetta[i] = 1;
            }

            // Истинные параметры
            double[] t = { 1, 0.5, -1, 0.25, -1 };
            double u_average = 0;
            for (int i = 0; i < n; i++)
            {
                u[i] = 0;
                for (int j = 0; j < m; j++)
                    u[i] += t[j] * F[j](x1[i], x2[i]);
                u_average += u[i];
            }
            u_average /= n;
            // Мощность сигнала
            double w = 0;
            for (int i = 0; i < n; i++)
                w += Math.Pow(u[i] - u_average, 2);
            w /= n - 1;
            // Дисперсия (уже в квадрате)
            double sigma = 0.5 * w;
            // Генерация вектора ошибок (нормальное распределение)
            Random ran = new Random();
            double dSumm = 0;
            for (int i = 0; i < n; i++)
            {
                dSumm = 0;
                for (int j = 0; j <= 12; j++)
                {
                    double R = ran.NextDouble();
                    dSumm = dSumm + R;
                }
                // т.к. сгенерированы числа с m = n/2, sigma = sqrt(n/12), n = 12,
                // то сначала переводим числа в N(0,1) ((dSumm - 6)/1)
                // а потом в N(0, sigma^2) формулой x = z * sigma^2 + m, m = 0
                eps[i] = Math.Round((sigma * (dSumm - 6)), 3);
                // Выходные данные
                //y[i] = eps[i] + u[i];
                y[i] = u[i];
            }
        }
        public void GaussMethod(double[] massive, double mu, double sigma, int num)
        {
            double dSumm = 0, dRandValue = 0;
            Random ran = new Random();
            for (int n = 0; n <= num; n++)
            {
                dSumm = 0;
                for (int i = 0; i <= 12; i++)
                {
                    double R = ran.NextDouble();
                    dSumm = dSumm + R;
                }
                dRandValue = Math.Round((mu + sigma * (dSumm - 6)), 3);
                massive[n] = dRandValue;
            }

        }
        double interfaces.Functional.Value(double[] parameters)
        {
            for (int i = 0; i < m; i++)
                tetta[i] = parameters[i];
            double sum = 0;
            for (int i = 0; i < n; i++)
            {
                u[i] = 0;
                for (int j = 0; j < m; j++)
                    u[i] += tetta[j] * F[j](x1[i], x2[i]);
                sum += Math.Pow(y[i] - u[i], 2);
            }
            return sum;
        }
        (double min, double max)[] interfaces.Functional.Range
        {
            get
            {
                (double min, double max)[] res = new(double min, double max)[m];
                for (int i = 0; i < m; i++)
                    res[i] = (minTetta[i], maxTetta[i]);
                return res;
            }
        }
        double interfaces.FunctionalWithDiff.DfDp(int k, double[] parameters)
        {
            double sum = 0;
            for (int i = 0; i < n; i++)
            {
                u[i] = 0;
                for (int j = 0; j < m; j++)
                    u[i] += tetta[j] * F[j](x1[i], x2[i]);
                sum += 2 * (u[i] - y[i]) * F[k](x1[i], x2[i]);
            };
            return sum;
        }
    }
}