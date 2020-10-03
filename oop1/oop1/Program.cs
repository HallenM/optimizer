using System;
using System.Collections.Generic;
using System.IO;

namespace oop1
{
    class Program
    {
        enum Method { SIMPLE_SEARCH, GAUSS, MSG};
        enum Functional { EASY_F, SPLINE, MNK };
        static void Main(string[] args)
        {
            Method m = Method.SIMPLE_SEARCH;
            Functional f = Functional.SPLINE;

            interfaces.Optimizer s = null;
            interfaces.FunctionalWithDiff F = null;
            optimizer.CurrentProgress pr = new optimizer.CurrentProgress();

            switch (m)
            {
                case Method.SIMPLE_SEARCH: s = new optimizer.SimpleRandomSearch(); break;
                case Method.GAUSS: s = new optimizer.Gauss(); break;
                case Method.MSG: s = new optimizer.MSG(); break;
            }

            switch (f)
            {
                case Functional.EASY_F: F = new Functionals.EasyFunction(); break;
                case Functional.SPLINE: F = new Functionals.Spline(); break;
                case Functional.MNK: F = new Functionals.LeastSquares(); break;
            }

            s.functional = F;
            s.Eps = 1e-6;
            DateTime userTime = new DateTime(); // год/месяц/день/час/минута/секунда
            userTime = DateTime.Now.Date.Add(new TimeSpan(0, 0, 120)); // час / минута / секунда
            s.MaxTime = userTime;

            int len = F.Range.Length;
            double[] initial = new double[len];
            double[] resultParameters;
            for (int i = 0; i < len; i++)
                initial[i] = (F.Range[i].max + F.Range[i].min) / 2;
            resultParameters = s.Optimize(initial, pr);
            double val = F.Value(resultParameters);

            Console.Write("\nLast value\nResidual: \t" + val.ToString() + "\n\nParameters value:\n");
            for (int i = 0; i < resultParameters.Length; i++)
                Console.Write(resultParameters[i].ToString() + "\n");
            Console.Write("Для продолжения нажмите Enter...\n");
            Console.ReadLine();

            using (StreamWriter file = new StreamWriter("result.txt", true))
            {
                file.WriteLine();
                file.WriteLine("Residual: \t" + val.ToString());
                file.WriteLine("Parameters value:");
                for (int i = 0; i < resultParameters.Length; i++)
                    file.WriteLine(resultParameters[i].ToString());
                file.WriteLine("=============================================================\n");
            }

        }
    }
}
