using System;

namespace interfaces
{
    interface Functional
    {
        // Функция  возвращающая значение функционала
        double Value (double[] parameters);
        // Диапазоны изменения параметров
        (double min, double max)[] Range { get; }
    }
    interface FunctionalWithDiff:Functional
    {//
        double DfDp(int i, double[] parameters);
    }
    
    interface Optimizer
    {
        FunctionalWithDiff functional { set; }
        double[] Optimize(double[] initial,
            IProgress<(double[] current,double residual,int progresslen,int progressval )> progress);
        double Eps { set;}
        DateTime MaxTime { set; }
    }

}

