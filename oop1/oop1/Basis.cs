using System.Collections.Generic;
using System;
namespace basis
{
    delegate double SplineFunc1D(double ksi, double h);
    class Basis 
    {
        public SplineFunc1D[] psi_ = new SplineFunc1D[4];
        public Basis()
        {
            psi_[0] = (ksi, h) => 1 - 3 * Math.Pow(ksi, 2) + 2 * Math.Pow(ksi, 3);
            psi_[1] = (ksi, h) => h * (ksi - 2 * Math.Pow(ksi, 2) + Math.Pow(ksi, 3));
            psi_[2] = (ksi, h) => 3 * Math.Pow(ksi, 2) - 2 * Math.Pow(ksi, 3);
            psi_[3] = (ksi, h) => h * (-Math.Pow(ksi, 2) + Math.Pow(ksi, 3));
        }
    }        
}