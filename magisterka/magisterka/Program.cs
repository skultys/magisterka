using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace magisterka
{
    class Program
    {
        static void Main(string[] args)
        {
            System.Globalization.CultureInfo cultureInfo = new System.Globalization.CultureInfo("en-US");
            QgAlgorithm alg = new QgAlgorithm();
            alg.Execute();
        }
    }
}
