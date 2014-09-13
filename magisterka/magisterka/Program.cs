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
            QgAlgorithm alg = new QgAlgorithm();
            UserInterface.RunTest(alg);
        }
    }
}
