using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace magisterka
{
    class Program
    {
        static void Main(string[] args)
        {
            QapData d = new QapData(2);
            d.Distance[0][0] = new double();
            d.Distance[0][1] = new double();
            d.Distance[1][0] = new double();
            d.Distance[1][1] = new double();
            Console.WriteLine(d.Distance[0][0]);
            Console.WriteLine(d.Distance[0][1]);
            Console.WriteLine(d.Distance[1][0]);
            Console.WriteLine(d.Distance[1][1]);
            Console.ReadKey();
        }
    }
}
