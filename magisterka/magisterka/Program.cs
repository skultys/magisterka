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
            double[,] flow = { { 1.0, 2.0 }, { 3.0, 4.0 } };
            double[,] flow2 = new double[2,2];
            for (int i = 0; i < 2; i++)
            {
                for (int j = 0; j < 2; j++)
                {
                    flow2[i, j] = flow[i, j];
                }
            }
            flow2[1, 1] = 10;

            QapData.Instance.setQapData(flow, flow2);

            int psize = 10;

            Population pop = new Population(1, psize);

            Console.ReadKey();
        }
    }
}
