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
            Random rand = new Random();

            int psize = 15;

            double[,] distance = new double[psize, psize];

            double[,] flow = new double[psize, psize];

            for (int i = 0; i < psize; i++)
            {
                for (int j = 0; j < psize; j ++)
                {
                    distance[i, j] = rand.Next(psize); //rand.NextDouble();
                    flow[i, j] = rand.Next(psize); //rand.NextDouble();
                }
            }

            QgAlgorithm alg = new QgAlgorithm(distance, flow, 1.0, 0.05, 0.02, 1000, 20, psize);
            alg.Execute();

            Console.ReadKey();
        }
    }
}
