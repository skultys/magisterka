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
                    distance[i, j] = rand.NextDouble();
                    flow[i, j] = rand.NextDouble();
                }
            }

            QgAlgorithm alg = new QgAlgorithm(distance, flow, 0.2, 0.5, 0.02, 1000, 3, psize);
            alg.InitRandomPopulation();

            foreach (Solution sol in alg.Population)
            {
                sol.PrintSolution();
            }

            Console.WriteLine();

            alg.cxOperator.Execute((Solution)alg.Population[0], (Solution)alg.Population[1]);

            Console.ReadKey();
        }
    }
}
