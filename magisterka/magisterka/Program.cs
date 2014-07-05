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

            int psize = 10;

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

            ((Solution)alg.Population[0]).PrintSolution();
            ((Solution)alg.Population[1]).PrintSolution();


            var OxResult = alg.oxOperator.Execute((Solution)alg.Population[0], (Solution)alg.Population[1]);

            //((Solution)OxResult[0]).PrintSolution();
            //((Solution)OxResult[1]).PrintSolution();
            //Console.WriteLine(((Solution)alg.Population[0]).Goal);
            //Console.WriteLine(((Solution)alg.Population[1]).Goal);

            //Console.WriteLine(OxResult[0].Goal);
            //Console.WriteLine(OxResult[1].Goal);

            Console.ReadKey();
        }
    }
}
