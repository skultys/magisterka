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

            QgAlgorithm alg = new QgAlgorithm(distance, flow, 0.2, 0.5, 0.02, 1000, 5, psize);
            alg.InitRandomPopulation();

            foreach (Solution sol in alg.Population)
            {
                Console.WriteLine(sol.Goal);
            }

            /*((Population)alg.Population).SortAscending();

            Console.WriteLine();
            foreach (Solution sol in alg.Population)
            {
                Console.WriteLine(sol.Goal);
            }

            ((Population)alg.Population).SortDescending();

            Console.WriteLine();
            foreach (Solution sol in alg.Population)
            {
                Console.WriteLine(sol.Goal);
            }
            //((Solution)alg.Population[0]).PrintSolution();
            //((Solution)alg.Population[1]).PrintSolution();*/


            //var pmxResult = alg.pmxOperator.Execute((Solution)alg.Population[0], (Solution)alg.Population[1]);

            //((Solution)pmxResult[0]).PrintSolution();
            //((Solution)pmxResult[1]).PrintSolution();
            //Console.WriteLine(((Solution)alg.Population[0]).Goal);
            //Console.WriteLine(((Solution)alg.Population[1]).Goal);

            //Console.WriteLine(OxResult[0].Goal);
            //Console.WriteLine(OxResult[1].Goal);

            Console.ReadKey();
        }
    }
}
