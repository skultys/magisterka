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

            int psize = 2;

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
            alg.GetBestSolution();

            Solution best = alg.best;

            int bits = (int)(Math.Log(psize, 2.0) + 1);

            foreach (Solution sol in alg.Population)
            {
                for (int i = 0; i < psize; i++)
                {
                    /*for (int j = 0; j < bits; j++)
                    {
                        Console.Write(sol[i][j].observedState);
                    }*/
                    Console.Write(sol.permutation[i]);
                    Console.Write(" ");
                }
                Console.Write("\t goal: " + sol.Goal);
                Console.WriteLine();
            }

            Console.WriteLine("BEST");

            for (int i = 0; i < psize; i++)
            {
                /*for (int j = 0; j < bits; j++)
                {
                    //Console.WriteLine(best[i][j].Alpha + " " + best[i][j].Beta);
                    Console.Write(best[i][j].observedState);
                }*/
                Console.Write(best.permutation[i]);
                Console.Write(" ");
            }
            Console.Write("\t goal: " + best.Goal);

            Console.WriteLine();
            Console.WriteLine();

            alg.rotOperator.Execute(alg.Population, best);

            foreach (Solution sol in alg.Population)
            {
                for (int i = 0; i < psize; i++)
                {
                    /*for (int j = 0; j < bits; j++)
                    {
                        //Console.WriteLine(sol[i][j].Alpha + " " + sol[i][j].Beta);
                        Console.Write(sol[i][j].observedState);
                    }*/
                    Console.Write(sol.permutation[i]);
                    Console.Write(" ");
                }
                Console.Write("\t goal: " + sol.Goal);
                Console.WriteLine();
            }

            Console.ReadKey();
        }
    }
}
