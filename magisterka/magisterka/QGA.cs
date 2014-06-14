using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using OptimisationClassLibrary;
using System.Collections;

namespace magisterka
{
    class Qbit
    {
        private static Random rand = new Random();

        public double Alpha { set; get; }
        public double Beta { set; get; }

        public Qbit()
        {
            this.Alpha = 1.0/Math.Sqrt(2.0);
            this.Beta = 1.0/Math.Sqrt(2.0);
        }

        public Qbit(double alpha, double beta)
        {
            this.Alpha = alpha;
            this.Beta = beta;
        }

        public int EvaluateState()
        {
            double treshold = rand.NextDouble();
            double alphaSquare = this.Alpha * this.Alpha;
            if (alphaSquare > treshold)
            {
                return 1;
            }
            else
            {
                return 0;
            }
        }

        public void ExecuteRotationGate(double theta)
        {
            this.Alpha = Math.Cos(theta) * this.Alpha - Math.Sin(theta) * this.Beta;
            this.Beta = Math.Sin(theta) * this.Alpha + Math.Cos(theta) * this.Beta;
        }

        public void ExecuteNotGate()
        {
            double temp = this.Alpha;
            this.Alpha = this.Beta;
            this.Beta = temp;
        }
    }

    class Chromosome
    {
        List<Qbit> genes = null;
        int decodedValue;

        public Chromosome()
        {
            genes = new List<Qbit>();
            decodedValue = 0;
        }
        
        public Chromosome(int size)
        {
            genes = new List<Qbit>();
            decodedValue = 0;
            for (int i = 0; i < size; i++)
            {
                Qbit qbit = new Qbit();
                this.genes.Add(qbit);
            }
        }

        public int DecodeChromosome()
        {
            int chromosomeValue = 0;
            for (int i = this.genes.Count - 1; i >= 0; i--)
            {
                chromosomeValue += this[i].EvaluateState() * (int)Math.Pow(2.0, i);
            }
            decodedValue = chromosomeValue;
            return chromosomeValue;
        }

        public Qbit this[int index]
        {
            set
            {
                this.genes[index] = value;
            }
            get
            {
                return this.genes[index];
            }
        }
    }

    class Solution : ISolution
    {
        List<Chromosome> chromosomes = null;
        List<int> permutation = null;

        public Solution(int size)
        {
            this.chromosomes = new List<Chromosome>();
            this.permutation = new List<int>();
            int chromosomeSize = (int)(Math.Log(size, 2.0) + 1);
            for (int i = 0; i < size; i++)
            {
                Chromosome chromosome = new Chromosome(chromosomeSize);
                this.chromosomes.Add(chromosome);
                this.permutation.Add(0);
            }
            toPermutation();
        }

        public Chromosome this[int index]
        {
            set
            {
                this.chromosomes[index] = value;
            }
            get
            {
                return this.chromosomes[index];
            }
        }

        private void toPermutation()
        {
            List<Tuple<int, int> > tupleList = new List<Tuple<int, int> >();

            for (int i = 0; i < this.chromosomes.Count; i++)
            {
                Tuple<int, int> tuple = new Tuple<int, int>(i, this.chromosomes[i].DecodeChromosome());
                tupleList.Add(tuple);
            }

            tupleList.Sort((a, b) =>
                {
                    if (a.Item2 < b.Item2) return -1;
                    else
                    {
                        if (a.Item2 == b.Item2 && a.Item1 < b.Item1) return -1;
                        else return 1;
                    }
                }
            );

            foreach (var el in tupleList)
            {
                Console.WriteLine(el);
            }
            Console.WriteLine();

            for (int i = 0; i < tupleList.Count; i++)
            {
                tupleList[i] = new Tuple<int, int>(tupleList[i].Item1, i + 1);
            }

            foreach (var tuple in tupleList)
            {
                int val = tuple.Item2;
                this.permutation[tuple.Item1] = val;
            }

            foreach (int val in this.permutation)
            {
                Console.WriteLine(val);
            }

            Console.WriteLine();
            Console.WriteLine();
        }

        public double Goal { set; get; }
    }

    class Population : IPopulation {
        int popSize;
        int currPopSize;
        List<Solution> solutions = null;

        public Population()
        {
            popSize = currPopSize = 0;
            solutions = new List<Solution>();
        }

        public Population(int solutionsCount, int problemSize)
        {
            solutions = new List<Solution>();
            this.popSize = this.currPopSize = solutionsCount;
            for (int i = 0; i < solutionsCount; i++)
            {
                Solution solution = new Solution(problemSize);
                solutions.Add(solution);
            }
        }

        public int Size
        {
            get
            {
                return popSize;
            }
        }

        public int Add(ISolution sol)
        {
            return 0;
        }

        public ISolution this[int index]
        {
            get
            {
                return this.solutions[index];
            }
            set
            {
                this.solutions[index] = (Solution)value;
            }
        }

        public int CurrentSize
        {
            get
            {
                return currPopSize;
            }
        }

        public void Remove(int idx)
        {
            this.solutions.RemoveAt(idx);
        }

        public IEnumerator GetEnumerator()
        {
            foreach (Solution solution in this.solutions)
            {
                if (solution == null)
                {
                    break;
                }
                yield return solution;
            }
        }
    }

    class MutationOperator : IMutationOperator
    {
        public ISolution Execute(IPopulation population)
        {
            return new Solution(1);
        }
    }

    class CrossoverOperator : ICrossoverOperator
    {
        public ISolution[] Execute(IPopulation population)
        {
            return new Solution[0];
        }
    }

    class CatastropheOperator : IEvolutionaryOperator
    {
        public ISolution Execute(IPopulation population)
        {
            return new Solution(1);
        }
    }

    public sealed class QapData
    {
        private static QapData m_oInstance = null;
        private double[,] Distance = null;
        private double[,] Flow = null;

        public static QapData Instance
        {
            get
            {
                if (m_oInstance == null)
                {
                    m_oInstance = new QapData();
                }
                return m_oInstance;
            }
        }

        public void setQapData(double[,] distance, double[,] flow)
        {
            this.Distance = distance;
            this.Flow = flow;
        }

        public double[,] getDistance()
        {
            return this.Distance;
        }

        public double[,] getFlow()
        {
            return this.Flow;
        }

        private QapData()
        {
        }
    }

    class QGA : IEvolutionAlgorithm
    {
        CrossoverOperator crossOperator;
        MutationOperator mutOperator;
        CatastropheOperator catOperator;
        QapData qapData;
        bool isStopped;

        public QGA(QapData qapData)
        {
            crossOperator = new CrossoverOperator();
            mutOperator = new MutationOperator();
            catOperator = new CatastropheOperator();
            this.qapData = qapData;
            this.isStopped = false;
        }

        public IPopulation Population { get; set; }

        public double EvaluateSolution(ISolution solution)
        {
            return solution.Goal;
        }

        public void InitRandomPopulation()
        {
        }

        public void Execute()
        {
        }

        public bool Stop()
        {
            return this.isStopped;
        }

        public ISolution GetBestSolution()
        {
            return new Solution(0);
        }
    }
}
