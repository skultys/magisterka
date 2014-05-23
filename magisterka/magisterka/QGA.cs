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
        int qbitState;
        public double Alpha { set; get; }
        public double Beta { set; get; }

        public Qbit()
        {
            qbitState = 0;
            this.Alpha = 1.0/Math.Sqrt(2.0);
            this.Beta = 1.0/Math.Sqrt(2.0);
        }

        public Qbit(double alpha, double beta)
        {
            qbitState = 0;
            this.Alpha = alpha;
            this.Beta = beta;
        }

        public int EvaluateState(double treshold)
        {
            double alphaSquare = this.Alpha * this.Alpha;
            if (alphaSquare > treshold)
            {
                qbitState = 1;
                return 1;
            }
            else
            {
                qbitState = 0;
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
        List<Qbit> genes;
        int decodedValue;

        public Chromosome()
        {
            genes = new List<Qbit>();
            decodedValue = 0;
        }
        
        public Chromosome(int size)
        {
            genes = new List<Qbit>(size);
            decodedValue = 0;
        }

        public int DecodeChromosome(double treshold)
        {
            int chromosomeValue = 0;
            for (int i = 0; i < this.genes.Count; i++)
            {
                chromosomeValue += this[i].EvaluateState(treshold) * (int)Math.Pow(2.0, i);
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
        List<int> orderList;
        List<Chromosome> chromosomeList;

        public Solution()
        {
            orderList = new List<int>();
            List<Chromosome> chromosomeList = new List<Chromosome>();
        }

        public Solution(int size)
        {
            orderList = new List<int>(size);
            int chromosomeSize = (int)(Math.Log(size, 2.0) + 1);
            for (int i = 0; i < size; i++)
            {
                Chromosome chromosome = new Chromosome(chromosomeSize);
                this.chromosomeList.Add(chromosome);
            }
        }

        public Chromosome this[int index]
        {
            set
            {
                this.chromosomeList[index] = value;
            }
            get
            {
                return this.chromosomeList[index];
            }
        }

        public double Goal { set; get; }
    }

    class Population : IPopulation {
        int popSize;
        int currPopSize;
        List<Solution> solutions;

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
            return new Solution();
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
            return new Solution();
        }
    }

    class QapData
    {
        int size;

        public double[][] Distance;
        public double[][] Flow;

        public QapData(int QapSize)
        {
            this.size = QapSize;
            this.Distance = new double[QapSize][];
            this.Flow = new double[QapSize][];
            for (int i = 0; i < QapSize; i++)
            {
                this.Distance[i] = new double[QapSize];
                this.Flow[i] = new double[QapSize];
            }
        }
    }

    class QGA : IEvolutionAlgorithm
    {
        CrossoverOperator crossOperator;
        MutationOperator mutOperator;
        CatastropheOperator catOperator;
        QapData qapData;

        public QGA(QapData qapData)
        {
            crossOperator = new CrossoverOperator();
            mutOperator = new MutationOperator();
            catOperator = new CatastropheOperator();
            this.qapData = qapData;
        }

        public IPopulation Population { get; set; }

        public double EvaluateSolution(ISolution solution)
        {
            double solValue = 0.0;
            for (int i = 0; i < this.qapData.Distance.Length - 1; i++)
            {
                //solValue += solution
            }
            return solValue;
        }

        public void InitRandomPopulation()
        {
        }

        public void Execute()
        {
        }

        public bool Stop()
        {
            return false;
        }

        public ISolution GetBestSolution()
        {
            return new Solution(0);
        }
    }
}
