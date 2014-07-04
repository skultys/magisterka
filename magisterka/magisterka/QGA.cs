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
        public int observedState;

        public double Alpha { set; get; }
        public double Beta { set; get; }

        public Qbit()
        {
            this.Alpha = 1.0/Math.Sqrt(2.0);
            this.Beta = 1.0/Math.Sqrt(2.0);
        }

        public Qbit(Qbit anotherQbit)
        {
            this.Alpha = anotherQbit.Alpha;
            this.Beta = anotherQbit.Beta;
            this.observedState = anotherQbit.observedState;
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
                this.observedState = 1;
                return 1;
            }
            else
            {
                this.observedState = 0;
                return 0;
            }
        }

        public void ExecuteRotationGate(Qbit best)
        {
            double theta;
            double alphaTimesBeta = this.Alpha * this.Beta;
            double angle = 0.0;
            int sign = 0;
            if (this.observedState == 1 && best.observedState == 0)
            {
                if (alphaTimesBeta > 0) sign = -1;
                else if (alphaTimesBeta < 0) sign = 1;
                else if (this.Alpha == 0)
                {
                    double d = rand.NextDouble();
                    if (d > 0.5) sign = 1;
                    else sign = -1;
                }
                angle = 0.5 * Math.PI;
            }
            else if (this.observedState == 1 && best.observedState == 1)
            {
                if (alphaTimesBeta > 0) sign = 1;
                else if (alphaTimesBeta < 0) sign = -1;
                else if (this.Beta == 0)
                {
                    double d = rand.NextDouble();
                    if (d > 0.5) sign = 1;
                    else sign = -1;
                }
                angle = 0.2 * Math.PI;
            };

            theta = angle * sign;

            double tempAlpha = this.Alpha;
            this.Alpha = Math.Cos(theta) * this.Alpha - Math.Sin(theta) * this.Beta;
            this.Beta =  Math.Sin(theta) * tempAlpha + Math.Cos(theta) * this.Beta;
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

        public Chromosome(Chromosome anotherChromosome)
        {
            this.genes = new List<Qbit>();

            for (int i = 0; i < anotherChromosome.genes.Count; i++)
            {
                Qbit qbit = new Qbit(anotherChromosome.genes[i]);
                this.genes.Add(qbit);
            }

            this.decodedValue = anotherChromosome.decodedValue;
        }

        public int DecodeChromosome()
        {
            int chromosomeValue = 0;
            for (int i = this.genes.Count - 1; i >= 0; i--)
            {
                chromosomeValue += this[i].EvaluateState() * (int)Math.Pow(2.0, i);
            }
            this.decodedValue = chromosomeValue;
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

        public int Value
        {
            get
            {
                return this.decodedValue;
            }
        }
    }

    class Solution : ISolution
    {
        List<Chromosome> chromosomes = null;

        //zmienic na private
        public List<int> permutation = null;
        double goal;
        int solSize;

        public int Size {
            get
            {
                return this.solSize;
            }
        }

        public Solution(int size)
        {
            this.chromosomes = new List<Chromosome>();
            this.permutation = new List<int>();
            int chromosomeSize = (int)(Math.Log(size, 2.0) + 1);
            this.goal = 0.0;
            this.solSize = size;
            for (int i = 0; i < size; i++)
            {
                Chromosome chromosome = new Chromosome(chromosomeSize);
                this.chromosomes.Add(chromosome);
                this.permutation.Add(0);
            }
            toPermutation();
        }

        public Solution(Solution anotherSolution)
        {
            this.chromosomes = new List<Chromosome>();
            this.permutation = new List<int>();
            for (int i = 0; i < anotherSolution.chromosomes.Count; i++)
            {
                Chromosome chr = new Chromosome(anotherSolution.chromosomes[i]);
                this.chromosomes.Add(chr);
            }

            for (int i = 0; i < anotherSolution.permutation.Count; i++)
            {
                this.permutation.Add(anotherSolution.permutation[i]);
            }

            this.goal = anotherSolution.goal;
            this.solSize = anotherSolution.solSize;
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

        public void toPermutation()
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

            for (int i = 0; i < tupleList.Count; i++)
            {
                tupleList[i] = new Tuple<int, int>(tupleList[i].Item1, i + 1);
            }

            foreach (var tuple in tupleList)
            {
                int val = tuple.Item2;
                this.permutation[tuple.Item1] = val;
            }

            this.goal = 0.0;
            for (int i = 0; i < this.permutation.Count - 1; i++)
            {
                int first = this.permutation[i] - 1;
                int second = this.permutation[i + 1] - 1;
                double flowValue = (QapData.Instance.getFlow())[first, second];
                double distanceValue = (QapData.Instance.getDistance())[first, second];
                this.goal += (flowValue + distanceValue);
            }
        }

        public double Goal
        {
            set
            {
                this.goal = value;
            }
            get
            {
                return this.goal;
            }
        }
    }

    class Population : IPopulation {
        int popSize;
        int problemSize;
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
            this.problemSize = problemSize;
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
                return this.popSize;
            }
        }

        public int ProblemSize
        {
            get
            {
                return this.problemSize;
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
        private static Random rand = new Random();

        public int solSize { get; set; }

        int bitsInSol;

        public double MutationProbability { get; set; }

        public MutationOperator(double probability, int solSize)
        {
            this.MutationProbability = probability;
            this.solSize = solSize;
            this.bitsInSol = (int)(Math.Log(solSize, 2.0) + 1);
        }

        public ISolution Execute(IPopulation population)
        {
            double ifMutate;
            foreach(Solution sol in population)
            {
                ifMutate = rand.NextDouble();
                if (ifMutate <= this.MutationProbability)
                {
                    int selectedChromosome = rand.Next(0, this.solSize - 1);
                    int selectedQbit = rand.Next(0, this.bitsInSol - 1);

                    sol[selectedChromosome][selectedQbit].ExecuteNotGate();
                    sol.toPermutation();
                }
            }
            return new Solution(1);
        }
    }

    class OxCrossoverOperator : ICrossoverOperator
    {
        public double CrossoverProbability { get; set; }

        public OxCrossoverOperator(double probability = 0.0)
        {
            this.CrossoverProbability = probability;
        }

        public ISolution[] Execute(IPopulation population)
        {
            return new Solution[0];
        }
    }


    class PmxCrossoverOperator : ICrossoverOperator
    {
        public double CrossoverProbability { get; set; }

        public PmxCrossoverOperator(double probability = 0.0)
        {
            this.CrossoverProbability = probability;
        }

        public ISolution[] Execute(IPopulation population)
        {
            return new Solution[0];
        }
    }

    class CxCrossoverOperator : ICrossoverOperator
    {
        public double CrossoverProbability { get; set; }

        public CxCrossoverOperator(double probability = 0.0)
        {
            this.CrossoverProbability = probability;
        }

        public ISolution[] Execute(IPopulation population)
        {
            Population pop = population as Population;
            return new Solution[0];
        }

        public ISolution[] Execute(Solution parentOne, Solution parentTwo)
        {
            Solution[] childrenArray = new Solution[2];
            int size = parentOne.Size;
            if (parentOne.Size == 1)
            {
                childrenArray[0] = new Solution(parentOne);
                childrenArray[1] = new Solution(parentTwo);
                return childrenArray;
            }
            else
            {
                List<Tuple<Chromosome, int, Chromosome, int, int>> tupleList = new List<Tuple<Chromosome, int, Chromosome, int, int>>();
                Chromosome[] childOne = new Chromosome[size];
                Chromosome[] childTwo = new Chromosome[size];
                int[] permutationOne = new int[size];
                int[] permutationTwo = new int[size];
                for (int i = 0; i < size; i++)
                {
                    Tuple<Chromosome, int, Chromosome, int, int> tuple = new Tuple<Chromosome, int, Chromosome, int, int>(parentOne[i], parentOne.permutation[i], parentTwo[i], parentTwo.permutation[i], i);
                    tupleList.Add(tuple);
                }
                bool swap = false;
                while (tupleList.Count != 0)
                {
                    int stopCondition = tupleList[0].Item2;
                    int currentValue = tupleList[0].Item2;
                    int tupleIndex;
                    do
                    {
                        tupleIndex = tupleList.FindIndex(t => t.Item2 == currentValue);
                        currentValue = tupleList[tupleIndex].Item4;
                        if (!swap)
                        {
                            childOne[tupleList[tupleIndex].Item5] = new Chromosome(tupleList[tupleIndex].Item1);
                            childTwo[tupleList[tupleIndex].Item5] = new Chromosome(tupleList[tupleIndex].Item3);
                            permutationOne[tupleList[tupleIndex].Item5] = tupleList[tupleIndex].Item2;
                            permutationTwo[tupleList[tupleIndex].Item5] = tupleList[tupleIndex].Item4;
                        }
                        else
                        {
                            childOne[tupleList[tupleIndex].Item5] = new Chromosome(tupleList[tupleIndex].Item3);
                            childTwo[tupleList[tupleIndex].Item5] = new Chromosome(tupleList[tupleIndex].Item1);
                            permutationOne[tupleList[tupleIndex].Item5] = tupleList[tupleIndex].Item4;
                            permutationTwo[tupleList[tupleIndex].Item5] = tupleList[tupleIndex].Item2;
                        }
                        tupleList.RemoveAt(tupleIndex);
                    }
                    while (currentValue != stopCondition);
                    swap = !swap;
                }
                for (int i = 0; i < size; i++)
                {
                    Console.Write(permutationOne[i] + "\t");
                }
                Console.WriteLine();
                Console.WriteLine();
                for (int i = 0; i < size; i++)
                {
                    Console.Write(permutationTwo[i] + "\t");
                }
            }
            return childrenArray;
        }
    }

    class CatastropheOperator : IEvolutionaryOperator
    {
        public double CastastropheProbability { get; set; }

        public void Execute(IPopulation population)
        {
        }
    }

    class RotationGateOperator : IEvolutionaryOperator
    {
        public int solSize { get; set; }

        int bitsInSol;

        public RotationGateOperator(int solSize)
        {
            this.solSize = solSize;
            this.bitsInSol = (int)(Math.Log(solSize, 2.0) + 1);
        }

        public void Execute(IPopulation population, Solution best)
        {
            foreach (Solution sol in population)
            {
                for (int i = 0; i < this.solSize; i++)
                {
                    for (int j = 0; j < this.bitsInSol; j++)
                    {
                        sol[i][j].ExecuteRotationGate(best[i][j]);
                    }
                }
                sol.toPermutation();
            }
        }
    }

    public sealed class QapData
    {
        private static QapData instance = null;
        private double[,] Distance = null;
        private double[,] Flow = null;

        public static QapData Instance
        {
            get
            {
                if (instance == null)
                {
                    instance = new QapData();
                }
                return instance;
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

    class QgAlgorithm : IEvolutionAlgorithm
    {
        PmxCrossoverOperator crossOperator = null;
        //zmienic na private
        public CxCrossoverOperator cxOperator = null;
        public MutationOperator mutOperator = null;
        CatastropheOperator catOperator = null;
        public RotationGateOperator rotOperator = null;
        bool isStopped;
        int iterations;
        int popSize;
        int problemSize;
        public Solution best = null;

        public QgAlgorithm(double[,] distance, double[,] flow, double crossProb, double mutProb, double catProb, int iterations, int popSize, int problemSize)
        {
            crossOperator = new PmxCrossoverOperator(crossProb);

            cxOperator = new CxCrossoverOperator();

            mutOperator = new MutationOperator(mutProb, problemSize);

            rotOperator = new RotationGateOperator(problemSize);

            catOperator = new CatastropheOperator();
            catOperator.CastastropheProbability = catProb;

            this.isStopped = false;
            this.iterations = iterations;
            this.popSize = popSize;
            this.problemSize = problemSize;

            QapData.Instance.setQapData(distance, flow);
        }

        public IPopulation Population { get; set; }

        public double EvaluateSolution(ISolution solution)
        {
            return solution.Goal;
        }

        public void InitRandomPopulation()
        {
            this.Population = new Population(this.popSize, this.problemSize);
        }

        public void Execute()
        {
            for (int i = 0; i < this.iterations; i++)
            {
            }
            this.isStopped = true;
        }

        public bool Stop()
        {
            return this.isStopped;
        }

        public ISolution GetBestSolution()
        {
            double currentGoal = -1.0;
            foreach (Solution sol in this.Population)
            {
                if (currentGoal == -1 || sol.Goal < currentGoal)
                {
                    currentGoal = sol.Goal;
                    this.best = new Solution(sol);
                }
            }
            return this.best;
        }
    }
}
