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
        private List<Qbit> genes = null;
        private int decodedValue;

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
        private List<Chromosome> chromosomes = null;
        public List<int> permutation = null;
        int solSize;

        public int Size
        {
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
            this.Goal = 0.0;
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

            this.Goal = anotherSolution.Goal;
            this.solSize = anotherSolution.solSize;
        }

        public Solution(Chromosome[] chromosomes)
        {
            this.solSize = chromosomes.Length;
            this.Goal = 0.0;
            this.chromosomes = new List<Chromosome>(chromosomes);
            this.permutation = new List<int>();
            for (int i = 0; i < this.Size; i++)
            {
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

            this.Goal = 0.0;
            for (int i = 0; i < this.permutation.Count - 1; i++)
            {
                int first = this.permutation[i] - 1;
                int second = this.permutation[i + 1] - 1;
                double flowValue = (QapData.Instance.getFlow())[first, second];
                double distanceValue = (QapData.Instance.getDistance())[first, second];
                this.Goal += (flowValue + distanceValue);
            }
        }

        public double Goal { get; set; }

        public void PrintSolution()
        {
            foreach (int i in this.permutation)
            {
                Console.Write(i + "  ");
            }
            Console.WriteLine();
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
            this.solutions = new List<Solution>();
            this.popSize = this.currPopSize = solutionsCount;
            this.problemSize = problemSize;
            for (int i = 0; i < solutionsCount; i++)
            {
                Solution solution = new Solution(problemSize);
                this.solutions.Add(solution);
            }
        }

        public Population(Solution[] solutions)
        {
            this.solutions = new List<Solution>(solutions);
            this.popSize = this.currPopSize = solutions.Length;
            this.problemSize = solutions[0].Size;
        }

        public void SortAscending()
        {
            this.solutions.Sort((a, b) => a.Goal.CompareTo(b.Goal));
        }

        public void SortDescending()
        {
            this.solutions.Sort((a, b) => b.Goal.CompareTo(a.Goal));
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
            Solution solToReturn = null;
            foreach(Solution sol in population)
            {
                ifMutate = rand.NextDouble();
                if (ifMutate <= this.MutationProbability)
                {
                    int selectedChromosome = rand.Next(0, this.solSize - 1);
                    int selectedQbit = rand.Next(0, this.bitsInSol - 1);

                    sol[selectedChromosome][selectedQbit].ExecuteNotGate();
                    sol.toPermutation();
                    solToReturn = new Solution(sol);
                    break;
                }
            }
            return solToReturn;
        }
    }

    class OxCrossoverOperator : ICrossoverOperator
    {
        public double CrossoverProbability { get; set; }

        private static Random rand = new Random();

        public OxCrossoverOperator(double probability = 0.0)
        {
            this.CrossoverProbability = probability;
        }

        public ISolution[] Execute(IPopulation population)
        {
            Solution[] crossedPopulation = new Solution[population.Size];
            List<Solution> chosenSolutions = new List<Solution>();

            foreach (Solution sol in population) {
                double ifChosen = rand.NextDouble();
                if(ifChosen < this.CrossoverProbability)
                {
                    Solution tempSol = new Solution(sol);
                    chosenSolutions.Add(tempSol);
                }
            }

            if (chosenSolutions.Count != population.Size)
            {
                int stop = population.Size - chosenSolutions.Count;
                for (int i = 0; i < stop; i++)
                {
                    int whichSol = rand.Next(chosenSolutions.Count - 1);
                    int position = rand.Next(chosenSolutions.Count - 1);
                    Solution tempSol = new Solution(chosenSolutions[whichSol]);
                    chosenSolutions.Insert(position, tempSol);
                }
            }

            for (int i = 0; i < population.Size; i += 2)
            {
                Solution[] temp = new Solution[2];
                temp = Execute(chosenSolutions[i], chosenSolutions[i + 1]);
                crossedPopulation[i] = new Solution(temp[0]);
                crossedPopulation[i + 1] = new Solution(temp[1]);
            }
            return crossedPopulation;
        }

        public Solution[] Execute(Solution parentOne, Solution parentTwo)
        {
            int size = parentOne.Size;

            //int[] permutationOne = new int[size];
            //int[] permutationTwo = new int[size];

            int leftBound = rand.Next(size - 1);
            int rightBound = rand.Next(size - 1);
            Solution[] childrenArray = new Solution[2];
            if (leftBound > rightBound)
            {
                int temp = rightBound;
                rightBound = leftBound;
                leftBound = temp;
            }

            //Console.WriteLine("left " + leftBound);
            //Console.WriteLine("right " + rightBound);
            //Console.WriteLine();

            Chromosome[] childOne = new Chromosome[size];
            Chromosome[] childTwo = new Chromosome[size];

            for (int i = leftBound; i <= rightBound; i++)
            {
                childOne[i] = new Chromosome(parentOne[i]);
                childTwo[i] = new Chromosome(parentTwo[i]);
                //permutationOne[i] = parentOne.permutation[i];
                //permutationTwo[i] = parentTwo.permutation[i];
            }

            int currentIndex = (rightBound  + 1) % size;
            int anotherParentIndex = currentIndex;
            int toAdd;
            bool skip;
            while (currentIndex != leftBound)
            {
                skip = false;
                toAdd = parentTwo.permutation[anotherParentIndex];
                for (int i = leftBound; i <= rightBound; i++)
                {
                    if (toAdd == parentOne.permutation[i])
                    {
                        anotherParentIndex = (anotherParentIndex + 1) % size;
                        skip = true;
                        break;
                    }
                }
                if (!skip)
                {
                    //permutationOne[currentIndex] = parentTwo.permutation[anotherParentIndex];
                    childOne[currentIndex] = new Chromosome(parentTwo[anotherParentIndex]);
                    anotherParentIndex = (anotherParentIndex + 1) % size;
                    currentIndex = (currentIndex + 1) % size;
                }
            }
                    
            currentIndex = (rightBound + 1) % size;
            anotherParentIndex = currentIndex;
            while (currentIndex != leftBound)
            {
                skip = false;
                toAdd = parentOne.permutation[anotherParentIndex];
                for (int i = leftBound; i <= rightBound; i++)
                {
                    if (toAdd == parentTwo.permutation[i])
                    {
                        anotherParentIndex = (anotherParentIndex + 1) % size;
                        skip = true;
                        break;
                    }
                }
                if (!skip)
                {
                    //permutationTwo[currentIndex] = parentOne.permutation[anotherParentIndex];
                    childTwo[currentIndex] = new Chromosome(parentOne[anotherParentIndex]);
                    anotherParentIndex = (anotherParentIndex + 1) % size;
                    currentIndex = (currentIndex + 1) % size;
                }
            }

            /*for (int i = 0; i < size; i++)
            {
                Console.Write(permutationOne[i] + " ");
            }

            Console.WriteLine();
            for (int i = 0; i < size; i++)
            {
                Console.Write(permutationTwo[i] + " ");
            }*/

            childrenArray[0] = new Solution(childOne);
            childrenArray[1] = new Solution(childTwo);
            return childrenArray;
        }
    }


    class PmxCrossoverOperator : ICrossoverOperator
    {
        public double CrossoverProbability { get; set; }

        private static Random rand = new Random();

        public PmxCrossoverOperator(double probability = 0.0)
        {
            this.CrossoverProbability = probability;
        }

        public ISolution[] Execute(IPopulation population)
        {
            Solution[] crossedPopulation = new Solution[population.Size];
            List<Solution> chosenSolutions = new List<Solution>();

            foreach (Solution sol in population)
            {
                double ifChosen = rand.NextDouble();
                if (ifChosen < this.CrossoverProbability)
                {
                    Solution tempSol = new Solution(sol);
                    chosenSolutions.Add(tempSol);
                }
            }

            if (chosenSolutions.Count != population.Size)
            {
                int stop = population.Size - chosenSolutions.Count;
                for (int i = 0; i < stop; i++)
                {
                    int whichSol = rand.Next(chosenSolutions.Count - 1);
                    int position = rand.Next(chosenSolutions.Count - 1);
                    Solution tempSol = new Solution(chosenSolutions[whichSol]);
                    chosenSolutions.Insert(position, tempSol);
                }
            }

            for (int i = 0; i < population.Size; i += 2)
            {
                Solution[] temp = new Solution[2];
                temp = Execute(chosenSolutions[i], chosenSolutions[i + 1]);
                crossedPopulation[i] = new Solution(temp[0]);
                crossedPopulation[i + 1] = new Solution(temp[1]);
            }
            return crossedPopulation;
        }

        public Solution[] Execute(Solution parentOne, Solution parentTwo)
        {
            int size = parentOne.Size;

            //int[] permutationOne = new int[size];
            //int[] permutationTwo = new int[size];

            int leftBound = rand.Next(size - 1);
            int rightBound = rand.Next(size - 1);
            Solution[] childrenArray = new Solution[2];
            List<Tuple<int, int>> MappingArray = new List<Tuple<int, int>>();

            if (leftBound > rightBound)
            {
                int temp = rightBound;
                rightBound = leftBound;
                leftBound = temp;
            }

            Chromosome[] childOne = new Chromosome[size];
            Chromosome[] childTwo = new Chromosome[size];

            for (int i = leftBound; i <= rightBound; i++)
            {
                childOne[i] = new Chromosome(parentTwo[i]);
                childTwo[i] = new Chromosome(parentOne[i]);
                //permutationOne[i] = parentTwo.permutation[i];
                //permutationTwo[i] = parentOne.permutation[i];
                Tuple<int, int> MappingTuple = new Tuple<int, int>(parentOne.permutation[i], parentTwo.permutation[i]);
                MappingArray.Add(MappingTuple);
            }

            bool bondingStopped = false;
            //sklejanie mapowan
            for (int i = 0; i < MappingArray.Count; i++)
            {
                if (MappingArray[i].Item1 == MappingArray[i].Item2)
                {
                    MappingArray.RemoveAt(i);
                    i--;
                }
            }
            while (!bondingStopped)
            {
                bondingStopped = true;
                for (int i = 0; i < MappingArray.Count; i++)
                {
                    for (int j = 0; j < MappingArray.Count; j++)
                    {
                        if (MappingArray[i].Item2 == MappingArray[j].Item1)
                        {
                            bondingStopped = false;
                            MappingArray[i] = new Tuple<int, int>(MappingArray[i].Item1, MappingArray[j].Item2);
                            MappingArray.RemoveAt(j);
                            break;
                        }
                    }
                }
            }

            bool mappingFound = false;

            for (int i = 0; i < size; i++)
            {
                if (i == leftBound)
                {
                    i = rightBound;
                    continue;
                }
                for (int j = 0; j < MappingArray.Count; j++)
                {
                    if (parentOne.permutation[i] == MappingArray[j].Item2)
                    {
                        for (int k = 0; k < size; k++)
                        {
                            if (parentTwo.permutation[k] == MappingArray[j].Item1)
                            {
                                childOne[i] = new Chromosome(parentTwo[k]);
                                //permutationOne[i] = parentTwo.permutation[k];
                                mappingFound = true;
                                break;
                            }
                        }
                    }
                    if (mappingFound) break;
                }
                if(!mappingFound)
                {
                    childOne[i] = new Chromosome(parentOne[i]);
                    //permutationOne[i] = parentOne.permutation[i];
                }
                mappingFound = false;
            }

            mappingFound = false;
            for (int i = 0; i < size; i++)
            {
                if (i == leftBound)
                {
                    i = rightBound;
                    continue;
                }
                for (int j = 0; j < MappingArray.Count; j++)
                {
                    if (parentTwo.permutation[i] == MappingArray[j].Item1)
                    {
                        for (int k = 0; k < size; k++)
                        {
                            if (parentOne.permutation[k] == MappingArray[j].Item2)
                            {
                                childTwo[i] = new Chromosome(parentOne[k]);
                                //permutationTwo[i] = parentOne.permutation[k];
                                mappingFound = true;
                                break;
                            }
                        }
                    }
                    if (mappingFound) break;
                }
                if (!mappingFound)
                {
                    childTwo[i] = new Chromosome(parentTwo[i]);
                    //permutationTwo[i] = parentTwo.permutation[i];
                }
                mappingFound = false;
            }

            /*Console.WriteLine("lef: " + leftBound + " right: " + rightBound);
            for (int i = 0; i < size; i++)
            {
                Console.Write(permutationOne[i] + "  ");
            }
            Console.WriteLine();

            for (int i = 0; i < size; i++) 
            {
                Console.Write(permutationTwo[i] + "  ");
            }
            Console.WriteLine();*/

            childrenArray[0] = new Solution(childOne);
            childrenArray[1] = new Solution(childTwo);
            return childrenArray;
        }
    }

    class CxCrossoverOperator : ICrossoverOperator
    {
        public double CrossoverProbability { get; set; }
        private static Random rand = new Random();

        public CxCrossoverOperator(double probability = 0.0)
        {
            this.CrossoverProbability = probability;
        }

        public ISolution[] Execute(IPopulation population)
        {
            Solution[] crossedPopulation = new Solution[population.Size];
            List<Solution> chosenSolutions = new List<Solution>();

            foreach (Solution sol in population)
            {
                double ifChosen = rand.NextDouble();
                if (ifChosen < this.CrossoverProbability)
                {
                    Solution tempSol = new Solution(sol);
                    chosenSolutions.Add(tempSol);
                }
            }

            if (chosenSolutions.Count != population.Size)
            {
                int stop = population.Size - chosenSolutions.Count;
                for (int i = 0; i < stop; i++)
                {
                    int whichSol = rand.Next(chosenSolutions.Count - 1);
                    int position = rand.Next(chosenSolutions.Count - 1);
                    Solution tempSol = new Solution(chosenSolutions[whichSol]);
                    chosenSolutions.Insert(position, tempSol);
                }
            }

            for (int i = 0; i < population.Size; i += 2)
            {
                Solution[] temp = new Solution[2];
                temp = Execute(chosenSolutions[i], chosenSolutions[i + 1]);
                crossedPopulation[i] = new Solution(temp[0]);
                crossedPopulation[i + 1] = new Solution(temp[1]);
            }
            return crossedPopulation;
        }

        public Solution[] Execute(Solution parentOne, Solution parentTwo)
        {
            Solution[] childrenArray = new Solution[2];
            int size = parentOne.Size;
            List<Tuple<Chromosome, int, Chromosome, int, int>> tupleList = new List<Tuple<Chromosome, int, Chromosome, int, int>>();
            Chromosome[] childOne = new Chromosome[size];
            Chromosome[] childTwo = new Chromosome[size];
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
                    }
                    else
                    {
                        childOne[tupleList[tupleIndex].Item5] = new Chromosome(tupleList[tupleIndex].Item3);
                        childTwo[tupleList[tupleIndex].Item5] = new Chromosome(tupleList[tupleIndex].Item1);
                    }
                    tupleList.RemoveAt(tupleIndex);
                }
                while (currentValue != stopCondition);
                swap = !swap;
            }
            childrenArray[0] = new Solution(childOne);
            childrenArray[1] = new Solution(childTwo);
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

    class SelectionOperator
    {
        private List<double> distribution = null;
        private static Random rand = new Random();
        public double CrossoverProbability { get; set; }
        private double bigNumber;

        public SelectionOperator(double probability, int solSize)
        {
            this.bigNumber = 0.0;
            double maxFlow = 0.0;
            double maxDistance = 0.0;
            this.CrossoverProbability = probability;
            for (int i = 0; i < solSize; i++)
            {
                for (int j = 0; j < solSize; j++)
                {
                    if (maxFlow < (QapData.Instance.getFlow())[i, j]) maxFlow = (QapData.Instance.getFlow())[i, j];
                    if (maxDistance < (QapData.Instance.getDistance())[i, j]) maxFlow = (QapData.Instance.getDistance())[i, j];
                }
            }
            bigNumber = (maxDistance + maxFlow) * solSize;
        }

        Population RouletteMethod(Population population)
        {
            int size = population.Size;
            Solution[] chosenSolutions = new Solution[size];
            this.distribution = new List<double>();

            double fitSum = 0.0;
            double currentFitSum = 0.0;

            foreach (Solution sol in population)
            {
                fitSum += (this.bigNumber - sol.Goal);
            }

            foreach (Solution sol in population)
            {
                currentFitSum += (this.bigNumber - sol.Goal);
                double probability = currentFitSum / fitSum;
                distribution.Add(probability);
            }

            this.distribution[size - 1] = 1.0;

            List<Solution> solutionsToMix = new List<Solution>();

            for (int i = 0; i < size; i++)
            {
                double toChoose = rand.NextDouble();
                for (int j = 0; j < size; j++)
                {
                    if (toChoose >= distribution[j])
                    {
                        chosenSolutions[i] = new Solution((Solution)population[i]);
                        break;
                    }
                }
            }

            return new Population(chosenSolutions);
        }

        Solution[] RankingMethod(Population population)
        {
            return new Solution[0];
        }
    }

    class QgAlgorithm : IEvolutionAlgorithm
    {
        PmxCrossoverOperator crossOperator = null;
        //zmienic na private
        public OxCrossoverOperator oxOperator = null;
        public CxCrossoverOperator cxOperator = null;
        public PmxCrossoverOperator pmxOperator = null;
        public MutationOperator mutOperator = null;
        private CatastropheOperator catOperator = null;
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

            oxOperator = new OxCrossoverOperator();

            pmxOperator = new PmxCrossoverOperator();

            mutOperator = new MutationOperator(mutProb, problemSize);

            rotOperator = new RotationGateOperator(problemSize);

            catOperator = new CatastropheOperator();
            catOperator.CastastropheProbability = catProb;

            this.isStopped = false;
            this.iterations = iterations;
            this.popSize = popSize;
            if (popSize % 2 != 0) popSize += 1;
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
