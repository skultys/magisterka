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
        private int observedState;
        public double Alpha { set; get; }
        public double Beta { set; get; }

        public int ObservedState
        {
            get
            {
                return this.observedState;
            }
        }

        public Qbit()
        {
            this.Alpha = 1.0 / Math.Sqrt(2.0);
            this.Beta = 1.0 / Math.Sqrt(2.0);
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

        public int  EvaluateState()
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

        public void ExecuteRotationGate(double theta)
        {
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
            PermutationValue = 0;
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
            this.PermutationValue = anotherChromosome.PermutationValue;
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

        public int PermutationValue { get; set; }
    }

    class Solution : ISolution
    {
        private List<Chromosome> chromosomes = null;
        //private List<int> permutation = null;
        int solSize;

        public int Size
        {
            get
            {
                return this.solSize;
            }
        }

        /*public List<int> Permutation
        {
            get
            {
                return this.permutation;
            }
        }*/

        public Solution(int size)
        {
            this.chromosomes = new List<Chromosome>();
            //this.permutation = new List<int>();
            int chromosomeSize = (int)(Math.Log(size, 2.0) + 1);
            this.Goal = 0.0;
            this.solSize = size;
            for (int i = 0; i < size; i++)
            {
                Chromosome chromosome = new Chromosome(chromosomeSize);
                this.chromosomes.Add(chromosome);
                //this.permutation.Add(0);
            }
            toPermutation();
        }

        public Solution(Solution anotherSolution)
        {
            this.chromosomes = new List<Chromosome>();
            //this.permutation = new List<int>();
            for (int i = 0; i < anotherSolution.chromosomes.Count; i++)
            {
                Chromosome chr = new Chromosome(anotherSolution.chromosomes[i]);
                this.chromosomes.Add(chr);
            }

            /*for (int i = 0; i < anotherSolution.permutation.Count; i++)
            {
                this.permutation.Add(anotherSolution.permutation[i]);
            }*/

            this.Goal = anotherSolution.Goal;
            this.solSize = anotherSolution.solSize;
        }

        public Solution(Chromosome[] chromosomes)
        {
            this.solSize = chromosomes.Length;
            this.Goal = 0.0;
            this.chromosomes = new List<Chromosome>();
            foreach (Chromosome chr in chromosomes)
            {
                Chromosome temp = new Chromosome(chr);
                this.chromosomes.Add(temp);
            }
            /*this.permutation = new List<int>();
            for (int i = 0; i < this.Size; i++)
            {
                this.permutation.Add(chromosomes[i].PermutationValue);
            }*/

            EvaluateGoal();
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
            List<int> permutation = new List<int>();

            for (int i = 0; i < this.chromosomes.Count; i++)
            {
                Tuple<int, int> tuple = new Tuple<int, int>(i, this.chromosomes[i].DecodeChromosome());
                tupleList.Add(tuple);
                permutation.Add(0);
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
                permutation[tuple.Item1] = val;
            }

            for (int i = 0; i < this.chromosomes.Count; i++)
            {
                this.chromosomes[i].PermutationValue = permutation[i];
            }

            EvaluateGoal();
        }

        void EvaluateGoal()
        {
            this.Goal = 0.0;
            double[,] flowM = QapData.Instance.getFlow();
            double[,] distanceM = QapData.Instance.getDistance();


            /*this[0].PermutationValue = 3;
            this[1].PermutationValue = 10;
            this[2].PermutationValue = 11;
            this[3].PermutationValue = 2;
            this[4].PermutationValue = 12;
            this[5].PermutationValue = 5;
            this[6].PermutationValue = 6;
            this[7].PermutationValue = 7;
            this[8].PermutationValue = 8;
            this[9].PermutationValue = 1;
            this[10].PermutationValue = 4;
            this[11].PermutationValue = 9;*/

            for (int i = 0; i < this.chromosomes.Count; i++)
            {
                for (int j = 0; j < this.chromosomes.Count; j++)
                {
                    this.Goal += distanceM[i, j] * flowM[this.chromosomes[i].PermutationValue - 1, this.chromosomes[j].PermutationValue - 1];
                }
            }
        }

        public double Goal { get; set; }

        public void PrintSolution()
        {
            foreach (Chromosome chr in this.chromosomes)
            {
                Console.Write(chr.PermutationValue + "  ");
            }
            Console.WriteLine();
        }
    }

    class Population : IPopulation
    {
        private int popSize;
        private int problemSize;
        private int currPopSize;
        private List<Solution> solutions = null;

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

        public Population(Population anotherPopulation)
        {
            this.solutions = new List<Solution>();

            foreach (Solution sol in anotherPopulation)
            {
                Solution tempSol = new Solution(sol);
                this.solutions.Add(tempSol);
                this.popSize = anotherPopulation.Size;
                this.problemSize = anotherPopulation.problemSize;
                this.currPopSize = anotherPopulation.currPopSize;
            }
        }

        public Population(ISolution[] solutions)
        {
            this.solutions = new List<Solution>();
            foreach (Solution sol in solutions)
            {
                Solution temp = new Solution(sol);
                this.solutions.Add(temp);
            }
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
            
            if (chosenSolutions.Count == 0)
            {
                for (int i = 0; i < population.Size; i++)
                {
                    int index = rand.Next(population.Size - 1);
                    Solution tempSol = new Solution((Solution)population[index]);
                    chosenSolutions.Add(tempSol);
                }
            }
            else if (chosenSolutions.Count != population.Size)
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
                //toAdd = parentTwo.Permutation[anotherParentIndex];
                toAdd = parentTwo[anotherParentIndex].PermutationValue;
                for (int i = leftBound; i <= rightBound; i++)
                {
                    //if (toAdd == parentOne.Permutation[i])
                    if (toAdd == parentOne[i].PermutationValue)
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
                //toAdd = parentOne.Permutation[anotherParentIndex];
                toAdd = parentOne[anotherParentIndex].PermutationValue;
                for (int i = leftBound; i <= rightBound; i++)
                {
                    //if (toAdd == parentTwo.Permutation[i])
                    if (toAdd == parentTwo[i].PermutationValue)
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
            if (chosenSolutions.Count == 0)
            {
                for (int i = 0; i < population.Size; i++)
                {
                    int index = rand.Next(population.Size - 1);
                    Solution tempSol = new Solution((Solution)population[index]);
                    chosenSolutions.Add(tempSol);
                }
            }
            else if (chosenSolutions.Count != population.Size)
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
                //Tuple<int, int> MappingTuple = new Tuple<int, int>(parentOne.Permutation[i], parentTwo.Permutation[i]);
                Tuple<int, int> MappingTuple = new Tuple<int, int>(parentOne[i].PermutationValue, parentTwo[i].PermutationValue);
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
                    //if (parentOne.Permutation[i] == MappingArray[j].Item2)
                    if (parentOne[i].PermutationValue == MappingArray[j].Item2)
                    {
                        for (int k = 0; k < size; k++)
                        {
                            //if (parentTwo.Permutation[k] == MappingArray[j].Item1)
                            if (parentTwo[k].PermutationValue == MappingArray[j].Item1)
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
                    //if (parentTwo.Permutation[i] == MappingArray[j].Item1)
                    if (parentTwo[i].PermutationValue == MappingArray[j].Item1)
                    {
                        for (int k = 0; k < size; k++)
                        {
                            //if (parentOne.Permutation[k] == MappingArray[j].Item2)
                            if (parentOne[k].PermutationValue == MappingArray[j].Item2)
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

            if (chosenSolutions.Count == 0)
            {
                for (int i = 0; i < population.Size; i++)
                {
                    int index = rand.Next(population.Size - 1);
                    Solution tempSol = new Solution((Solution)population[index]);
                    chosenSolutions.Add(tempSol);
                }
            }
            else if (chosenSolutions.Count != population.Size)
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
                //Tuple<Chromosome, int, Chromosome, int, int> tuple = new Tuple<Chromosome, int, Chromosome, int, int>(parentOne[i], parentOne.Permutation[i], parentTwo[i], parentTwo.Permutation[i], i);
                Tuple<Chromosome, int, Chromosome, int, int> tuple = new Tuple<Chromosome, int, Chromosome, int, int>(parentOne[i], parentOne[i].PermutationValue, parentTwo[i], parentTwo[i].PermutationValue, i);
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
        private int solSize;

        private int bitsInSol;

        public RotationGateOperator(int solSize)
        {
            this.solSize = solSize;
            this.bitsInSol = (int)(Math.Log(solSize, 2.0) + 1);
        }

        public void Execute(IPopulation population, Solution best)
        {
            foreach (Solution sol in population)
            {
                if()
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

        public Population RouletteMethod(IPopulation population)
        {
            population.SortAscending();
            int size = population.Size;
            double bigNumber = population[population.Size - 1].Goal;
            Solution[] chosenSolutions = new Solution[size];
            this.distribution = new List<double>();

            double fitSum = 0.0;
            double currentFitSum = 0.0;

            foreach (Solution sol in population)
            {
                fitSum += (bigNumber - sol.Goal);
            }

            foreach (Solution sol in population)
            {
                currentFitSum += (bigNumber - sol.Goal);
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
                    if (toChoose <= distribution[j])
                    {
                        chosenSolutions[i] = new Solution((Solution)population[j]);
                        break;
                    }
                }
            }

            return new Population(chosenSolutions);
        }

        public Solution[] RankingMethod(Population population)
        {
            return new Solution[0];
        }
    }

    class QgAlgorithm : IEvolutionAlgorithm
    {
        private OxCrossoverOperator oxOperator = null;
        private CxCrossoverOperator cxOperator = null;
        private PmxCrossoverOperator pmxOperator = null;
        private MutationOperator mutOperator = null;
        private CatastropheOperator catOperator = null;
        private RotationGateOperator rotOperator = null;
        private SelectionOperator selOPerator = null;
        bool isStopped;
        int iterations;
        int popSize;
        int problemSize;
        public Solution best = null;

        public QgAlgorithm(double[,] distance, double[,] flow, double crossProb, double mutProb, double catProb, int iterations, int popSize, int problemSize)
        {
            QapData.Instance.setQapData(distance, flow);

            cxOperator = new CxCrossoverOperator(crossProb);
            oxOperator = new OxCrossoverOperator(crossProb);
            pmxOperator = new PmxCrossoverOperator(crossProb);
            mutOperator = new MutationOperator(mutProb, problemSize);
            rotOperator = new RotationGateOperator(problemSize);
            catOperator = new CatastropheOperator();
            catOperator.CastastropheProbability = catProb;
            selOPerator = new SelectionOperator();
            this.isStopped = false;
            this.iterations = iterations;
            this.popSize = popSize;
            if (popSize % 2 != 0) popSize += 1;
            this.problemSize = problemSize;
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
            InitRandomPopulation();
            //double avarage;
            for (int i = 1; i < this.iterations; i++)
            {
                GetBestSolution();
                //avarage = 0.0;
                //foreach (Solution sol in this.Population) avarage += sol.Goal;
                //Console.WriteLine("Srednia wartosc startowa w iteracji nr        " + i + " : " + avarage / this.popSize);

                //this.Population = new Population(this.selOPerator.RouletteMethod(this.Population));
                //avarage = 0.0;
                //foreach (Solution sol in this.Population) avarage += sol.Goal;
                //Console.WriteLine("Srednia wartosc po selekcji w iteracji nr     " + i + " : " + avarage / this.popSize);

                //this.Population = new Population(this.pmxOperator.Execute(this.Population));
                //avarage = 0.0;
                //foreach (Solution sol in this.Population) avarage += sol.Goal;
                //Console.WriteLine("Srednia wartosc po krzyzowaniu w iteracji nr  " + i + " : " + avarage / this.popSize);

                //this.mutOperator.Execute(this.Population);

                this.rotOperator.Execute(this.Population, this.best);
                /*avarage = 0.0;
                foreach (Solution sol in this.Population) avarage += sol.Goal;
                Console.WriteLine("Srednia wartosc po rotacji w iteracji nr      " + i + " : " + avarage / this.popSize);*/
                //Console.WriteLine();
                //Console.ReadKey();
            }

            this.isStopped = true;
            GetBestSolution();
            this.best.PrintSolution();
            Console.WriteLine(this.best.Goal);
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

    class testOne
    {
        public void execute(testTwo two)
        {
            two.Dupa = 1;
        }
    }

    class testTwo
    {
        public int Dupa { set; get; }
        public testTwo()
        {
            this.Dupa = 666;
        }
    }
}
