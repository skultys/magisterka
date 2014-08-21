using System;
using System.Collections.Generic;
using System.Text;
using OptimisationClassLibrary;
using System.Collections;
using System.IO;

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

        public bool ExecuteRotationGate(Qbit best)
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
            this.Beta = Math.Sin(theta) * tempAlpha + Math.Cos(theta) * this.Beta;
            if (theta != 0.0) return true;
            else return false;
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
            int chromosomeSize = (int)(Math.Log(size, 2.0) + 1);
            this.Goal = 0.0;
            this.solSize = size;
            for (int i = 0; i < size; i++)
            {
                Chromosome chromosome = new Chromosome(chromosomeSize);
                this.chromosomes.Add(chromosome);
            }
            toPermutation();
        }

        public Solution(Solution anotherSolution)
        {
            this.chromosomes = new List<Chromosome>();
            for (int i = 0; i < anotherSolution.chromosomes.Count; i++)
            {
                Chromosome chr = new Chromosome(anotherSolution.chromosomes[i]);
                this.chromosomes.Add(chr);
            }

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

        public void PrintSolution(System.IO.StreamWriter file)
        {
            foreach (Chromosome chr in this.chromosomes)
            {
                file.Write(chr.PermutationValue + "  ");
            }
            file.WriteLine();
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

        public Population(IPopulation anotherPopulation)
        {
            this.solutions = new List<Solution>();

            foreach (Solution sol in anotherPopulation)
            {
                Solution tempSol = new Solution(sol);
                this.solutions.Add(tempSol);
                this.popSize = anotherPopulation.Size;
                this.problemSize = anotherPopulation.Size;
                this.currPopSize = anotherPopulation.CurrentSize;
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

            int leftBound = rand.Next(size - 1);
            int rightBound = rand.Next(size - 1);
            Solution[] childrenArray = new Solution[2];
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
                childOne[i] = new Chromosome(parentOne[i]);
                childTwo[i] = new Chromosome(parentTwo[i]);
            }

            int currentIndex = (rightBound  + 1) % size;
            int anotherParentIndex = currentIndex;
            int toAdd;
            bool skip;
            while (currentIndex != leftBound)
            {
                skip = false;
                toAdd = parentTwo[anotherParentIndex].PermutationValue;
                for (int i = leftBound; i <= rightBound; i++)
                {
                    if (toAdd == parentOne[i].PermutationValue)
                    {
                        anotherParentIndex = (anotherParentIndex + 1) % size;
                        skip = true;
                        break;
                    }
                }
                if (!skip)
                {
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
                toAdd = parentOne[anotherParentIndex].PermutationValue;
                for (int i = leftBound; i <= rightBound; i++)
                {
                    if (toAdd == parentTwo[i].PermutationValue)
                    {
                        anotherParentIndex = (anotherParentIndex + 1) % size;
                        skip = true;
                        break;
                    }
                }
                if (!skip)
                {
                    childTwo[currentIndex] = new Chromosome(parentOne[anotherParentIndex]);
                    anotherParentIndex = (anotherParentIndex + 1) % size;
                    currentIndex = (currentIndex + 1) % size;
                }
            }

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
                Tuple<int, int> MappingTuple = new Tuple<int, int>(parentOne[i].PermutationValue, parentTwo[i].PermutationValue);
                MappingArray.Add(MappingTuple);
            }

            bool bondingStopped = false;
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
                    if (parentOne[i].PermutationValue == MappingArray[j].Item2)
                    {
                        for (int k = 0; k < size; k++)
                        {
                            if (parentTwo[k].PermutationValue == MappingArray[j].Item1)
                            {
                                childOne[i] = new Chromosome(parentTwo[k]);
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
                    if (parentTwo[i].PermutationValue == MappingArray[j].Item1)
                    {
                        for (int k = 0; k < size; k++)
                        {
                            if (parentOne[k].PermutationValue == MappingArray[j].Item2)
                            {
                                childTwo[i] = new Chromosome(parentOne[k]);
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
                }
                mappingFound = false;
            }

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
            bool changed = false;
            foreach (Solution sol in population)
            {
                if (best.Goal < sol.Goal)
                {
                    for (int i = 0; i < this.solSize; i++)
                    {
                        for (int j = 0; j < this.bitsInSol; j++)
                        {
                            changed = changed || sol[i][j].ExecuteRotationGate(best[i][j]);
                        }
                    }
                    if (changed)
                    {
                        sol.toPermutation();
                        changed = false;
                    }
                    else
                    {
                        int b;
                    }
                }
                else
                {
                    int a;
                }
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
        private double EtaMin;
        public double EtaMax { get; set; }
        private double biggestGoal;

        public SelectionOperator(double etaMax = 2.0)
        {
            this.EtaMax = etaMax;
            if (this.EtaMax > 2.0 || this.EtaMax < 1.0) this.EtaMax = 2.0;
            this.EtaMin = 2.0 - this.EtaMax;
            this.biggestGoal = -100000.0;
        }

        public Population RouletteMethod(IPopulation population)
        {
            population.SortAscending();
            int size = population.Size;
            double bigNumber = population[population.Size - 1].Goal;
            if (this.biggestGoal == -100000.0 || this.biggestGoal < bigNumber) this.biggestGoal = bigNumber;
            Solution[] chosenSolutions = new Solution[size];
            this.distribution = new List<double>();

            double fitSum = 0.0;
            double currentFitSum = 0.0;

            foreach (Solution sol in population)
            {
                fitSum += (this.biggestGoal - sol.Goal);
            }

            foreach (Solution sol in population)
            {
                currentFitSum += (this.biggestGoal - sol.Goal);
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

        public Population RankingMethod(IPopulation population)
        {
            population.SortAscending();
            this.distribution = new List<double>();
            Solution[] chosenSolutions = new Solution[population.Size];
            double[] probabilities = new double[population.Size];
            double probabilitySum = 0.0;
            double currentSum = 0.0;

            for (int i = 0; i < population.Size; i++)
            {
                double probability = (1.0 / population.Size) * (this.EtaMax - ((this.EtaMax - this.EtaMin) * (Convert.ToDouble(i) / (population.Size - 1.0))));
                probabilities[i] = probability;
            }

            foreach (double prob in probabilities) probabilitySum += prob;

            for (int i = 0; i < population.Size; i++)
            {
                currentSum += probabilities[i];
                double toAdd = currentSum / probabilitySum;
                this.distribution.Add(toAdd);
            }
            this.distribution[population.Size - 1] = 1.0;

            for (int i = 0; i < population.Size; i++)
            {
                double toChoose = rand.NextDouble();
                for (int j = 0; j < population.Size; j++)
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
    }

    class QgAlgorithm : IEvolutionAlgorithm
    {
        private CrossoverOperatorChosen crossOperChosen;
        private SelectionMethodChosen selMethodChosen;
        private OxCrossoverOperator oxOperator = null;
        private CxCrossoverOperator cxOperator = null;
        private PmxCrossoverOperator pmxOperator = null;
        private MutationOperator mutOperator = null;
        private RotationGateOperator rotOperator = null;
        private SelectionOperator selOPerator = null;
        bool isStopped;
        int iterations;
        int popSize;
        int problemSize;
        double crossoverProbMax;
        double crossoverProbMin;
        double mutationProb;
        public Solution best = null;
        private bool isSet;
        private bool saveEnabled;
        private bool rotGateEnabled;
        private string fileName;
        private string instance;
        private Population initialPopulationCopy;

        public QgAlgorithm()
        {
            this.isStopped = false;
            this.isSet = false;
            this.saveEnabled = false;
            string path = Directory.GetCurrentDirectory() + "\\Results";
            if (!Directory.Exists(path)) Directory.CreateDirectory(path);

            path = Directory.GetCurrentDirectory() + "\\QAPLib";
            if (!Directory.Exists(path)) Directory.CreateDirectory(path);
        }

        public bool SetParameters(bool newPopulation, bool onlyFileName)
        {
            string enteredValue;
            double doubleValue;
            int intValue;
            bool OK = false;
            ConsoleKeyInfo cki;

            if (!onlyFileName)
            {
                if (newPopulation)
                {
                    while (!OK)
                    {
                        try
                        {
                            Console.Clear();
                            Console.WriteLine("Wybierz instancje testowa:");
                            Console.WriteLine();
                            string path = Directory.GetCurrentDirectory() + "\\QAPLib";
                            string[] filePaths = Directory.GetFiles(path, "*.dat");
                            if (filePaths.Length == 0)
                            {
                                Console.Clear();
                                Console.WriteLine("Folder jest pusty!");
                                Console.ReadKey();
                                return false;
                            }
                            for (int i = 0; i < filePaths.Length; i++)
                            {
                                filePaths[i] = Path.GetFileName(filePaths[i]);
                                Console.WriteLine((i + 1) + ". " + filePaths[i]);
                            }
                            enteredValue = Console.ReadLine();
                            intValue = Convert.ToInt32(enteredValue);
                            if (intValue < 1 || intValue > filePaths.Length)
                            {
                                Console.Clear();
                                Console.WriteLine("Zla wartosc!");
                                Console.ReadKey();
                            }
                            else
                            {
                                try
                                {
                                    intValue--;
                                    path = path + "\\" + filePaths[intValue];
                                    string[] testData = System.IO.File.ReadAllLines(@path);

                                    int length = Convert.ToInt32(testData[0]);

                                    double[,] distance = new double[length, length];

                                    double[,] flow = new double[length, length];

                                    List<string> lines = new List<string>();

                                    foreach (string s in testData)
                                    {
                                        string s2 = s.Trim();
                                        if (s2 != string.Empty) lines.Add(s2);
                                    }

                                    for (int i = 1; i <= length; i++)
                                    {
                                        int index2;
                                        int j = 0;
                                        while (j < length)
                                        {
                                            index2 = lines[i].IndexOf(" ", 0);
                                            if (index2 < 0) index2 = lines[i].Length;
                                            string number = lines[i].Substring(0, index2);
                                            distance[i - 1, j] = Convert.ToDouble(number);
                                            lines[i] = (lines[i].Substring(index2)).Trim();

                                            index2 = lines[lines.Count - i].IndexOf(" ", 0);
                                            if (index2 < 0) index2 = lines[lines.Count - i].Length;
                                            number = lines[lines.Count - i].Substring(0, index2);
                                            flow[length - i, j] = Convert.ToDouble(number);
                                            lines[lines.Count - i] = (lines[lines.Count - i].Substring(index2)).Trim();

                                            j++;
                                        }
                                    }
                                    QapData.Instance.setQapData(distance, flow);
                                    this.problemSize = length;
                                    this.rotOperator = new RotationGateOperator(this.problemSize);
                                    this.instance = filePaths[intValue];
                                    OK = true;
                                }
                                catch (Exception e)
                                {
                                    Console.Clear();
                                    Console.WriteLine("Nie mozna otworzyc pliku!");
                                    Console.ReadKey();
                                    return false;
                                }
                            }
                        }
                        catch (FormatException e)
                        {
                            Console.Clear();
                            Console.WriteLine("Zla wartosc!");
                            Console.ReadKey();
                        }
                        catch (Exception e)
                        {
                            Console.Clear();
                            Console.WriteLine("Nie mozna otworzyc pliku!");
                            Console.ReadKey();
                            return false;
                        }
                    }


                    OK = false;
                    while (!OK)
                    {
                        try
                        {
                            Console.Clear();
                            Console.WriteLine("Podaj rozmiar populacji: ");
                            enteredValue = Console.ReadLine();
                            intValue = Convert.ToInt32(enteredValue);
                            if (intValue > 0)
                            {
                                this.popSize = intValue;
                                if (popSize % 2 != 0) popSize += 1;
                                OK = true;
                            }
                            else
                            {
                                Console.Clear();
                                Console.WriteLine("Zla wartosc.");
                                Console.ReadKey();
                            }
                        }
                        catch (Exception e)
                        {
                            Console.Clear();
                            Console.WriteLine("Zla wartosc.");
                            Console.ReadKey();
                        }
                    }
                    OK = false;
                }

                while (!OK)
                {
                    try
                    {
                        Console.Clear();
                        Console.WriteLine("Podaj ilosc iteracji: ");
                        enteredValue = Console.ReadLine();
                        intValue = Convert.ToInt32(enteredValue);
                        if (intValue > 0)
                        {
                            this.iterations = intValue;
                            OK = true;
                        }
                        else
                        {
                            Console.Clear();
                            Console.WriteLine("Zla wartosc.");
                            Console.ReadKey();
                        }
                    }
                    catch (Exception e)
                    {
                        Console.Clear();
                        Console.WriteLine("Zla wartosc.");
                        Console.ReadKey();
                    }
                }

                OK = false;
                while (!OK)
                {
                    Console.Clear();
                    Console.WriteLine("Wybierz metode selekcji: ");
                    Console.WriteLine("1. Ruletkowa");
                    Console.WriteLine("2. Rankingowa");
                    cki = Console.ReadKey();
                    if (cki.Key == ConsoleKey.D1)
                    {
                        this.selMethodChosen = SelectionMethodChosen.Roulette;
                        this.selOPerator = new SelectionOperator();
                        OK = true;
                    }
                    else if (cki.Key == ConsoleKey.D2)
                    {
                        this.selMethodChosen = SelectionMethodChosen.Ranking;
                        do
                        {
                            try
                            {
                                Console.Clear();
                                Console.WriteLine("Podaj parametr eta (1,0 - 2,0): ");
                                enteredValue = Console.ReadLine();
                                doubleValue = Convert.ToDouble(enteredValue);
                                if (doubleValue >= 1.0 && doubleValue <= 2.0)
                                {
                                    this.selOPerator = new SelectionOperator(doubleValue);
                                    OK = true;
                                }
                                else
                                {
                                    Console.Clear();
                                    Console.WriteLine("Zla wartosc.");
                                    Console.ReadKey();
                                }
                            }
                            catch (Exception e)
                            {
                                Console.Clear();
                                Console.WriteLine("Zla wartosc.");
                                Console.ReadKey();
                            }
                        }
                        while (!OK);
                    }
                }

                OK = false;
                while (!OK)
                {
                    try
                    {
                        Console.Clear();
                        Console.WriteLine("Podaj maksymalne prawdopodobienstwo krzyzowania (0,0 - 1,0): ");
                        enteredValue = Console.ReadLine();
                        doubleValue = Convert.ToDouble(enteredValue);
                        if (doubleValue >= 0.0 && doubleValue <= 1.0)
                        {
                            this.crossoverProbMax = doubleValue;
                            OK = true;
                        }
                        else
                        {
                            Console.Clear();
                            Console.WriteLine("Zla wartosc.");
                            Console.ReadKey();
                        }
                    }
                    catch (Exception e)
                    {
                        Console.Clear();
                        Console.WriteLine("Zla wartosc.");
                        Console.ReadKey();
                    }
                }

                OK = false;
                while (!OK)
                {
                    try
                    {
                        Console.Clear();
                        Console.WriteLine("Podaj minimalne prawdopodobienstwo krzyzowania (0,0 - " + this.crossoverProbMax + "): ");
                        enteredValue = Console.ReadLine();
                        doubleValue = Convert.ToDouble(enteredValue);
                        if (doubleValue >= 0.0 && doubleValue <= 1.0 && doubleValue <= this.crossoverProbMax)
                        {
                            this.crossoverProbMin = doubleValue;
                            OK = true;
                        }
                        else
                        {
                            Console.Clear();
                            Console.WriteLine("Zla wartosc.");
                            Console.ReadKey();
                        }
                    }
                    catch (Exception e)
                    {
                        Console.Clear();
                        Console.WriteLine("Zla wartosc.");
                        Console.ReadKey();
                    }
                }

                OK = false;
                while (!OK)
                {
                    Console.Clear();
                    Console.WriteLine("Wybierz operator krzyzowania: ");
                    Console.WriteLine("1. CX");
                    Console.WriteLine("2. OX");
                    Console.WriteLine("3. PMX");
                    cki = Console.ReadKey();
                    if (cki.Key == ConsoleKey.D1)
                    {
                        this.crossOperChosen = CrossoverOperatorChosen.CX;
                        this.cxOperator = new CxCrossoverOperator(this.crossoverProbMax);
                        OK = true;
                    }
                    else if (cki.Key == ConsoleKey.D2)
                    {
                        this.crossOperChosen = CrossoverOperatorChosen.OX;
                        this.oxOperator = new OxCrossoverOperator(this.crossoverProbMax);
                        OK = true;
                    }
                    else if (cki.Key == ConsoleKey.D3)
                    {
                        this.crossOperChosen = CrossoverOperatorChosen.PMX;
                        this.pmxOperator = new PmxCrossoverOperator(this.crossoverProbMax);
                        OK = true;
                    }
                    else
                    {
                        Console.Clear();
                        Console.WriteLine("Zla wartosc.");
                        Console.ReadKey();
                    }
                }

                OK = false;
                while (!OK)
                {
                    try
                    {
                        Console.Clear();
                        Console.WriteLine("Podaj prawdopodobienstwo mutacji (0.0 - 1.0): ");
                        enteredValue = Console.ReadLine();
                        doubleValue = Convert.ToDouble(enteredValue);
                        if (doubleValue >= 0.0 && doubleValue <= 1.0)
                        {
                            this.mutOperator = new MutationOperator(doubleValue, this.problemSize);
                            this.mutationProb = doubleValue;
                            OK = true;
                        }
                        else
                        {
                            Console.Clear();
                            Console.WriteLine("Zla wartosc.");
                            Console.ReadKey();
                        }
                    }
                    catch (Exception e)
                    {
                        Console.Clear();
                        Console.WriteLine("Zla wartosc.");
                        Console.ReadKey();
                    }
                }

                cki = new ConsoleKeyInfo('1', ConsoleKey.D3, false, false, false);
                while (cki.Key != ConsoleKey.D1 && cki.Key != ConsoleKey.D2)
                {
                    Console.Clear();
                    Console.WriteLine("Czy chcesz uzyc operator bramki rotacyjnej?");
                    Console.WriteLine("1. Tak.");
                    Console.WriteLine("2. Nie.");
                    cki = Console.ReadKey();
                    if (cki.Key == ConsoleKey.D1)
                    {
                        this.rotGateEnabled = true;
                    }
                    else if (cki.Key == ConsoleKey.D2)
                    {
                        this.rotGateEnabled = false;
                    }
                }

            }
            OK = false;
            while (!OK)
            {
                Console.Clear();
                Console.WriteLine("Chcesz zapisac rezultaty dzialania algorytmu w pliku?");
                Console.WriteLine("1. Tak");
                Console.WriteLine("2. Nie");
                cki = Console.ReadKey();
                if (cki.Key == ConsoleKey.D1)
                {
                    bool wrongName = true;
                    try
                    {
                        while (wrongName)
                        {
                            this.fileName = "";
                            Console.Clear();
                            Console.WriteLine("Podaj nazwe pliku:");
                            this.fileName = Console.ReadLine();
                            if (this.fileName != "")
                            {
                                this.fileName = Directory.GetCurrentDirectory() + "\\Results\\" + this.fileName + ".txt";
                                this.saveEnabled = true;
                                wrongName = false;
                                OK = true;
                            }
                        }
                    }
                    catch (Exception e)
                    {
                        Console.Clear();
                        Console.WriteLine("Zla nazwa pliku!");
                        Console.ReadKey();
                    }
                }
                else if (cki.Key == ConsoleKey.D2)
                {
                    this.saveEnabled = false;
                    OK = true;
                }
            }


            Console.Clear();
            Console.WriteLine("Parametry pomyslnie ustawione!");
            Console.ReadKey();

            return true;
        }

        public IPopulation Population { get; set; }

        public double EvaluateSolution(ISolution solution)
        {
            return solution.Goal;
        }

        public void InitRandomPopulation()
        {
            if (isSet) this.Population = new Population(this.popSize, this.problemSize);
            else
            {
                Console.WriteLine("Parametry nie sa ustawione!");
                Console.ReadKey();
                Console.Clear();
            }
        }

        public void Execute()
        {
            bool exit = false;
            while (!exit)
            {
                if (isSet)
                {
                    Console.Clear();
                    Console.WriteLine("Algorytm rozpoczal dzialanie. Trwa poszukiwanie rozwiazania...");
                    GetBestSolution();
                    Solution bestSol = new Solution(this.best);
                    int bestIteration = 0;
                    double bestCurrentGoal = GetBestSolution().Goal;
                    Solution globalBest = new Solution(this.best);
                    double avarage = 0.0;
                    double previousAvarage = 0.0;
                    foreach (Solution sol in this.Population) avarage += sol.Goal;
                    avarage /= this.popSize;
                    int cntr = 0;

                    if (saveEnabled)
                    {
                        System.IO.StreamWriter file = new System.IO.StreamWriter(this.fileName);
                        file.WriteLine("Plik utworzono\t\t\t\t\t" + DateTime.Now);
                        file.WriteLine("Wybrana instancja testowa:\t\t\t" + this.instance);
                        file.WriteLine("Rozmiar problemu:\t\t\t\t" + this.problemSize);
                        file.WriteLine("Ilosc iteracji:\t\t\t\t\t" + this.iterations);
                        if (this.selMethodChosen == SelectionMethodChosen.Roulette)
                        {
                            file.WriteLine("Metoda selekcji:\t\t\t\tRuletkowa");
                        }
                        else
                        {
                            file.WriteLine("Metoda selekcji:\t\t\t\tRankingowa");
                            file.WriteLine("Parametr eta MAX:\t\t\t\t" + this.selOPerator.EtaMax);
                        }

                        if (this.crossOperChosen == CrossoverOperatorChosen.CX)
                        {
                            file.WriteLine("Operator krzyzowania:\t\t\t\tCX");
                        }
                        else if (this.crossOperChosen == CrossoverOperatorChosen.OX)
                        {
                            file.WriteLine("Operator krzyzowania:\t\t\t\tOX");
                        }
                        else
                        {
                            file.WriteLine("Operator krzyzowania:\t\t\t\tPMX");
                        }
                        file.WriteLine("Maksymalne prawdopodobienstwo krzyzowania:\t" + this.crossoverProbMax);
                        file.WriteLine("Minimalne prawdopodobienstwo krzyzowania:\t" + this.crossoverProbMax);
                        file.WriteLine("Prawdopodobienstwo mutacji:\t\t\t" + this.mutationProb);
                        file.WriteLine();
                        file.WriteLine("Nr iteracji | Najlepsza wartosc funkcji celu | Srednia wartosc funkcji celu");
                        file.WriteLine();
                        file.WriteLine("0\t" + this.best.Goal + "\t" + avarage);

                        for (int i = 1; i < this.iterations; i++)
                        {
                            double crossProb = this.crossoverProbMax / (1.0 + Convert.ToDouble(i)/Convert.ToDouble(this.iterations));
                            if (crossProb < this.crossoverProbMin) crossProb = this.crossoverProbMin;
                            previousAvarage = avarage;
                            if (this.selMethodChosen == SelectionMethodChosen.Roulette)
                            {
                                this.Population = new Population(this.selOPerator.RouletteMethod(this.Population));
                            }
                            else
                            {
                                this.Population = new Population(this.selOPerator.RankingMethod(this.Population));
                            }

                            if (this.crossOperChosen == CrossoverOperatorChosen.CX)
                            {
                                this.cxOperator.CrossoverProbability = crossProb;
                                this.Population = new Population(this.cxOperator.Execute(this.Population));
                            }
                            else if (this.crossOperChosen == CrossoverOperatorChosen.OX)
                            {
                                this.oxOperator.CrossoverProbability = crossProb;
                                this.Population = new Population(this.oxOperator.Execute(this.Population));
                            }
                            else
                            {
                                this.pmxOperator.CrossoverProbability = crossProb;
                                this.Population = new Population(this.pmxOperator.Execute(this.Population));
                            }

                            this.mutOperator.Execute(this.Population);

                            if (this.rotGateEnabled) this.rotOperator.Execute(this.Population, this.best);

                            avarage = 0.0;
                            foreach (Solution sol in this.Population) avarage += sol.Goal;
                            avarage /= this.popSize;

                            if (avarage > previousAvarage) cntr++;

                            GetBestSolution();

                            if (bestCurrentGoal > this.best.Goal)
                            {
                                bestCurrentGoal = this.best.Goal;
                                bestIteration = i;
                                globalBest = new Solution(this.best);
                            }
                            file.WriteLine(i + "\t" + this.best.Goal + " \t" + avarage);
                        }
                        file.WriteLine();
                        file.WriteLine("Znalezione rozwiazanie:");
                        this.best.PrintSolution(file);

                        file.WriteLine();
                        file.WriteLine("Wartosc funkcji celu:");
                        file.WriteLine(this.best.Goal);

                        file.WriteLine();
                        file.WriteLine("Najlepsze znalezione rozwiazanie:");
                        globalBest.PrintSolution(file);

                        file.WriteLine();
                        file.WriteLine("Wartosc funkcji celu najlepszego rozwiazania:");
                        file.WriteLine(globalBest.Goal);

                        file.WriteLine("Znalezione w " + bestIteration + " iteracji.");

                        file.WriteLine();
                        if (this.rotGateEnabled)
                        {
                            file.WriteLine("bramka rotacyjna: WYKORZYSTANA");
                        }
                        else
                        {
                            file.WriteLine("bramka rotacyjna: NIE WYKORZYSTANA");
                        }

                        file.Close();
                    }
                    else
                    {
                        for (int i = 1; i < this.iterations; i++)
                        {
                            double crossProb = this.crossoverProbMax / (1.0 + Convert.ToDouble(i) / Convert.ToDouble(this.iterations));
                            if (crossProb < this.crossoverProbMin) crossProb = this.crossoverProbMin;
                            if (this.selMethodChosen == SelectionMethodChosen.Roulette)
                            {
                                this.Population = new Population(this.selOPerator.RouletteMethod(this.Population));
                            }
                            else
                            {
                                this.Population = new Population(this.selOPerator.RankingMethod(this.Population));
                            }

                            if (this.crossOperChosen == CrossoverOperatorChosen.CX)
                            {
                                this.cxOperator.CrossoverProbability = crossProb;
                                this.Population = new Population(this.cxOperator.Execute(this.Population));
                            }
                            else if (this.crossOperChosen == CrossoverOperatorChosen.OX)
                            {
                                this.oxOperator.CrossoverProbability = crossProb;
                                this.Population = new Population(this.oxOperator.Execute(this.Population));
                            }
                            else
                            {
                                this.pmxOperator.CrossoverProbability = crossProb;
                                this.Population = new Population(this.pmxOperator.Execute(this.Population));
                            }

                            this.mutOperator.Execute(this.Population);

                            if (this.rotGateEnabled) this.rotOperator.Execute(this.Population, this.best);

                            GetBestSolution();

                            if (bestCurrentGoal > this.best.Goal)
                            {
                                bestCurrentGoal = this.best.Goal;
                                bestIteration = i;
                                globalBest = new Solution(this.best);
                            }
                        }
                    }

                    GetBestSolution();

                    this.isStopped = true;

                    Console.Clear();
                    Console.WriteLine("Znalezione rozwiazanie:");
                    this.best.PrintSolution();

                    Console.WriteLine();
                    Console.WriteLine("Wartosc funkcji celu:");
                    Console.WriteLine(this.best.Goal);

                    Console.WriteLine();
                    Console.WriteLine("Najlepsze znalezione rozwiazanie");
                    globalBest.PrintSolution();

                    Console.WriteLine();
                    Console.WriteLine("Wartosc funkcji celu najlepszego rozwiazania:");
                    Console.WriteLine(globalBest.Goal);

                    Console.WriteLine("Znalezione w " + bestIteration + " iteracji.");
                    Console.ReadKey();
                    this.isSet = false;
                }
                else
                {
                    bool OK = false;
                    ConsoleKeyInfo cki;
                    if (this.isStopped == false)
                    {
                        while (!OK)
                        {
                            Console.Clear();
                            Console.WriteLine("Czy chcesz ustawic parametru algorytmu?");
                            Console.WriteLine("1. Tak.");
                            Console.WriteLine("2. Nie. Wyjdz z programu.");
                            cki = Console.ReadKey();
                            if (cki.Key == ConsoleKey.D2)
                            {
                                exit = true;
                                OK = true;
                            }
                            else if (cki.Key == ConsoleKey.D1)
                            {
                                this.isSet = SetParameters(true, false);
                                if (this.isSet)
                                {
                                    InitRandomPopulation();
                                    this.initialPopulationCopy = new Population(this.Population);
                                    OK = true;
                                }
                            }
                        }
                    }
                    else
                    {
                        while (!OK)
                        {
                            Console.Clear();
                            Console.WriteLine("Czy chcesz uruchomic algorytm ponownie?");
                            Console.WriteLine("1. Tak.");
                            Console.WriteLine("2. Nie. Wyjdz z programu.");
                            cki = Console.ReadKey();
                            if (cki.Key == ConsoleKey.D1)
                            {
                                OK = false;
                                while (!OK)
                                {
                                    Console.Clear();
                                    Console.WriteLine("1. Uruchom dla tej samej populacji poczatkowej i tych samych parametrow.");
                                    Console.WriteLine("2. Ustaw nowe parametry i uruchom dla tej samej populacji poczatkowej.");
                                    Console.WriteLine("3. Ustaw nowe parametry i wygeneruj nowa populacje poczatkowa.");
                                    cki = Console.ReadKey();
                                    if (cki.Key == ConsoleKey.D1)
                                    {
                                        this.SetParameters(false, true);
                                        this.Population = new Population(this.initialPopulationCopy);
                                        this.isSet = true;
                                        OK = true;
                                    }
                                    else if (cki.Key == ConsoleKey.D2)
                                    {
                                        this.isSet = SetParameters(false, false);
                                        if (this.isSet)
                                        {
                                            this.Population = new Population(this.initialPopulationCopy);
                                        }
                                        OK = true;
                                    }
                                    else if (cki.Key == ConsoleKey.D3)
                                    {
                                        this.isSet = SetParameters(true, false);
                                        if (this.isSet)
                                        {
                                            InitRandomPopulation();
                                            this.initialPopulationCopy = new Population(this.Population);
                                        }
                                        OK = true;
                                    }
                                }
                            }
                            else if (cki.Key == ConsoleKey.D2)
                            {
                                exit = true;
                                OK = true;
                            }
                        }
                    }
                }
            }
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
    public enum CrossoverOperatorChosen
    {
        CX,
        OX,
        PMX
    }
    public enum SelectionMethodChosen
    {
        Roulette,
        Ranking
    }
}
