using System;
using System.Collections.Generic;
using System.Text;
using OptimisationClassLibrary;
using System.Collections;
using System.IO;
using Microsoft.Office.Interop.Excel;

namespace magisterka
{
    public class Qbit
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
            //this.observedState = -1;
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
            //this.observedState = -1;
        }

        public int  EvaluateState()
        {
            double treshold = rand.NextDouble();
            double alphaSquare = this.Alpha * this.Alpha;
            if (treshold > alphaSquare)
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
            this.Beta = Math.Sin(theta) * tempAlpha + Math.Cos(theta) * this.Beta;
            //this.observedState = -1;
        }

        public void ExecuteNotGate()
        {
            double temp = this.Alpha;
            this.Alpha = this.Beta;
            this.Beta = temp;
            //this.observedState = -1;
        }
    }

    public class Chromosome
    {
        private List<Qbit> genes = null;
        private int decodedValue;
        private int size;

        /*public Chromosome()
        {
            genes = new List<Qbit>();
            decodedValue = 0;
        }*/
        
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
            this.size = size;
        }

        public int Size
        {
            get
            {
                return this.size;
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
                //if (this[i].ObservedState == -1) this[i].EvaluateState();
                //chromosomeValue += this[i].ObservedState * (int)Math.Pow(2.0, i);
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

    public class Solution : ISolution
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

            // HAD12
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

            //LIPA50b
            //for (int i = 0; i < this.Size; i++) this[i].PermutationValue = i + 1;

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

    public class Population : IPopulation
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

    public class MutationOperator : IMutationOperator
    {
        private static Random rand = new Random();

        //public int solSize { get; set; }

        //int bitsInSol;

        public double MutationProbability { get; set; }

        public MutationOperator(double probability = 0.0)
        {
            this.MutationProbability = probability;
            //this.solSize = solSize;
            //this.bitsInSol = (int)(Math.Log(solSize, 2.0) + 1);
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
                    int selectedChromosome = rand.Next(0, sol.Size - 1);
                    int bitsInChromosome = (int)(Math.Log(sol.Size, 2.0) + 1);
                    int selectedQbit = rand.Next(0, bitsInChromosome - 1);

                    sol[selectedChromosome][selectedQbit].ExecuteNotGate();
                    sol.toPermutation();
                    solToReturn = new Solution(sol);
                    break;
                }
            }
            return solToReturn;
        }
    }

    public class OxCrossoverOperator : ICrossoverOperator
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
            Solution[] childrenArray = new Solution[2];

            int leftBound = rand.Next(Convert.ToInt32(0.6 * size));
            int rightBound = leftBound + Convert.ToInt32(0.4 * size);
            if (rightBound < size - 1)
            {
                rightBound = rand.Next(rightBound, size - 1);
            }
            else if (rightBound >= size)
            {
                rightBound = size - 1;
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
            /*Console.Clear();
            Console.WriteLine(leftBound);
            Console.WriteLine(rightBound);
            Console.WriteLine();
            for (int i = 0; i < size; i++)
            {
                Console.Write(parentOne[i].PermutationValue + " ");
            }
            Console.WriteLine();
            for (int i = 0; i < size; i++)
            {
                Console.Write(parentTwo[i].PermutationValue + " ");
            }*/

            childrenArray[0] = new Solution(childOne);
            childrenArray[1] = new Solution(childTwo);

            /*Console.WriteLine();
            for (int i = 0; i < size; i++)
            {
                Console.Write(childrenArray[0][i].PermutationValue + " ");
            }
            Console.WriteLine();
            for (int i = 0; i < size; i++)
            {
                Console.Write(childrenArray[1][i].PermutationValue + " ");
            }
            Console.ReadKey();*/
            return childrenArray;
        }
    }


    public class PmxCrossoverOperator : ICrossoverOperator
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

            int leftBound = rand.Next(Convert.ToInt32(0.6 * size));
            int rightBound = leftBound + Convert.ToInt32(0.4 * size);
            if (rightBound < size - 1)
            {
                rightBound = rand.Next(rightBound, size - 1);
            }
            else if (rightBound >= size)
            {
                rightBound = size - 1;
            }
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
            /*Console.Clear();
            Console.WriteLine(leftBound);
            Console.WriteLine(rightBound);
            Console.WriteLine();
            for (int i = 0; i < size; i++)
            {
                Console.Write(parentOne[i].PermutationValue + " ");
            }
            Console.WriteLine();
            for (int i = 0; i < size; i++)
            {
                Console.Write(parentTwo[i].PermutationValue + " ");
            }*/

            childrenArray[0] = new Solution(childOne);
            childrenArray[1] = new Solution(childTwo);

            /*Console.WriteLine();
            for (int i = 0; i < size; i++)
            {
                Console.Write(childrenArray[0][i].PermutationValue + " ");
            }
            Console.WriteLine();
            for (int i = 0; i < size; i++)
            {
                Console.Write(childrenArray[1][i].PermutationValue + " ");
            }
            Console.ReadKey();*/
            return childrenArray;
        }
    }

    public class CxCrossoverOperator : ICrossoverOperator
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

    public class CatastropheOperator : IEvolutionaryOperator
    {
        public double CastastropheProbability { get; set; }

        public void Execute(IPopulation population)
        {
        }
    }

    public class RotationGateOperator : IEvolutionaryOperator
    {
        private int solSize;

        private int bitsInSol;

        private static Random rand = new Random();

        public RotationGateOperator()
        {
        }

        public void ExecuteOriginal(IPopulation population, Solution best)
        {
            double theta;
            double alphaTimesBeta;
            double angle = 0.0;
            double sign = 0.0;
            this.solSize = best.Size;
            this.bitsInSol = (int)(Math.Log(this.solSize, 2.0) + 1);
            foreach (Solution sol in population)
            {
                if (sol.Goal > best.Goal)
                {
                    for (int i = 0; i < this.solSize; i++)
                    {
                        for (int j = 0; j < this.bitsInSol; j++)
                        {
                            theta = 0.0;
                            alphaTimesBeta = sol[i][j].Alpha * sol[i][j].Beta;
                            angle = 0.0;
                            sign = 0.0;
                            if (sol[i][j].ObservedState == 1 && best[i][j].ObservedState == 0)
                            //if (best[i][j].ObservedState == 0)
                            {
                                if (alphaTimesBeta > 0.0) sign = -1.0;
                                else if (alphaTimesBeta < 0.0) sign = 1.0;
                                else if (sol[i][j].Alpha == 0.0)
                                {
                                    double d = rand.NextDouble();
                                    if (d > 0.5) sign = 1;
                                    else sign = -1.0;
                                }
                                angle = 0.5 * Math.PI;
                            }
                            else if (sol[i][j].ObservedState == 1 && best[i][j].ObservedState == 1)
                            //else if (best[i][j].ObservedState == 1)
                            {
                                if (alphaTimesBeta > 0.0) sign = 1.0;
                                else if (alphaTimesBeta < 0.0) sign = -1.0;
                                else if (sol[i][j].Beta == 0.0)
                                {
                                    double d = rand.NextDouble();
                                    if (d > 0.5) sign = 1.0;
                                    else sign = -1.0;
                                }
                                angle = 0.2 * Math.PI;
                            };
                            //if (sol[i][j].ObservedState != sol[i][j].ObservedState) angle = 0.5 * Math.PI;
                            //else if (sol[i][j].ObservedState == best[i][j].ObservedState) angle = 0.2 * Math.PI;

                            theta = angle * sign;

                            if (theta != 0.0) sol[i][j].ExecuteRotationGate(theta);
                        }
                    }
                    sol.toPermutation(); 
                }
            }
        }

        public void ExecuteModified(IPopulation population, Solution best)
        {
            double sign = 0.0;
            double theta = 0.0;
            double alphaTimesBeta = 0.0;
            double smallAngle = 0.001 * Math.PI;
            double bigAngle = 0.08 * Math.PI;
            this.solSize = best.Size;
            this.bitsInSol = (int)(Math.Log(this.solSize, 2.0) + 1);
            foreach (Solution sol in population)
            {
                for (int i = 0; i < this.solSize; i++)
                {
                    for (int j = 0; j < this.bitsInSol; j++)
                    {
                        theta = 0.0;
                        alphaTimesBeta = sol[i][j].Alpha * sol[i][j].Beta;
                        if (sol[i][j].ObservedState == 0 && best[i][j].ObservedState == 1 && sol.Goal > best.Goal)
                        {
                            if (alphaTimesBeta > 0) sign = 1.0;
                            else if (alphaTimesBeta < 0) sign = -1.0;
                            else if (sol[i][j].Alpha == 0.0) sign = 0.0;
                            else if (sol[i][j].Beta == 0.0)
                            {
                                double d = rand.NextDouble();
                                if (d >= 0.5) sign = 1.0;
                                else sign = -1.0;
                            }
                            theta = bigAngle * sign;
                        }
                        else if (sol[i][j].ObservedState == 0 && best[i][j].ObservedState == 1 && sol.Goal < best.Goal)
                        {
                            if (alphaTimesBeta > 0) sign = -1.0;
                            else if (alphaTimesBeta < 0) sign = 1.0;
                            else if (sol[i][j].Beta == 0.0) sign = 0.0;
                            else if (sol[i][j].Alpha == 0.0)
                            {
                                double d = rand.NextDouble();
                                if (d >= 0.5) sign = 1.0;
                                else sign = -1.0;
                            }
                            theta = smallAngle * sign;
                        }
                        else if (sol[i][j].ObservedState == 1 && best[i][j].ObservedState == 0 && sol.Goal > best.Goal)
                        {
                            if (alphaTimesBeta > 0) sign = -1.0;
                            else if (alphaTimesBeta < 0) sign = 1.0;
                            else if (sol[i][j].Beta == 0.0) sign = 0.0;
                            else if (sol[i][j].Alpha == 0.0)
                            {
                                double d = rand.NextDouble();
                                if (d >= 0.5) sign = 1.0;
                                else sign = -1.0;
                            }
                            theta = bigAngle * sign;
                        }
                        else if (sol[i][j].ObservedState == 1 && best[i][j].ObservedState == 0 && sol.Goal < best.Goal)
                        {
                            if (alphaTimesBeta > 0) sign = 1.0;
                            else if (alphaTimesBeta < 0) sign = -1.0;
                            else if (sol[i][j].Alpha == 0.0) sign = 0.0;
                            else if (sol[i][j].Beta == 0.0)
                            {
                                double d = rand.NextDouble();
                                if (d >= 0.5) sign = 1.0;
                                else sign = -1.0;
                            }
                            theta = smallAngle * sign;
                        }
                        if (theta != 0.0) sol[i][j].ExecuteRotationGate(theta);
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
        private int size;

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

        public void setQapData(double[,] distance, double[,] flow, int size)
        {
            this.Distance = distance;
            this.Flow = flow;
            this.size = size;
        }

        public double[,] getDistance()
        {
            return this.Distance;
        }

        public double[,] getFlow()
        {
            return this.Flow;
        }

        public int getSize() 
        {
            return this.size;
        }

        private QapData()
        {
        }
    }

    public class SelectionOperator
    {
        private List<double> distribution = null;
        private static Random rand = new Random();
        private double etaMin;
        private double etaMax;
        public double EtaMax
        {
            get
            {
                return this.etaMax;
            }
            set
            {
                if (value > 2.0 || value < 1.0) this.etaMax = 2.0;
                this.etaMax = value;
                this.etaMin = 2.0 - value;
            }
        }

        public SelectionOperator(double etaMax = 2.0)
        {
            this.etaMax = etaMax;
            if (this.etaMax > 2.0 || this.etaMax < 1.0) this.etaMax = 2.0;
            this.etaMin = 2.0 - this.etaMax;
        }

        public Population RouletteMethod(IPopulation population)
        {
            int size = population.Size;
            List<double> goals = new List<double>();
            for (int i = 0; i < size; i++) goals.Add(population[i].Goal);
            goals.Sort();
            double biggestGoal = goals[size - 1];
            Solution[] chosenSolutions = new Solution[size];
            this.distribution = new List<double>();

            double fitSum = 0.0;
            double currentFitSum = 0.0;

            foreach (Solution sol in population)
            {
                fitSum += (biggestGoal - sol.Goal);
            }

            foreach (Solution sol in population)
            {
                currentFitSum += (biggestGoal - sol.Goal);
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
                double probability = (1.0 / population.Size) * (this.EtaMax - ((this.etaMax - this.etaMin) * (Convert.ToDouble(i) / (population.Size - 1.0))));
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

    public class QgAlgorithm : IEvolutionAlgorithm
    {
        private CrossoverOperatorChosen crossOperChosen;
        private SelectionMethodChosen selMethodChosen;
        private RotationGateVersion rotVersion;
        private OxCrossoverOperator oxOperator = null;
        private CxCrossoverOperator cxOperator = null;
        private PmxCrossoverOperator pmxOperator = null;
        private MutationOperator mutOperator = null;
        private RotationGateOperator rotOperator = null;
        private SelectionOperator selOPerator = null;
        private bool isStopped;
        private int iterations;
        private int popSize;
        private int problemSize;
        private int bestIeration;
        private double crossoverProbMax;
        private double crossoverProbMin;
        private double mutationProb;
        private Solution best = null;
        private Solution initialBest = null;
        private Population initialPopulationCopy;

        public QgAlgorithm()
        {
            this.isStopped = false;
            this.cxOperator = new CxCrossoverOperator();
            this.oxOperator = new OxCrossoverOperator();
            this.pmxOperator = new PmxCrossoverOperator();
            this.selOPerator = new SelectionOperator();
            this.rotOperator = new RotationGateOperator();
            this.mutOperator = new MutationOperator();
        }

        public bool setParameters(
            int popSize,
            int iterations,
            SelectionMethodChosen selMethod,
            double eta,
            CrossoverOperatorChosen crossMethod,
            double crossMax,
            double crossMin,
            double mutProb,
            RotationGateVersion rotVersion,
            bool initNew)
        {
            this.problemSize = QapData.Instance.getSize();
            this.popSize = popSize;
            this.iterations = iterations;
            this.bestIeration = 0;
            this.selMethodChosen = selMethod;
            if (eta != -1.0) this.selOPerator.EtaMax = eta;
            this.crossOperChosen = crossMethod;
            this.crossoverProbMax = crossMax;
            this.crossoverProbMin = crossMin;
            this.mutationProb = mutProb;
            this.mutOperator.MutationProbability = mutProb;
            this.rotVersion = rotVersion;
            if (initNew)
            {
                InitRandomPopulation();
            }
            else this.Population = new Population(this.initialPopulationCopy);
            return true;
        }

        public IPopulation Population { get; set; }

        public double EvaluateSolution(ISolution solution)
        {
            return solution.Goal;
        }

        public void InitRandomPopulation()
        {
            this.Population = new Population(this.popSize, this.problemSize);
            this.initialPopulationCopy = new Population(this.Population);
        }

        public void Execute()
        {
            this.bestIeration = 0;
            double crossProb = this.crossoverProbMax;
            this.best = new Solution(GetBestFromPopulation());
            Solution popBest = new Solution(this.best);
            this.initialBest = new Solution(this.best);
            for (int i = 1; i < this.iterations; i++)
            {
                crossProb = this.crossoverProbMax / (1.0 + Convert.ToDouble(i) / Convert.ToDouble(this.iterations));
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

                if (this.rotVersion == RotationGateVersion.Modified) this.rotOperator.ExecuteModified(this.Population, popBest);
                else if (this.rotVersion == RotationGateVersion.Original) this.rotOperator.ExecuteOriginal(this.Population, popBest);

                popBest = new Solution(GetBestFromPopulation());
                if (this.best.Goal > popBest.Goal)
                {
                    this.best = new Solution(popBest);
                    this.bestIeration = i;
                }
            }

            this.isStopped = true;
        }

        public bool Stop()
        {
            return this.isStopped;
        }

        public ISolution GetBestSolution()
        {
            return this.best;
        }

        public Solution GetBestInitialSoultion()
        {
            return this.initialBest;
        }

        public int GetBestIteration()
        {
            return this.bestIeration;
        }

        private Solution GetBestFromPopulation()
        {
            Solution bestSol = null;
            double currentGoal = -1.0;
            foreach (Solution sol in this.Population)
            {
                if (currentGoal == -1 || sol.Goal < currentGoal)
                {
                    currentGoal = sol.Goal;
                    bestSol = new Solution(sol);
                }
            }
            return bestSol;
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

    public enum RotationGateVersion
    {
        None,
        Original,
        Modified
    }

    public enum TestedParameter
    {
        selectionOperator,
        crossoverProb,
        crossoverOperator,
        mutationProb,
        rotationOperator
    }

    public class UserInterface
    {
        public static void RunTest(QgAlgorithm algorithm)
        {
            Application excelApp = null;
            Workbook wb = null;
            Worksheet ws = null;
            string path = Directory.GetCurrentDirectory() + "\\Results";
            if (!Directory.Exists(path)) Directory.CreateDirectory(path);

            path = Directory.GetCurrentDirectory() + "\\QAPLib";
            if (!Directory.Exists(path))
            {
                Directory.CreateDirectory(path);
                Console.Clear();
                Console.WriteLine("Utworzono folder z instancjami testowymi");
                Console.WriteLine("Umiesc w nim instancje testowe i uruchom program ponownie.");
                Console.ReadKey();
                return;
            }

            ConsoleKeyInfo cki;
            double bestInitialGoal = 0.0;
            double bestGoal = 0.0;
            double avGoal = 0.0;
            double avIter = 0.0;
            int bestIter = 0;
            int popSize = 0;
            int iterations = 0;
            int repeatTest = 0;
            int instances = 0;
            Solution bestSol = null;
            TestedParameter testedParameter = TestedParameter.selectionOperator;
            int valuesCount = 0;
            string inputString;
            string[] instancePaths = new string[0];
            SelectionMethodChosen[] selOperArray = new SelectionMethodChosen[0];
            double[] etaArray = new double[0];
            double[] maxCrossProbArray = new double[0];
            double[] minCrossProbArray = new double[0];
            CrossoverOperatorChosen[] crossOperArray = new CrossoverOperatorChosen[0];
            double[] mutProbArray = new double[0];
            RotationGateVersion[] rotOperArray = new RotationGateVersion[0];
            string templatePath = null;


            bool OK = false;

            path = Directory.GetCurrentDirectory() + "\\QAPLib";
            string[] filePaths = Directory.GetFiles(path, "*.dat");
            if (filePaths.Length == 0)
            {
                Console.Clear();
                Console.WriteLine("Folder z instancjami testowymi jest pusty!");
                Console.ReadKey();
                return;
            }

            bool exit = false;
            excelApp = new Application();
            while (!exit)
            {
                while (!OK)
                {
                    Console.Clear();
                    Console.WriteLine("Wybierz parametr algorytmu, ktory chcesz testowac:");
                    Console.WriteLine("1. Operator selekcji.");
                    Console.WriteLine("2. Prawdopodobienstwo krzyzowania.");
                    Console.WriteLine("3. Operator krzyzowania");
                    Console.WriteLine("4. Prawdopodbienstwo mutacji");
                    Console.WriteLine("5. Bramke kwantowa");
                    cki = Console.ReadKey();
                    if (cki.Key == ConsoleKey.D1)
                    {
                        testedParameter = TestedParameter.selectionOperator;
                        OK = true;
                    }
                    else if (cki.Key == ConsoleKey.D2)
                    {
                        testedParameter = TestedParameter.crossoverProb;
                        OK = true;
                    }
                    else if (cki.Key == ConsoleKey.D3)
                    {
                        testedParameter = TestedParameter.crossoverOperator;
                        OK = true;
                    }
                    else if (cki.Key == ConsoleKey.D4)
                    {
                        testedParameter = TestedParameter.mutationProb;
                        OK = true;
                    }
                    else if (cki.Key == ConsoleKey.D5)
                    {
                        testedParameter = TestedParameter.rotationOperator;
                        OK = true;
                    }
                }

                OK = false;
                while (!OK)
                {
                    try
                    {
                        Console.Clear();
                        Console.WriteLine("Ile razy powtorzyc test dla wybranego parametru?");
                        inputString = Console.ReadLine();
                        repeatTest = Convert.ToInt32(inputString);
                        if (repeatTest > 0)
                        {
                            OK = true;
                        }
                    }
                    catch (Exception e)
                    {
                    }
                }

                OK = false;
                while (!OK)
                {
                    try
                    {
                        Console.Clear();
                        Console.WriteLine("Ile razy wybrany parametr ma sie zmieniac?");
                        inputString = Console.ReadLine();
                        valuesCount = Convert.ToInt32(inputString);
                        if (valuesCount > 0)
                        {
                            selOperArray = new SelectionMethodChosen[valuesCount];
                            etaArray = new double[valuesCount];
                            maxCrossProbArray = new double[valuesCount];
                            minCrossProbArray = new double[valuesCount];
                            crossOperArray = new CrossoverOperatorChosen[valuesCount];
                            mutProbArray = new double[valuesCount];
                            rotOperArray = new RotationGateVersion[valuesCount];
                            OK = true;
                        }
                    }
                    catch (Exception e)
                    {
                    }
                }

                for (int i = 0; i < valuesCount; i++)
                {
                    OK = false;
                    try
                    {
                        while (!OK)
                        {
                            if (testedParameter == TestedParameter.selectionOperator)
                            {
                                Console.Clear();
                                Console.WriteLine("Podaj operator selekcji nr " + (i + 1) + ":");
                                Console.WriteLine("1. Selekcja ruletkowa.");
                                Console.WriteLine("2. Selekcja rankingowa");
                                cki = Console.ReadKey();
                                if (cki.Key == ConsoleKey.D1)
                                {
                                    selOperArray[i] = SelectionMethodChosen.Roulette;
                                    etaArray[i] = -1.0;
                                    OK = true;
                                }
                                else if (cki.Key == ConsoleKey.D2)
                                {
                                    while (!OK)
                                    {
                                        try
                                        {
                                            Console.Clear();
                                            Console.WriteLine("Podaj parametr eta (1,0 - 2,0) nr " + (i + 1) + ": ");
                                            inputString = Console.ReadLine();
                                            etaArray[i] = Convert.ToDouble(inputString);
                                            if (etaArray[i] >= 1.0 && etaArray[i] <= 2.0)
                                            {
                                                selOperArray[i] = SelectionMethodChosen.Ranking;
                                                OK = true;
                                            }
                                        }
                                        catch (Exception e)
                                        {
                                        }
                                    }
                                }
                            }
                            else if (testedParameter == TestedParameter.crossoverProb)
                            {
                                while (!OK)
                                {
                                    try
                                    {
                                        Console.Clear();
                                        Console.WriteLine("Podaj maksymalne prawdopodobienstwo krzyzowania (0,0 - 1,0) nr " + (i + 1) + ":");
                                        inputString = Console.ReadLine();
                                        maxCrossProbArray[i] = Convert.ToDouble(inputString);
                                        if (maxCrossProbArray[i] >= 0.0 && maxCrossProbArray[i] <= 1.0)
                                        {
                                            while (!OK)
                                            {
                                                try
                                                {
                                                    Console.Clear();
                                                    Console.WriteLine("Podaj minimalne prawdopodobienstwo krzyzowania (0,0 - " + maxCrossProbArray[i] + ") nr " + (i + 1) + ": ");
                                                    inputString = Console.ReadLine();
                                                    minCrossProbArray[i] = Convert.ToDouble(inputString);
                                                    if (minCrossProbArray[i] >= 0.0 && minCrossProbArray[i] <= maxCrossProbArray[i])
                                                    {
                                                        OK = true;
                                                    }
                                                }
                                                catch (Exception e)
                                                {
                                                }
                                            }
                                        }
                                    }
                                    catch (Exception e)
                                    {
                                    }
                                }
                            }
                            else if (testedParameter == TestedParameter.crossoverOperator)
                            {
                                while (!OK)
                                {
                                    Console.Clear();
                                    Console.WriteLine("Podaj operator krzyzowania nr " + (i + 1) + ":");
                                    Console.WriteLine("1. CX.");
                                    Console.WriteLine("2. OX.");
                                    Console.WriteLine("3. PMX.");
                                    cki = Console.ReadKey();
                                    if (cki.Key == ConsoleKey.D1)
                                    {
                                        crossOperArray[i] = CrossoverOperatorChosen.CX;
                                        OK = true;
                                    }
                                    else if (cki.Key == ConsoleKey.D2)
                                    {
                                        crossOperArray[i] = CrossoverOperatorChosen.OX;
                                        OK = true;
                                    }
                                    else if (cki.Key == ConsoleKey.D3)
                                    {
                                        crossOperArray[i] = CrossoverOperatorChosen.PMX;
                                        OK = true;
                                    }
                                }
                            }
                            else if (testedParameter == TestedParameter.mutationProb)
                            {
                                while (!OK)
                                {
                                    try
                                    {
                                        Console.Clear();
                                        Console.WriteLine("Podaj prawdopodobienstwo mutacji (0,0 - 1,0) nr " + (i + 1) + ": ");
                                        inputString = Console.ReadLine();
                                        mutProbArray[i] = Convert.ToDouble(inputString);
                                        if (mutProbArray[i] >= 0.0 && mutProbArray[i] <= 1.0)
                                        {
                                            OK = true;
                                        }
                                    }
                                    catch (Exception e)
                                    {
                                    }
                                }
                            }
                            else if (testedParameter == TestedParameter.rotationOperator)
                            {
                                while (!OK)
                                {
                                    Console.Clear();
                                    Console.WriteLine("Podaj wersje operatora bramki kwantowej nr " + (i + 1) + ":");
                                    Console.WriteLine("1. Brak operatora.");
                                    Console.WriteLine("2. Oryginalna.");
                                    Console.WriteLine("3. Zmodyfikowana.");
                                    cki = Console.ReadKey();
                                    if (cki.Key == ConsoleKey.D1)
                                    {
                                        rotOperArray[i] = RotationGateVersion.None;
                                        OK = true;
                                    }
                                    else if (cki.Key == ConsoleKey.D2)
                                    {
                                        rotOperArray[i] = RotationGateVersion.Original;
                                        OK = true;
                                    }
                                    else if (cki.Key == ConsoleKey.D3)
                                    {
                                        rotOperArray[i] = RotationGateVersion.Modified;
                                        OK = true;
                                    }
                                }
                            }
                        }
                    }
                    catch (Exception e)
                    {
                    }
                }

                OK = false;
                while (!OK)
                {
                    try
                    {
                        Console.Clear();
                        Console.WriteLine("Dla ilu instancji testowych przeprowadzic testy?");
                        inputString = Console.ReadLine();
                        instances = Convert.ToInt32(inputString);
                        if (instances > 0)
                        {
                            instancePaths = new string[instances];
                            for (int i = 0; i < instances; i++)
                            {
                                while (!OK)
                                {
                                    try
                                    {
                                        Console.Clear();
                                        Console.WriteLine("Wybierz instancje testowa nr " + (i + 1) + ":");
                                        Console.WriteLine();
                                        for (int j = 0; j < filePaths.Length; j++)
                                        {
                                            filePaths[j] = Path.GetFileName(filePaths[j]);
                                            Console.WriteLine((j + 1) + ". " + filePaths[j]);
                                        }
                                        inputString = Console.ReadLine();
                                        int fileNumber = Convert.ToInt32(inputString);
                                        if (fileNumber > 0 && fileNumber <= filePaths.Length)
                                        {
                                            instancePaths[i] = filePaths[fileNumber - 1];
                                            OK = true;
                                        }
                                    }
                                    catch (Exception e)
                                    {
                                    }
                                }
                                OK = false;
                            }
                            OK = true;
                        }
                    }
                    catch (Exception e)
                    {
                    }
                }

                Console.Clear();
                Console.WriteLine("Ustawianie pozostalych parametrow, niezmiennych w kazdym z testow.");
                Console.ReadKey();

                OK = false;
                while (!OK)
                {
                    try
                    {
                        Console.Clear();
                        Console.WriteLine("Podaj rozmiar populacji: ");
                        inputString = Console.ReadLine();
                        popSize = Convert.ToInt32(inputString);
                        if (popSize > 0)
                        {
                            if (popSize % 2 != 0) popSize += 1;
                            OK = true;
                        }
                    }
                    catch (Exception e)
                    {
                    }
                }
                OK = false;

                while (!OK)
                {
                    try
                    {
                        Console.Clear();
                        Console.WriteLine("Podaj ilosc iteracji: ");
                        inputString = Console.ReadLine();
                        iterations = Convert.ToInt32(inputString);
                        if (iterations > 0)
                        {
                            OK = true;
                        }
                    }
                    catch (Exception e)
                    {
                    }
                }
                OK = false;

                while (!OK && testedParameter != TestedParameter.selectionOperator)
                {
                    Console.Clear();
                    Console.WriteLine("Wybierz metode selekcji: ");
                    Console.WriteLine("1. Ruletkowa");
                    Console.WriteLine("2. Rankingowa");
                    cki = Console.ReadKey();
                    if (cki.Key == ConsoleKey.D1)
                    {
                        for (int i = 0; i < valuesCount; i++)
                        {
                            selOperArray[i] = SelectionMethodChosen.Roulette;
                            etaArray[i] = -1.0;
                        }
                        OK = true;
                    }
                    else if (cki.Key == ConsoleKey.D2)
                    {
                        try
                        {
                            while (!OK)
                            {
                                Console.Clear();
                                Console.WriteLine("Podaj parametr eta (1,0 - 2,0): ");
                                inputString = Console.ReadLine();
                                etaArray[0] = Convert.ToDouble(inputString);
                                if (etaArray[0] >= 1.0 && etaArray[0] <= 2.0)
                                {
                                    for (int i = 1; i < valuesCount; i++)
                                    {
                                        selOperArray[i] = SelectionMethodChosen.Ranking;
                                        etaArray[i] = etaArray[0];
                                    }
                                    OK = true;
                                }
                            }
                        }
                        catch (Exception e)
                        {
                        }
                    }
                }
                OK = false;

                while (!OK && testedParameter != TestedParameter.crossoverProb)
                {
                    try
                    {
                        Console.Clear();
                        Console.WriteLine("Podaj maksymalne prawdopodobienstwo krzyzowania (0,0 - 1,0): ");
                        inputString = Console.ReadLine();
                        maxCrossProbArray[0] = Convert.ToDouble(inputString);
                        if (maxCrossProbArray[0] >= 0.0 && maxCrossProbArray[0] <= 1.0)
                        {
                            for (int i = 1; i < valuesCount; i++)
                            {
                                maxCrossProbArray[i] = maxCrossProbArray[0];
                            }
                            while (!OK)
                            {
                                try
                                {
                                    Console.Clear();
                                    Console.WriteLine("Podaj minimalne prawdopodobienstwo krzyzowania (0,0 - " + maxCrossProbArray[0] + "): ");
                                    inputString = Console.ReadLine();
                                    minCrossProbArray[0] = Convert.ToDouble(inputString);
                                    if (minCrossProbArray[0] >= 0.0 && minCrossProbArray[0] <= maxCrossProbArray[0])
                                    {
                                        for (int i = 1; i < valuesCount; i++)
                                        {
                                            minCrossProbArray[i] = minCrossProbArray[0];
                                        }
                                        OK = true;
                                    }
                                }
                                catch (Exception e)
                                {
                                }
                            }
                            OK = true;
                        }
                    }
                    catch (Exception e)
                    {
                    }
                }
                OK = false;

                while (!OK && testedParameter != TestedParameter.crossoverOperator)
                {
                    Console.Clear();
                    Console.WriteLine("Wybierz operator krzyzowania: ");
                    Console.WriteLine("1. CX");
                    Console.WriteLine("2. OX");
                    Console.WriteLine("3. PMX");
                    cki = Console.ReadKey();
                    if (cki.Key == ConsoleKey.D1)
                    {
                        for (int i = 0; i < valuesCount; i++)
                        {
                            crossOperArray[i] = CrossoverOperatorChosen.CX;
                        }
                        OK = true;
                    }
                    else if (cki.Key == ConsoleKey.D2)
                    {
                        for (int i = 0; i < valuesCount; i++)
                        {
                            crossOperArray[i] = CrossoverOperatorChosen.OX;
                        }
                        OK = true;
                    }
                    else if (cki.Key == ConsoleKey.D3)
                    {
                        for (int i = 0; i < valuesCount; i++)
                        {
                            crossOperArray[i] = CrossoverOperatorChosen.PMX;
                        }
                        OK = true;
                    }
                }
                OK = false;

                while (!OK && testedParameter != TestedParameter.mutationProb)
                {
                    try
                    {
                        Console.Clear();
                        Console.WriteLine("Podaj prawdopodobienstwo mutacji (0.0 - 1.0): ");
                        inputString = Console.ReadLine();
                        mutProbArray[0] = Convert.ToDouble(inputString);
                        if (mutProbArray[0] >= 0.0 && mutProbArray[0] <= 1.0)
                        {
                            for (int i = 1; i < valuesCount; i++)
                            {
                                mutProbArray[i] = mutProbArray[0];
                            }
                            OK = true;
                        }
                    }
                    catch (Exception e)
                    {
                    }
                }
                OK = false;

                while (!OK && testedParameter != TestedParameter.rotationOperator)
                {
                    Console.Clear();
                    Console.WriteLine("Czy chcesz uzyc operator bramki rotacyjnej?");
                    Console.WriteLine("1. Tak.");
                    Console.WriteLine("2. Nie.");
                    cki = Console.ReadKey();
                    if (cki.Key == ConsoleKey.D1)
                    {
                        while (!OK)
                        {
                            Console.Clear();
                            Console.WriteLine("Ktora wersje bramki rotacyjnej chcesz uzyc?");
                            Console.WriteLine("1. Oryginalna");
                            Console.WriteLine("2. Zmodyfikowana");
                            cki = Console.ReadKey();
                            if (cki.Key == ConsoleKey.D1)
                            {
                                for (int i = 0; i < valuesCount; i++)
                                {
                                    rotOperArray[i] = RotationGateVersion.Original;
                                }
                                OK = true;
                            }
                            else if (cki.Key == ConsoleKey.D2)
                            {
                                for (int i = 0; i < valuesCount; i++)
                                {
                                    rotOperArray[i] = RotationGateVersion.Modified;
                                }
                                OK = true;
                            }
                        }
                    }
                    else if (cki.Key == ConsoleKey.D2)
                    {
                        for (int i = 0; i < valuesCount; i++)
                        {
                            rotOperArray[i] = RotationGateVersion.None;
                        }
                        OK = true;
                    }
                }

                OK = false;
                while (!OK)
                {
                    Console.Clear();
                    Console.WriteLine("Podaj szablon nazwy pliku, do ktorego beda zapisane dane.");
                    Console.WriteLine("Zostana utworzone pliki o nastepujacej nazwie: {szablon}_{parametr testowany}");
                    templatePath = Console.ReadLine();
                    if (templatePath.Length != 0)
                    {
                        OK = true;
                    }
                }

                //testowanie
                Console.Clear();
                Console.WriteLine("Rozpoczeto testy...");
                string resultsPath = "";
                bool initNew;
                excelApp.SheetsInNewWorkbook = instances;
                if (testedParameter == TestedParameter.selectionOperator)
                {
                    resultsPath = templatePath + "_selOper";
                }
                else if (testedParameter == TestedParameter.crossoverProb)
                {
                    resultsPath = templatePath + "_crossProb";
                }
                else if (testedParameter == TestedParameter.crossoverOperator)
                {
                    resultsPath = templatePath + "_crossOper";
                }
                else if (testedParameter == TestedParameter.mutationProb)
                {
                    resultsPath = templatePath + "_mutProb";
                }
                else if (testedParameter == TestedParameter.rotationOperator)
                {
                    resultsPath = templatePath + "_rotOper";
                }

                resultsPath = Directory.GetCurrentDirectory() + "\\Results\\" + resultsPath;
                if (File.Exists(resultsPath + ".xlsx"))
                {
                    File.Delete(resultsPath + ".xlsx");
                }
                else if (File.Exists(resultsPath + ".xls"))
                {
                    File.Delete(resultsPath + ".xls");
                }
                wb = excelApp.Workbooks.Add();
                wb.SaveAs(resultsPath);
                System.IO.FileInfo fileInfo = new System.IO.FileInfo(wb.FullName);
                fileInfo.IsReadOnly = false;
                int skip = 0;
                double[] bestvalues = {5426670, 7990, 1534, 3796, 116, 1652, 88700, 3683, 1210244, 253195, 360630};
                for (int i = 0; i < instances; i++)
                {
                    skip = 0;
                    ws = wb.Worksheets[i + 1];
                    ws.Name = instancePaths[i].Substring(0, instancePaths[i].Length - 4);
                    ws.Cells[1, 2] = "best init";
                    ws.Cells[1, 3] = "best";
                    ws.Cells[1, 4] = "av best";
                    ws.Cells[1, 5] = "best iter";
                    ws.Cells[1, 6] = "av iter";
                    ws.Cells[1, 7] = bestvalues[i];
                    UserInterface.ReadQapFile(instancePaths[i]);

                    ws.Cells[valuesCount + 3, 1] = "iter";
                    ws.Cells[valuesCount + 4, 1] = iterations;
                    ws.Cells[valuesCount + 3, 2] = "pop size";
                    ws.Cells[valuesCount + 4, 2] = popSize;
                    if (testedParameter == TestedParameter.crossoverProb)
                    {
                        skip = 2;
                    }
                    else
                    {
                        ws.Cells[valuesCount + 3, 3] = "min cross prob";
                        ws.Cells[valuesCount + 4, 3] = minCrossProbArray[0];
                        ws.Cells[valuesCount + 3, 4] = "max cross prob";
                        ws.Cells[valuesCount + 4, 4] = maxCrossProbArray[0];
                    }
                    if (testedParameter == TestedParameter.crossoverOperator)
                    {
                        skip = 1;
                    }
                    else
                    {
                        ws.Cells[valuesCount + 3, 5 - skip] = "cross oper";
                        if (crossOperArray[0] == CrossoverOperatorChosen.CX)
                        {
                            ws.Cells[valuesCount + 4, 5 - skip] = "CX";
                        }
                        else if (crossOperArray[0] == CrossoverOperatorChosen.OX)
                        {
                            ws.Cells[valuesCount + 4, 5 - skip] = "OX";
                        }
                        else
                        {
                            ws.Cells[valuesCount + 4, 5 - skip] = "PMX";
                        }
                    }
                    if (testedParameter == TestedParameter.mutationProb)
                    {
                        skip = 1;
                    }
                    else
                    {
                        ws.Cells[valuesCount + 3, 6 - skip] = "mut prob";
                        ws.Cells[valuesCount + 4, 6 - skip] = mutProbArray[0];
                    }
                    if (testedParameter == TestedParameter.rotationOperator)
                    {
                        skip = 1;
                    }
                    else
                    {
                        ws.Cells[valuesCount + 3, 7 - skip] = "rot gate";
                        if (rotOperArray[0] == RotationGateVersion.Modified)
                        {
                            ws.Cells[valuesCount + 4, 7 - skip] = "mod";
                        }
                        else if (rotOperArray[0] == RotationGateVersion.Original)
                        {
                            ws.Cells[valuesCount + 4, 7 - skip] = "orig";
                        }
                        else
                        {
                            ws.Cells[valuesCount + 4, 7 - skip] = "none";
                        }
                    }
                    if (testedParameter == TestedParameter.selectionOperator)
                    {
                        skip = 1;
                    }
                    else
                    {
                        ws.Cells[valuesCount + 3, 8 - skip] = "sel oper";
                        if (selOperArray[0] == SelectionMethodChosen.Ranking)
                        {
                            skip = -1;
                            ws.Cells[valuesCount + 4, 8 - skip] = "rank";
                            ws.Cells[valuesCount + 3, 9 - skip] = "eta";
                            ws.Cells[valuesCount + 4, 9 - skip] = etaArray[0];
                        }
                        else
                        {
                            ws.Cells[valuesCount + 4, 8 - skip] = "roull";
                        }
                    }
                    ws.Cells[valuesCount + 3, 9 - skip] = "test count";
                    ws.Cells[valuesCount + 4, 9 - skip] = repeatTest;
                    for (int j = 0; j < valuesCount; j++)
                    {
                        if (testedParameter == TestedParameter.selectionOperator)
                        {
                            if (selOperArray[j] == SelectionMethodChosen.Ranking)
                            {
                                ws.Cells[j + 2, 1] = "rank_eta_" + etaArray[j];
                            }
                            else
                            {
                                ws.Cells[j + 2, 1] = "roull";
                            }
                        }
                        else if (testedParameter == TestedParameter.crossoverProb)
                        {
                            ws.Cells[j + 2, 1] = "min_" + minCrossProbArray[j] + "_max_" + maxCrossProbArray[j];
                        }
                        else if (testedParameter == TestedParameter.crossoverOperator)
                        {
                            if (crossOperArray[j] == CrossoverOperatorChosen.CX)
                            {
                                ws.Cells[j + 2, 1] = "CX";
                            }
                            else if (crossOperArray[j] == CrossoverOperatorChosen.OX)
                            {
                                ws.Cells[j + 2, 1] = "OX";
                            }
                            else
                            {
                                ws.Cells[j + 2, 1] = "PMX";
                            }
                        }
                        else if (testedParameter == TestedParameter.mutationProb)
                        {
                            ws.Cells[j + 2, 1] = mutProbArray[j].ToString();
                        }
                        else if (testedParameter == TestedParameter.rotationOperator)
                        {
                            if (rotOperArray[j] == RotationGateVersion.Modified)
                            {
                                ws.Cells[j + 2, 1] = "mod";
                            }
                            else if (rotOperArray[j] == RotationGateVersion.Original)
                            {
                                ws.Cells[j + 2, 1] = "orig";
                            }
                            else
                            {
                                ws.Cells[j + 2, 1] = "none";
                            }
                        }

                        bestInitialGoal = 0.0;
                        bestGoal = 0.0;
                        avGoal = 0.0;
                        bestIter = 0;
                        avIter = 0;
                        for (int k = 0; k < repeatTest; k++)
                        {
                            if (j == 0 && k == 0)
                            {
                                initNew = true;
                            }
                            else
                            {
                                initNew = false;
                            }
                            if (algorithm.setParameters(popSize, iterations, selOperArray[j], etaArray[j], crossOperArray[j], maxCrossProbArray[j], minCrossProbArray[j], mutProbArray[j], rotOperArray[j], initNew))
                            {
                                initNew = false;
                                Console.Write("Instancja nr {0}, wartosc parametru nr {1}, test nr {2}...", (i + 1), (j + 1), (k + 1));
                                algorithm.Execute();
                                avGoal += algorithm.GetBestSolution().Goal;
                                avIter += algorithm.GetBestIteration();
                                if (k == 0)
                                {
                                    bestIter = algorithm.GetBestIteration();
                                    bestInitialGoal = algorithm.GetBestInitialSoultion().Goal;
                                    bestGoal = algorithm.GetBestSolution().Goal;
                                    bestSol = new Solution((Solution)algorithm.GetBestSolution());
                                }
                                else
                                {
                                    if (bestIter > algorithm.GetBestIteration())
                                    {
                                        bestIter = algorithm.GetBestIteration();
                                    }
                                    if (bestGoal > algorithm.GetBestSolution().Goal)
                                    {
                                        bestGoal = algorithm.GetBestSolution().Goal;
                                        bestSol = new Solution((Solution)algorithm.GetBestSolution());
                                    }
                                }
                                Console.Write("OK");
                                Console.WriteLine();
                            }
                            else
                            {
                                wb.Save();
                                wb.Close();
                                excelApp.Quit();
                                Console.Clear();
                                Console.WriteLine("Blad podczas testu");
                                while (System.Runtime.InteropServices.Marshal.ReleaseComObject(excelApp) != 0);
                                ws = null;
                                wb = null;
                                excelApp = null;
                                return;
                            }
                        }
                        avGoal /= Convert.ToDouble(repeatTest);
                        avIter /= Convert.ToDouble(repeatTest);
                        ws.Cells[j + 2, 2] = bestInitialGoal;
                        ws.Cells[j + 2, 3] = bestGoal;
                        ws.Cells[j + 2, 4] = avGoal;
                        ws.Cells[j + 2, 5] = bestIter;
                        ws.Cells[j + 2, 6] = avIter;
                        ws.Cells[j + 2, 7] = (bestGoal - bestvalues[i]) / bestvalues[i];
                        for (int k = 0; k < bestSol.Size; k++)
                        {
                            ws.Cells[j + 2, k + 8] = bestSol[k].PermutationValue;
                        }
                        wb.Save();
                    }
                }
                OK = false;
                Console.Clear();
                Console.WriteLine("Testy ukonczone pomyslnie.");
                Console.WriteLine("Uruchomic testy ponownie?");
                Console.WriteLine("1. Tak");
                Console.WriteLine("2. Nie. Wyjdz z programu.");
                cki = Console.ReadKey();
                if (cki.Key == ConsoleKey.D2)
                {
                    exit = true;
                    wb.Close();
                    excelApp.Quit();
                    while (System.Runtime.InteropServices.Marshal.ReleaseComObject(excelApp) != 0) ;
                    ws = null;
                    wb = null;
                    excelApp = null;
                }
            }
        }

        public static bool ReadQapFile(string path)
        {
            try
            {
                string fullPath = Directory.GetCurrentDirectory() + "\\QAPLib\\" + path;

                string[] testData = System.IO.File.ReadAllLines(@fullPath);

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
                QapData.Instance.setQapData(distance, flow, length);
                return true;
            }
            catch (Exception e)
            {
                return false;
            }
        }
    }

    public class CtrlcEvent
    {
        private Application excelApp;
        //private Workbook wb;
        public CtrlcEvent(Application app)
        {
            this.excelApp = app;
        }
        public void myHandler(object sender, ConsoleCancelEventArgs args)
        {
            this.excelApp.Quit();
            System.Runtime.InteropServices.Marshal.ReleaseComObject(excelApp);
        }
    }
}
