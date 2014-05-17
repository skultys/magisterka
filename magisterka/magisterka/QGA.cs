using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using OptimisationClassLibrary;
using System.Collections;

namespace magisterka
{
    class Solution : ISolution
    {
        public double Goal { set; get; }
    }

    class Population : IPopulation {
        public int Size
        {
            get
            {
                return 0;
            }
        }
        public int Add(ISolution sol)
        {
            return 0;
        }
        public ISolution this[int idex]
        {
            get
            {
                return new Solution();
            }
            set
            {
            }
        }
        public int CurrentSize
        {
            get
            {
                return 0;
            }
        }
        public void Remove(int idx)
        {
        }

        public IEnumerator GetEnumerator()
        {
            foreach (object o in this)
            {
                if (o == null)
                {
                    break;
                }
                yield return o;
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

    class QGA : IEvolutionAlgorithm
    {
        public IPopulation Population { get; set; }
        public double EvaluateSolution(ISolution solution)
        {
            return 0.0;
        }
        public void InitRandomPopulation() {
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
            return new Solution();
        }

    }
}
