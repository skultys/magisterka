
/////////////////////////////////////////////////////////////////////////////
// C# Version Copyright (c) 2011   Wojciech Chmiel                         //
// All rights reserved                                                     //
/////////////////////////////////////////////////////////////////////////////
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Collections;

namespace OptimisationClassLibrary
{
    public interface IEvolutionAlgorithm : IOptimisationAlgorithm
    {
       IPopulation Population{get;set;}
    }

    public interface IPopulation : IEnumerable
    {
        int Size { get; }
        //static IPopulation operator +(IPopulation pop, ISolution sol);
        int Add(ISolution sol);
        ISolution this[int index] { get; set; }

        /// <summary>
        /// Number of solution in population
        /// </summary>
        int CurrentSize { get; }

        /// <summary>
        /// Remove solution on index position
        /// </summary>
        /// <param name="idx"></param>
        void Remove(int idx);

        void SortAscending();
        void SortDescending();
    }


    public interface IOptimisationAlgorithm
    {
        void InitRandomPopulation();
        /// <summary>
        /// Run algorithm
        /// </summary>
        void Execute();
        /// <summary>
        /// Check stop
        /// </summary>
        /// <returns></returns>
        bool Stop();
        /// <summary>
        /// Return best solution 
        /// </summary>
        /// <returns></returns>
        ISolution GetBestSolution();

        /// <summary>
        /// Get solution criteria value
        /// </summary>
        /// <param name="solution"></param>
        /// <returns></returns>
        double EvaluateSolution(ISolution solution);
    }



    public interface IOptimisationAlgorithmPhyton
    {
    
    
    
    
    
    }
    






    public interface ISolution
    {
        double Goal { get; set; }
        int Size { get; }
    }


    interface IEvolutionaryOperator
    {
    }

    interface IMutationOperator : IEvolutionaryOperator
    {

        ISolution Execute(IPopulation population);
    }

    interface ICrossoverOperator : IEvolutionaryOperator
    {
        ISolution[] Execute(IPopulation population);
    }



  
}
