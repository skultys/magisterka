using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

namespace magisterka
{
    class Program
    {
        static void Main(string[] args)
        {
            string[] testData = System.IO.File.ReadAllLines(@"D:\Magisterka\had12.dat");

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

            QgAlgorithm alg = new QgAlgorithm(distance, flow, length);
            alg.Execute();
        }
    }
}
