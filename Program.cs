using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using BEGA.Classes;
using GeneticAlgorithm;

namespace BEGA {
    static class Globals {
        public static int Runs = 1;
        public static int MaxGens = 1000000;
        public static int ReportRate = MaxGens / 10;
        public static int MaxProblemSize = 100;
        public static int P = 90;
        public static int E = 15;
        public static int ConvergenceGen = 20000;
    }
    class Program {

        static void WriteToCsv(string file, string[] data) {
            string row = String.Join(", ", data);
            using (FileStream fs = new FileStream(file,FileMode.Append, FileAccess.Write))
            {
                using (StreamWriter sw = new StreamWriter(fs))
                {
                    sw.WriteLine(row);
                }
            }
        }

        static void Main(string[] args) {
            // Get the list of problems
            string path = System.IO.Path.GetDirectoryName(        System.Reflection.Assembly.GetExecutingAssembly().GetName().CodeBase );
            path = path.Replace("file:\\","");
            string instancePath = path + @"\instances\";
            string csvPath = path + @"\csv\";
            
            // Get all the files in instances
            string[] files = new string[0];
            try {
                files = Directory.GetFiles(instancePath);
            }
            catch (Exception e) {
                Console.Error.WriteLine($"[ERROR]: {e.ToString()}");
                Console.WriteLine("Exiting program");
                return;
            }

            // Get the TSP files
            List<string> instances = new List<string>();
            foreach (var f in files) {
                if (f.Contains(".tsp")) {
                    string[] parts = f.Split("\\");
                    string name = parts[parts.Length - 1].Replace(".tsp", "");
                    instances.Add(name);
                }
            }

            if (instances.Count == 0) {
                Console.WriteLine("No instances found");
                Console.Write("Exiting program");
                return;
            }
            Console.WriteLine($"Instances found: {instances.Count}");
            
            // Run each instance
            foreach (var instanceName in instances) {
                Instance instance = new Instance(instancePath, instanceName, 1, Globals.MaxProblemSize);
                if (!instance.IsComplete) {
                    continue;
                }
                if (Globals.Instance.Length > 0 && instance.Name.ToLower() != Globals.Instance.ToLower()) {
                    continue;
                }

                if (Globals.Instance.Length == 0 && instance.Dimension > Globals.MaxProblemSize) {
                    continue;
                }

                Console.WriteLine("-----------------------------------------------------------");
                Console.WriteLine($"I: {instance.Name} | N: {instance.Nodes.Length} | O: {instance.Optimum}");
                
                Stopwatch stopwatch = new Stopwatch();
                for (int r = 0; r < Globals.Runs; r++) {
                    var ga = new Classes.BEGA_MSGM(instance, Globals.P, Globals.E);
                    int lastImprovement = 0;
                    int lastCost = ga.P[0].Cost;
                    stopwatch.Reset();
                    stopwatch.Start();
                    for (int i = 1; i <= Globals.MaxGens; i++) {
                        ga.Evolve();
                        if (i == 1 || i % Globals.ReportRate == 0) {
                            double error = (ga.P[0].Cost - instance.Optimum) / (double) instance.Optimum;
                            var time = stopwatch.Elapsed.ToString("hh\\:mm\\:ss");
                            Console.WriteLine(
                                $"R: {r + 1}/{Globals.Runs} | G: {i.ToString("00000")} | C: {ga.P[0].Cost} | O: {instance.Optimum} | E: {error.ToString("0.0000")} | T: {time}");
                        }
                        if(ga.P[0].Cost < lastCost){
                            lastCost = ga.P[0].Cost;
                            lastImprovement = 0;
                        } else {
                            lastImprovement++;
                        }
                        if(lastImprovement > Globals.ConvergenceGen){
                            // End criteria met. End run.
                            break;
                        }
                    }
                    stopwatch.Stop();
                    double error = (ga.P[0].Cost - instance.Optimum) / (double) instance.Optimum;
                    var time = stopwatch.Elapsed.ToString("hh\\:mm\\:ss");
                    Console.WriteLine(
                        $"R: {r + 1}/{Globals.Runs} | G: {i.ToString("00000")} | C: {ga.P[0].Cost} | O: {instance.Optimum} | E: {error.ToString("0.0000")} | T: {time}");

                    // Compile routes into a string array
                    List<string> dataList = new List<string>();
                    dataList.Add(cost.ToString());
                    foreach (var route in routes) {
                        string routeString = String.Join(" ", route);
                        dataList.Add(routeString);
                    }
                    // Write to csv
                    WriteToCsv(csvPath + $"{name}.csv", dataList.ToArray());
                }
            }
            
            Console.WriteLine($"[INFO] Completed all instances");
        }
    }
}
