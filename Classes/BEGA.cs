using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using GeneticAlgorithm;

namespace BEGA.Classes {
    /// <summary>
    /// P: Population
    /// S_G: Diversity subpopulation
    /// S_E: Intensification subpopulation
    /// P_e: Elite subpopulation
    /// M_d: Minimum distance between each elite
    /// M_s: Mutation multiplier
    /// D_sl: Shift limit
    /// Problem: The instance of the problem currently being solved
    /// </summary>
    public class BEGA {
        public Genotype[] P { get; set; }
        public Genotype[] S_G { get; set; }
        public Genotype[] S_E { get; set; }
        public int P_e { get; set; }
        public double M_d { get; set; }
        public double D_sl { get; set; }
        public double M_s { get; set; }
        public Instance Problem;
        public int N { get; set; }

        public BEGA(){}
        public BEGA(Instance instance, int p = 90, int pE = 15, double eliteDistance = 3, double mutationMultipler = 5,
            double shiftLimit = 0.075) {
            this.P_e = pE;
            this.M_d = eliteDistance;
            this.M_s = mutationMultipler;
            this.D_sl = shiftLimit;
            this.Problem = instance;
            this.N = instance.Nodes.Length;
            
            // Initialise the population and sort
            this.P = this.InitialisePopulation(instance, p);
            this.P = this.P.OrderBy(c => c.Cost).ToArray();
        }

        public virtual void Evolve() {
            // calc SGM
            var sgm = BuildSgm(this.P);
            
            // Evolve the two subpopulation
            var intensificationPop = this.Intensification(this.P, this.P_e, sgm);
            var diversityPop = this.Diversity(this.P, this.P.Length - intensificationPop.Length, sgm);
            
            // Combine the two populations
            Genotype[] newPop = new Genotype[this.P.Length];
            try {
                Array.Copy(intensificationPop, 0, newPop, 0, intensificationPop.Length);
            }
            catch (Exception e) {
                Console.Error.Write(e);
                System.Environment.Exit(1);
            }
            try {
                Array.Copy(diversityPop, 0, newPop, intensificationPop.Length, diversityPop.Length);
            }
            catch (Exception e) {
                Console.Error.Write(e);
                System.Environment.Exit(1);
            }

            // Sort
            newPop = newPop.OrderBy(c => c.Cost).ToArray();
            this.P = newPop;
        }

        /// <summary>
        /// Perturbation matrices are used to produce individuals for the diversity and intensification methods. This
        /// method allows for the both the negative and positive perturbations to be made.
        ///
        /// According to Alg 3 and Alg 4 from Zhang et al. (2018)
        /// </summary>
        /// <param name="sgm">Similarity guide matrix (S_G)</param>
        /// <param name="isNegativePerturbation">Negative perturbation for diversity method. Set to false to calculate
        /// the positive perturbation for the intensification method</param>
        /// <returns></returns>
        public double[,] BuildPerturbationMatrix(double[,] sgm, double CA, bool isNegativePerturbation = true) {
            Random rand = new Random(Guid.NewGuid().GetHashCode());
            
            double[,] M_N = new double[sgm.GetLength(0),sgm.GetLength(1)];

            for (int i = 0; i < 1; i++) {
                M_N[i, 0] = -1;
            }
            
            
            // for (int j = 1; j < M_N.GetLength(1); j++) {
            Parallel.For(1, M_N.GetLength(1), j => {
                M_N[j, j] = -1;
                if (rand.NextDouble() < CA) {
                    int k = -1;
                    if (isNegativePerturbation) {
                        k = GetMaxIndex(sgm, j);
                    }
                    else {
                        k = GetMinIndex(sgm, j);
                    }

                    int l = rand.Next(1, sgm.GetLength(0));
                    while (k == l || l == j) {
                        l = rand.Next(1, sgm.GetLength(0));
                    }

                    double maxSvj = sgm[k, j];
                    double rng = rand.NextDouble();
                    double r_n = CA * rng * maxSvj;

                    double originalValue, newValue;
                    originalValue = sgm[k, j];
                    newValue = sgm[k, j] - r_n;

                    M_N[k, j] = sgm[k, j] - r_n;
                    M_N[l, j] = sgm[l, j] + r_n;
                }
                else {
                    for (int i = 0; i < sgm.GetLength(0); i++) {
                        M_N[i, j] = sgm[i, j];
                    }
                }
            });
            //}

            return M_N;
        }


        public virtual int[] NewTemporaryIndividual(double[,] matrix, int n) {
            int[] tempSeq = new int[n];
            for (int i = 0; i < tempSeq.Length; i++) {
                tempSeq[i] = GetNodeValue(matrix, i);
            }

            return GeneticAlgorithm.GeneticAlgorithmFunctions.FixGenotype(tempSeq);
        }
        public int GetMinIndex(double[,] sgm, int j) {
            double minValue = 2;
            int minIndex = -1;
            double lastVal = -1;
            List<int> zeroIndexes = new List<int>();
            for (int i = 1; i < sgm.GetLength(0); i++) {
                // Skip invalid values
                if (sgm[i, j] < 0) {
                    continue;
                }

                if (sgm[i, j] == 0) {
                    zeroIndexes.Add(i);
                    continue;
                }

                lastVal = sgm[i, j];
                if (sgm[i, j] < minValue) {
                    minValue = sgm[i, j];
                    minIndex = i;
                }
            }

            if (minIndex < 0 && zeroIndexes.Count > 0) {
                Random rand = new Random(Guid.NewGuid().GetHashCode());
                minIndex = zeroIndexes[rand.Next(0, zeroIndexes.Count)];
            }
            if (minIndex < 1) {
                Console.Error.Write(sgm.ToString());
                throw new Exception("BEGA.GetMaxIndex(): Invalid SGM");
            }
            return minIndex;
        }
        public int GetMaxIndex(double[,] sgm, int j) {
            double maxValue = -1;
            int maxIndex = -1;
            for (int i = 1; i < sgm.GetLength(0); i++) {
                if (sgm[i, j] > maxValue) {
                    maxValue = sgm[i, j];
                    maxIndex = i;
                }
            }

            if (maxIndex < 1) {
                Console.Error.Write(sgm.ToString());
                throw new Exception("BEGA.GetMaxIndex(): Invalid SGM");
            }
            return maxIndex;
        }
        
        
        
        /// <summary>
        /// Diversity method according to Zhang et al. (2018) in Algorithm 3.
        /// Some parts have been tweaked to adapt to ordered problem constraints
        /// </summary>
        /// <param name="pop"></param>
        public virtual Genotype[] Diversity(Genotype[] pop, int subPopSize, double[,] sgm) {
            Random rand = new Random(Guid.NewGuid().GetHashCode());
            Genotype[] newPop = new Genotype[pop.Length];
            
            
            // Get the control amplitude
            double CA = CalcControlAmplitude(pop);
            
            //----------------------------------------------------------------------------------------------------------
            // Step 1: According to S_G, compute each vector of np_SGM M_N from Alg. 3 (Zhang et al. 2018)
            //----------------------------------------------------------------------------------------------------------
            // var sgm = BuildSgm(pop);
            var M_N = BuildPerturbationMatrix(sgm, CA, true);
            
            
            // Create the temporary individuals that reflect the subpopulation
            int[][] temporaryIndividuals = new int[pop.Length][];
            Parallel.For(0, newPop.Length, ind => {
                temporaryIndividuals[ind] = this.NewTemporaryIndividual(M_N, this.N);
            });
            // Calculate the mutation probability
            double p_m = (1.0 + this.M_s * CA) / (double) N;
            
            //----------------------------------------------------------------------------------------------------------
            // Step 2: Breed temporary individuals by using np-SGM M_N and Algorithm 1, and execute crossover
            //----------------------------------------------------------------------------------------------------------
            Parallel.For(0, newPop.Length, ind => {
                int[] O = new int[this.N];
                for (int i = 0; i < O.Length; i++) {
                    O[i] = rand.NextDouble() <= 0.5 ? temporaryIndividuals[ind][i] : pop[ind].Sequence[i];
                }

                // Step 3: Execute mutation according to mutation probability pm in Eq. (17) and ms is multiplier factor
                O = GeneticAlgorithm.GeneticAlgorithmFunctions.Mutate(O, p_m);
                int cost = Problem.CalcCostRoute(O);
                newPop[ind] = new Genotype(O, cost);
            });
            
            // Sort and return
            newPop = newPop.OrderBy(c => c.Cost).ToArray();
            Genotype[] returnPop = new Genotype[subPopSize];
            Array.Copy(newPop, 0, returnPop, 0, returnPop.Length);
            return returnPop;
        }

        public Genotype[] GetElite(int elite, Genotype[] pop, double minDist) {
            List<Genotype> elites = new List<Genotype>();
            Genotype lastElite = null;
            int pIndex = 0;
            while (elites.Count < elite) {
                // Get the first elite
                if (elites.Count == 0) {
                    lastElite = pop[pIndex];
                    elites.Add(lastElite);
                    pIndex++;
                    continue;
                }
                
                // Get the rest of the elites
                if (GeneticAlgorithm.GeneticAlgorithmFunctions.HammingDis(lastElite.Sequence, pop[pIndex].Sequence) >= minDist) {
                    lastElite = pop[pIndex];
                    elites.Add(lastElite);
                }
                pIndex++;
                if (pIndex >= pop.Length) {
                    break;
                }
            }
            return elites.ToArray();
        }
        
        
        /// <summary>
        /// Intensification method according to Zhang et al. (2018) in Algorithm 3.
        /// Some parts have been tweaked to adapt to ordered problem constraints
        /// </summary>
        /// <param name="pop"></param>
        /// <param name="subPopSize"></param>
        /// <returns></returns>
        public virtual Genotype[] Intensification(Genotype[] pop, int subPopSize, double[,] sgm) {
            Random rand = new Random(Guid.NewGuid().GetHashCode());
            
            //----------------------------------------------------------------------------------------------------------
            // Step 1: Compute S_E and copy the best individuals into the next generation
            //----------------------------------------------------------------------------------------------------------
            Genotype[] subPop = GetElite(subPopSize, pop, M_d);
            Genotype[] newPop = new Genotype[subPop.Length];
            // var sgm = BuildSgm(subPop);
            var S_E = GetReferencePoint(sgm);

            //----------------------------------------------------------------------------------------------------------
            // Step 2: According to S_E, compute each vector of ppSGM M_P
            //----------------------------------------------------------------------------------------------------------
            // Get the control amplitude
            double CA = CalcControlAmplitude(pop);
            var M_P = BuildPerturbationMatrix(sgm, CA, false);

            //----------------------------------------------------------------------------------------------------------
            // Step 3: Breed temporary individuals by using pp=SGM and M_P and Alg 1. Similar to Alg 3, execute
            // crossover between temporary individuals and intensification subpopulation after updating. 
            //----------------------------------------------------------------------------------------------------------
            // Calc probabilty of mutation
            double p_m = 1 / (double)this.N;
            
            // Create the temporary individuals that reflect the subpopulation
            int[][] temporaryIndividuals = new int[pop.Length][];
            Parallel.For(0, subPop.Length, ind => {
                temporaryIndividuals[ind] = this.NewTemporaryIndividual(M_P, this.N);
            });
            
            // Crossover and mutation
            Parallel.For(0, subPop.Length, i => {
                if (i == 0) {
                    newPop[i] = subPop[i];
                }
                var candidate = temporaryIndividuals[i];
                var genotype = subPop[i];
                var seq = GeneticAlgorithm.GeneticAlgorithmFunctions.Crossover(candidate, genotype.Sequence);
                seq = GeneticAlgorithm.GeneticAlgorithmFunctions.Mutate(seq, p_m);
                var cost = Problem.CalcCostRoute(seq);
                newPop[i] = new Genotype(seq, cost);
            });
            
            // Sort and return
            newPop = subPop.OrderBy(c => c.Cost).ToArray();
            return newPop;

        }
        
        public static int GetNodeValue(double[,] matrix, int i) {
            i += 1;
            Random rand = new Random(Guid.NewGuid().GetHashCode());
            double sum = 0.0;
            for (int j = 1; j < matrix.GetLength(1); j++) {
                double val = rand.NextDouble();
                try {
                    if (val > sum && val < sum + matrix[i, j]) {
                        return j;
                    }

                    sum += matrix[i, j];
                }
                catch (Exception e) {
                    Console.Error.WriteLine(e);
                    System.Environment.Exit(1);
                }
            }

            return matrix.GetLength(1) - 1;
        }

        public double CalcControlAmplitude(Genotype[] pop) {
            double[,] sgm = BuildSgm(pop);
            int[] X_cp = GetReferencePoint(sgm);
            double D_l = LinearDiversityIndex(pop, X_cp);
            return D_l > this.D_sl ? D_l : this.D_sl;
        }
        
        public Genotype[] InitialisePopulation(Instance instance, int p) {
            Random rand = new Random(Guid.NewGuid().GetHashCode());
            Genotype[] pop = new Genotype[p];
            for (int i = 0; i < p; i++) {
                // Build a sequence
                int[] sequence = new int[this.Problem.Nodes.Length];
                for (int j = 1; j <= sequence.Length; j++) {
                    sequence[j - 1] = j;
                }
                
                // Shuffle it and get the cost
                sequence = sequence.OrderBy(x => rand.Next()).ToArray();
                int cost = instance.CalcCostRoute(sequence);
                
                Genotype genotype = new Genotype(sequence, cost);
                pop[i] = genotype;
            }
            return pop;
        }


        public virtual double[,] BuildSgm(Genotype[] pop) {
            int N = pop[0].Sequence.Length;
            double[,] sgm = new double[N + 1, N + 1];
            int[,] countMatrix = new int[N+1, N + 1];

            // Build the count table
            Parallel.For(0, countMatrix.GetLength(0), r => {
                for (int c = 0; c < countMatrix.GetLength(1); c++) {
                    try {
                        countMatrix[r, c] = 0;
                    }
                    catch (Exception e) {
                        Console.Error.WriteLine(e.ToString());
                    }
                }
            });

            for (int i = 0; i < pop.Length; i++) {
                int[] seq = pop[i].Sequence;
                for (int j = 0; j < seq.Length; j++) {
                    try {
                        countMatrix[seq[j], j + 1]++;
                    }
                    catch (Exception e) {
                        Console.Error.WriteLine(e.ToString());
                    }
                }
            }
            
            // Calculate the distribution
            for (int j = 0; j < countMatrix.GetLength(0); j++) {
                for (int i = 0; i < countMatrix.GetLength(1); i++) {
                    sgm[j, i] = countMatrix[j, i] / (double) pop.Length;
                }
            }

            return sgm;
        }

        public virtual int[] GetReferencePoint(double[,] sgm) {
            int[] seq = new int[sgm.GetLength(0)];

            Parallel.For(1, sgm.GetLength(1), c => {
                int maxN = -1;
                double maxValue = -1;

                for (int r = 1; r < sgm.GetLength(0); r++) {
                    if (sgm[r, c] > maxValue) {
                        maxValue = sgm[r, c];
                        maxN = r;
                    }
                }

                seq[c] = maxN;
            });

            int[] results = new int[seq.Length - 1];
            Array.Copy(seq,1, results, 0, seq.Length-1);
            return results;
        }

        public static double LinearDiversityIndex(Genotype[] pop, int[] reference) {
            double result = 0.0;
            double P = pop.Length;
            Parallel.For(0, pop.Length, i => {
                var t = pop[i];
                result += DistanceP(t.Sequence, reference);
            });
            return result / P;
        }

        public static double DistanceP(int[] x, int[] reference) {
            double N = x.Length;
            return GeneticAlgorithm.GeneticAlgorithmFunctions.HammingDis(x, reference) / N;
        }

    }
}