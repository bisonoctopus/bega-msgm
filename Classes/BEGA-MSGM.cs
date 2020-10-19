using System;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Threading.Tasks;
using GeneticAlgorithm;
using Markov;

namespace BEGA.Classes {
    public class BEGA_MSGM : BEGA {
        
        public BEGA_MSGM(Instance instance, int p = 90, int pE = 15, double eliteDistance = 5, double mutationMultipler = 5,
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
            
            Console.WriteLine($"[INFO] BEGA MSGM initiated | P {p} | E: {pE}");
        }
        
        /// <summary>
        /// Diversity method according to Zhang et al. (2018) in Algorithm 3.
        /// Some parts have been tweaked to adapt to ordered problem constraints
        /// </summary>
        /// <param name="pop"></param>
        public override Genotype[] Diversity(Genotype[] pop, int subPopSize, double[,] sgm) {
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
        
        
        /// <summary>
        /// Intensification method according to Zhang et al. (2018) in Algorithm 3.
        /// Some parts have been tweaked to adapt to ordered problem constraints
        /// </summary>
        /// <param name="pop"></param>
        /// <param name="subPopSize"></param>
        /// <returns></returns>
        public override Genotype[] Intensification(Genotype[] pop, int subPopSize, double[,] sgm) {
            Random rand = new Random(Guid.NewGuid().GetHashCode());
            
            //----------------------------------------------------------------------------------------------------------
            // Step 1: Compute S_E and copy the best individuals into the next generation
            //----------------------------------------------------------------------------------------------------------
            Genotype[] subPop = GetElite(subPopSize, pop, M_d);
            Genotype[] newPop = new Genotype[subPop.Length];
            var sgm = BuildSgm(subPop);
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
                // var candidate = temporaryIndividuals[i];
                var candidate = pop[GeneticAlgorithmFunctions.Select(pop.Length)].Sequence;
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
        
        public override void Evolve() {
            // calc SGM
            var sgm = BuildFitnessSgm(this.P);
            
            // Evolve the two subpopulation
            var intensificationPop = this.Intensification(this.P, this.P_e, sgm);
            var diversityPop = this.Diversity(this.P, this.P.Length - intensificationPop.Length, sgm);
            
            // Combine the two populations
            Genotype[] newPop = new Genotype[this.P.Length];
            for (int i = 0; i < intensificationPop.Length; i++) {
                newPop[i] = intensificationPop[i];
            }
            for (int i = 0; i < diversityPop.Length; i++) {
                newPop[i + intensificationPop.Length] = diversityPop[i];
            }
            
            
            // Sort and go to the next generation
            newPop = newPop.OrderBy(c => c.Cost).ToArray();
            if (this.P[0].Cost < newPop[0].Cost) {
                newPop[0] = this.P[0];
            }
            this.P = newPop;
        }
        
        

        public override int[] GetReferencePoint(double[,] sgm) {
            return this.NewTemporaryIndividual(sgm, sgm.GetLength(0) - 1);
        }

        public override int[] NewTemporaryIndividual(double[,] sgm, int n) {
            Random rand = new Random(Guid.NewGuid().GetHashCode());

            int[] seq = new int[n];
            var msgm = GeneticAlgorithmFunctions.CopyArray(sgm);
            seq[0] = 1;
            
            // Go through sgm and build new sequence
            for (int i = 1; i < seq.Length; i++) {
                int from = seq[i - 1];

                // Remove the possibility of going back to the FROM node
                var next = Markov.Markov.GetNext(from, msgm);
                seq[i] = next;
                Markov.Markov.RemoveColumn(from, msgm);
            }


            if (seq.Distinct().Count() != seq.Length || !GeneticAlgorithmFunctions.IsDistinct(seq)) {
                Console.Error.WriteLine("Error with sequence");
                Console.WriteLine(String.Join(", ", seq));
                System.Environment.Exit(1);
            }
            return seq;
        }


        public override double[,] BuildSgm(Genotype[] pop) {
            int[][] sequences = new int[pop.Length][];
            for (int i = 0; i < pop.Length; i++) {
                sequences[i] = pop[i].Sequence;
            }
            return Markov.Markov.TransitionMatrix(sequences);
        }

        public double[,] BuildFitnessSgm(Genotype[] pop) {
            double[,] fitnessMatrix = new double[this.N + 1, this.N + 1];
            
            // Set invalid cells to -1
            for (int i = 0; i < this.N + 1; i++) {
                for (int j = 0; j < this.N + 1; j++) {
                    if(i == 0 || j == 0 || i == j)
                    {
                        fitnessMatrix[i, j] = -1;
                    } else
                    {
                        fitnessMatrix[i, j] = 0;
                    }
                }
            }

            double totalFitness = 0;
            foreach (var ind in pop) {
                totalFitness += ind.Cost;
            }
            
            // Sum the transition fitness values
            foreach (var ind in pop) {
                var from = -1;
                var to = -1;
                var relativeFitness = ind.Cost / totalFitness;

                for (int i = 0; i < ind.Sequence.Length - 1; i++) {
                    from = ind.Sequence[i];
                    to = ind.Sequence[i + 1];
                    if (from > 0 && to > 0) {
                        fitnessMatrix[from, to] += relativeFitness;
                    }
                }
                // Get the loop back
                fitnessMatrix[to, ind.Sequence[0]] += relativeFitness;
            }

            return fitnessMatrix;
        }

        public double[,] BuildBalancedSgm(Genotype[] pop) {
            double p = pop.Length;

            double[,] fitnessMatrix = BuildFitnessSgm(pop);
            double[,] transitionMatrix = BuildSgm(pop);
            double[,] balancedMatrix = new double[N + 1, N + 1];
            
            // Set invalid cells to -1
            for (int i = 0; i < this.N + 1; i++) {
                for (int j = 0; j < this.N + 1; j++) {
                    if(i == 0 || j == 0 || i == j)
                    {
                        balancedMatrix[i, j] = -1;
                    } else
                    {
                        balancedMatrix[i, j] = 0;
                    }
                }
            }
            // Calculate the balanced matrix
            for (var i = 1; i < balancedMatrix.GetLength(0); i++) {
                for (var j = 1; j < balancedMatrix.GetLength(1); j++) {
                    if (fitnessMatrix[i, j] < 0 || transitionMatrix[i, j] < 0) {
                        balancedMatrix[i, j] = -1;
                    }
                    else {
                        balancedMatrix[i, j] = (fitnessMatrix[i, j] + transitionMatrix[i, j]) / (double) 2;
                    }
                }
            }

            return balancedMatrix;
        }
    }
}