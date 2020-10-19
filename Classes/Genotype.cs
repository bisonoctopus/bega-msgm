using System.Collections.Generic;

namespace BEGA.Classes {
    public class Genotype {
        
        public int[] Sequence;
        public int Cost;
        public double Fitness;
        public List<int>[] Routes;

        public Genotype(int[] sequence, int cost = -1) {
            this.Sequence = sequence;
            this.SetCost(cost);
            this.Routes = new List<int>[0];
        }

        public void SetGenotype(int[] seq, int cost = -1) {
            this.Sequence = seq;
            this.SetCost(cost);
        }

        public void SetCost(int cost) {
            this.Cost = cost;
            this.Fitness = 1 / (double) cost;
        }
    }
}