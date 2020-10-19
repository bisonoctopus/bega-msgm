using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;

namespace BEGA.Classes {
    public struct Node {
        public bool IsDepot;
        public int Id;
        public double X;
        public double Y;
        public double Demand;

        public Node(int id, double x, double y, double demand = -1, bool isDepot = false) {
            this.Id = id;
            this.X = x;
            this.Y = y;
            this.IsDepot = isDepot;
            this.Demand = demand;
        }

        public void SetDemand(double demand) {
            this.Demand = demand;
        }
    }
    
    public class Instance {
        public Node[] Nodes { get; private set; }
        public int[] Solution { get; private set; }
        public int[][] Routes { get; private set; }
        public int Optimum { get; private set; }
        public int[,] CostMatrix { get; private set; }
        public string Name { get; private set; }
        public string Comment { get; private set; }
        public string Type { get; private set; }
        public int Dimension { get; private set; }
        public double Capacity { get; private set; }
        public int K { get; private set; }
        public bool IsComplete { get; private set; }
        
        
        public Instance(string path, string instanceName, int offset, int maxDimensions = 1000) {
            // Get the problem data
            this.ReadFile($"{path}\\{instanceName}.tsp", maxDimensions);
            if (this.Dimension > maxDimensions) {
                this.IsComplete = false;
                return;
            }
                
            
            this.CalcCostMatrix(offset);
            
            // Get the solution data
            this.ReadSol($"{path}\\{instanceName}.opt.tour");
            this.Optimum = this.CalcCostRoute(this.Solution);
            this.IsComplete = true;
        }

        public int CalcCostRoute(int[] route) {
            if (route.Length == 0) {
                return -1;
            }
            int cost = 0;
            int fromNode = 0;
            int toNode = 0;
            int dist = 0;
            List<int> costs = new List<int>();
            for (int i = 0; i < route.Length - 1; i++) {
                fromNode = route[i];
                toNode = route[i + 1];
                dist = this.CostMatrix[fromNode, toNode];
                cost += dist;
                costs.Add(dist);
            }
            
            // Go back to base
            dist = this.CostMatrix[toNode, route[0]];
            cost += dist;
            costs.Add(dist);
            return cost;
        }
        
        public int CalcCostRoutes(List<List<int>> routes) {
            
            int totalCost = 0;
            foreach (var route in routes) {
                if (route.Count == 0) {
                    continue;
                }

                // Route is only 1 delivery, go to node and then back to depot
                if (route.Count == 1) {
                    totalCost += this.CostMatrix[0, route[0]] * 2;
                    continue;
                }
                int routeCost = 0;
                int fromNode = -1;
                int toNode = -1;
                for (int i = 0; i < route.Count - 1; i++) {
                    fromNode = route[i];
                    toNode = route[i + 1];
                    // Cost from depot
                    if (i == 0) {
                        routeCost += this.CostMatrix[0, fromNode];
                    }
                    
                    // Cost to next node
                    routeCost += this.CostMatrix[fromNode, toNode];
                }

                // Cost to depot
                routeCost += this.CostMatrix[toNode, 0];
                totalCost += routeCost;
            }

            return totalCost;
        }
        public int CalcCostRoutes(int[][] routes) {
            int totalCost = 0;
            foreach (var route in routes) {
                if (route.Length == 0) {
                    continue;
                }

                // Route is only 1 delivery, go to node and then back to depot
                if (route.Length == 1) {
                    totalCost += this.CostMatrix[0, route[0]] * 2;
                    continue;
                }
                int routeCost = 0;
                int fromNode = -1;
                int toNode = -1;
                for (int i = 0; i < route.Length - 1; i++) {
                    fromNode = route[i];
                    toNode = route[i + 1];
                    // Cost from depot
                    if (i == 0) {
                        routeCost += this.CostMatrix[0, fromNode];
                    }
                    
                    // Cost to next node
                    routeCost += this.CostMatrix[fromNode, toNode];
                }

                // Cost to depot
                routeCost += this.CostMatrix[toNode, 0];
                totalCost += routeCost;
            }

            return totalCost;
        }
        
        public int CalcDistance(Node a, Node b) {
            double x = a.X - b.X;
            double y = a.Y - b.Y;
            return (int)Math.Round(Math.Sqrt(Math.Pow(x, 2) + Math.Pow(y, 2)));
        }
        public void CalcCostMatrix(int offset = 0) {
            int[,] costMatrix = new int[Nodes.Length + offset,Nodes.Length + offset];

            for (int i = 0; i < costMatrix.GetLength(0); i++) {
                for (int j = 0; j < costMatrix.GetLength(1); j++) {
                    costMatrix[i, j] = -1;
                }
            }
            Parallel.For(0, this.Nodes.Length, i => {
                for (int j = 0; j < Nodes.Length; j++) {
                    // if (i == j || i < offset || j < offset) {
                    //     costMatrix[i, j] = -1;
                    //     continue;
                    // }
                    var fromNode = this.Nodes[i];
                    var toNode = this.Nodes[j];
                    
                    if (fromNode.Id == toNode.Id) {
                        costMatrix[fromNode.Id, toNode.Id] = -1;
                        continue;
                    }
                    int dist = this.CalcDistance(fromNode, toNode);
                    costMatrix[fromNode.Id, toNode.Id] = dist;
                    //costMatrix[i, j] = this.CalcDistance(this.Nodes[i], this.Nodes[j]);
                }
            });

            this.CostMatrix = costMatrix;
        }

        public string[] CleanRow(string line) {
            
            // Try splitting by tabs
            string[] row = line.Trim().Split('\t');
                
            // Try splitting by colon
            if (row.Length == 1) {
                row = line.Trim().Split(":");
            }
            // Finally try splitting by space
            if (row.Length == 1) {
                row = line.Trim().Split(" ");
            }
                
            // Get rid of empty elements
            List<string> data = new List<string>();
            foreach (var item in row) {
                if (item.Length > 0) {
                    data.Add(item);
                }
            }
            return data.ToArray();
        }
        public void ReadFile(string path, int maxDimension) {
            // Console.WriteLine(path);
            string[] lines = new string[0];
            try {
                lines = File.ReadAllLines(path);
            } catch (Exception e) {
                throw e;
            }

            bool isNodeCoordSection = false;
            bool isDemandSection = false;
            bool isDepotSection = false;
            
            List<Node> nodesList = new List<Node>();
            
            // Get the meta data first
            foreach (var line in lines) {
                string[] row = CleanRow(line);
                if (line.Length == 0) {
                    continue;
                }
                
                // Trim the elements of the row and change text to lower case
                for (int i = 0; i < row.Length; i++) {
                    row[i] = row[i].Trim();
                }
                string key = row[0].Trim().ToLower();
                
                // Make sure we're in the right sections
                if (key == "node_coord_section" || key == "display_data_section") {
                    isNodeCoordSection = true;
                    isDemandSection = false;
                    isDepotSection = false;
                    continue;
                }

                if (key == "demand_section") {
                    isNodeCoordSection = false;
                    isDemandSection = true;
                    isDepotSection = false;
                    continue;
                }

                if (key == "depot_section") {
                    isNodeCoordSection = false;
                    isDemandSection = false;
                    isDepotSection = true;
                    continue;
                }

                if (key == "eof") {
                    break;
                }
                
                // Get the instance meta data
                if (!isNodeCoordSection && !isDemandSection && !isDepotSection) {
                    string value = "";
                    if (row.Length > 1) {
                        value = row[1].Trim();
                    }
                    else {
                        value = row[0].Trim();
                    }

                    key = key.Replace(":","");
                    key = key.Trim();

                    if (key == "name") {
                        this.Name = value;
                        if (value.Contains("-k")) {
                            string[] k = value.Split(new string[] {"-k"}, StringSplitOptions.None);
                            this.K = Int32.Parse(k[1]);
                        }
                        else {
                            this.K = -1;
                        }
                    }

                    if (key == "comment") {
                        this.Comment = value;
                    }

                    if (key == "type") {
                        this.Type = value;
                    }

                    if (key == "dimension") {
                        this.Dimension = Int32.Parse(value);
                        if (this.Dimension > maxDimension) {
                            return;
                        }
                    }

                    if (key == "capacity") {
                        this.Capacity = Double.Parse(value);
                    }
                }
            }
            
            
            // Second iteration to get the node data
            Node[] nodesArr = new Node[this.Dimension];
            
            isNodeCoordSection = false;
            isDemandSection = false;
            isDepotSection = false;
            
            foreach (var line in lines) {
                string[] row = CleanRow(line);
                string key = row[0].Trim().ToLower();
                
                // Make sure we're in the right sections
                if (key == "node_coord_section") {
                    isNodeCoordSection = true;
                    isDemandSection = false;
                    isDepotSection = false;
                    continue;
                }
                
                if (key == "demand_section") {
                    isNodeCoordSection = false;
                    isDemandSection = true;
                    isDepotSection = false;
                    continue;
                }
                
                if (key == "depot_section") {
                    isNodeCoordSection = false;
                    isDemandSection = false;
                    isDepotSection = true;
                    continue;
                }

                if (key == "eof") {
                    break;
                }
                
                // Get the coord data
                if (isNodeCoordSection) {
                    int id = -1;
                    double x = double.NaN;
                    double y = double.NaN;
                    try {
                        id = Int32.Parse(key);
                        x = double.Parse(row[1]);
                        y = double.Parse(row[2]);
                    }
                    catch (Exception e) {
                        throw e;
                    }


                    Node temp = new Node(id, x, y);
                    //nodesList.Add(temp);
                    nodesArr[id - 1] = temp;
                }

                
                
                // Get the demand data
                if (isDemandSection) {
                    int id = Int32.Parse(key);
                    double demand = double.Parse(row[1]);
                    nodesArr[id - 1].Demand = demand;
                    int i = id - 1;
                    nodesArr[id - 1].SetDemand(demand);
                    // Console.WriteLine($"i: {i} | id: {id} | d: {nodesArr[id - 1].Demand}");
                }
                
                // Get the depot data
                if (isDepotSection) {
                    try {
                        int id = Int32.Parse(key);
                        if (id > 0) {
                            nodesArr[id - 1].IsDepot = true;
                        }
                    } catch(Exception e)
                    {
                        if (key == "") {
                            continue;
                        }
                        // Console.WriteLine(key);
                    }
                }
            }
            this.Nodes = nodesArr;
        }

        public void ReadSol(string path) {
            string[] lines = new string[0];
            try {
                lines = File.ReadAllLines(path);
            }
            catch (Exception e) {
                Console.WriteLine($"[INFO] No solution file for: {this.Name}");
                this.Optimum = 0;
            }

            List<int[]> routes = new List<int[]>();
            List<int> solution = new List<int>();
            bool isTourSection = false;
            
            foreach (var line in lines) {
                var row = this.CleanRow(line);
                
                // Skip all until we reach the tour section
                if (row[0].ToLower() == "tour_section") {
                    isTourSection = true;
                    continue;
                }
                if (!isTourSection) {
                    continue;
                }
                
                if (row[0] == "-1") {
                    // End of file
                    break;
                }
                if (row[0].ToLower().Contains("cost")) {
                    row = row[0].Split(' ');
                    this.Optimum = Int32.Parse(row[1]);
                    break;
                }
                
                // Create the routes vars
                List<int> routeList = new List<int>();
                string[] routeString = new string[0];

                if (row.Length == 1) {
                    string[] temp = row[0].Split(' ');
                    if (temp.Length == 1) {
                        routeString = temp;
                    }
                    else {
                        routeString = row[0].Split(' ');
                        routeString = routeString.Skip(1).Take(routeString.Length - 1).ToArray();
                    }
                }
                else {
                    routeString = row[1].Split(' ');    
                }
                
                for(int i = 0; i < routeString.Length; i++) {
                    try {
                        var node = Int32.Parse(routeString[i]);
                        routeList.Add(node);
                        solution.Add(node);
                    }
                    catch (Exception e) {
                        if (routeString[i] == "") {
                            continue;
                        }
                        Console.WriteLine(e.ToString());
                        Console.WriteLine(routeString[i]);
                    }
                }

                if (routeList.Count > 0) {
                    routes.Add(routeList.ToArray());
                }
            }

            this.Routes = routes.ToArray();
            this.Solution = solution.ToArray();
        }
    }
}