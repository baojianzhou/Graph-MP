package edu.albany.cs.tailApprox;

import edu.albany.cs.base.DisjointSet;
import edu.albany.cs.fastPCST.FastPCST;
import edu.albany.cs.headApprox.Digraph;
import edu.albany.cs.headApprox.DirectedEulerianCycle;
import org.apache.commons.lang3.ArrayUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

/**
 * Tail Approximation
 * <p>
 * Aglorithm 1 of
 * "http://people.csail.mit.edu/ludwigs/papers/icml15_graphsparsity.pdf" Title :
 * "A Nearly-Linear Time Framework for Graph-Structured Sparsity" Authors :
 * Ludwig Schmidt, Chinmay Hegde, and Piotr Indyk
 *
 * @author Baojian bzhou6@albany.edu
 */
public class TailApprox {

    /** Graph G is split into 3 parts edges, edgeCosts, and prizesPi. */
    /**
     * represents graph G
     */
    private final ArrayList<Integer[]> edges;
    /**
     * edge costs c(e) (in algorithm 1 c(e) = w(e) + B / s).
     */
    private final ArrayList<Double> edgeCostsc;
    /**
     * node prizes pi
     */
    private final ArrayList<Double> prizesPi;
    /**
     * the budget
     */
    private final double C;
    /**
     * total sparsity
     */
    private final int s;
    /**
     * the constant parameter, which is > 2.
     */
    private final double nu;
    /**
     * delta = min(1/2,1/nu)
     */
    private final double delta;
    /**
     * g (number of active clusters in PCST's forest)
     */
    private final int g;

    /**
     * additional variable. may not be useful.
     */
    private int[] trueSubGraph;
    private int verboseLevel = 0;

    /**
     * the result forest that is returned by Tail algorithm.
     */
    public F bestForest;
    /** This parameter "valid" to show the result is valid. */
    /**
     * But, it may not be useful.
     */
    public boolean valid;

    /**
     * A general input constructor
     *
     * @param edges       the edges of graph G
     * @param edgeCostsc  the edge costs corresponding to edges in G
     * @param prizesPi    the node prizes Pi
     * @param g           the number of connected components
     * @param costBudgetC the cost budget
     * @param nu          the constant parameter nu, which is > 2.
     * @param delta       the constant parameter delta, which is min(1/2,1/\nu).
     */
    public TailApprox(ArrayList<Integer[]> edges, ArrayList<Double> edgeCostsc, ArrayList<Double> prizesPi, int g,
                      double costBudgetC, double nu, double delta) {

        this.edges = edges;
        this.edgeCostsc = edgeCostsc;
        this.prizesPi = prizesPi;
        this.g = g;
        this.s = -1;
        this.C = costBudgetC;
        this.nu = nu;
        this.delta = delta;

        this.bestForest = run();

        /**
         * check equation 8 : the parameter of equation 8 is cT > 1 ; cT is
         * arbitrary and fixed constant.
         */
        double cT = Math.sqrt(1.0D + 3.0D / (this.nu - 2.0D)); //
        double[] b = new double[prizesPi.size()];
        for (int i = 0; i < prizesPi.size(); i++) {
            b[i] = prizesPi.get(i);
        }
        double[] bS = new double[b.length];
        double[] bSPrime = new double[b.length];
        for (int i = 0; i < b.length; i++) {
            if (bestForest.nodesInF.contains(i)) {
                bS[i] = b[i];
            } else {
                bS[i] = 0.0D;
            }
            if (ArrayUtils.contains(this.trueSubGraph, i)) {
                bSPrime[i] = b[i];
            } else {
                bSPrime[i] = 0.0D;
            }
        }

        if (this.checkEqu8Valid(cT, b, bS, bSPrime)) {
            this.valid = true;
            System.out.println("result is valid");
        } else {
            this.valid = false;
            System.out.println("result of Tail approximation is invalid");
            System.exit(0);
        }
    }

    /**
     * This constructor is for Algorithm 1 only.
     *
     * @param edges        the edges corresponding to WGM
     * @param edgeCostsW:  the edge costs corresponding to WGM
     * @param z            the real vector that we want to estimate.
     * @param s            the sparsity parameter.
     * @param g            the number of connected components.
     * @param B            the budget constraint, it is a real number.
     * @param trueSubGraph let it be null if you do not need it.
     */
    public TailApprox(ArrayList<Integer[]> edges, ArrayList<Double> edgeCostsW, double[] z, int s, int g, double B,
                      int[] trueSubGraph) {
        this.trueSubGraph = trueSubGraph;
        /** edges */
        this.edges = edges;
        /** edge costs c, let c(e) = w(e) + (B / s) */
        this.edgeCostsc = new ArrayList<Double>();
        for (double w : edgeCostsW) {
            edgeCostsc.add(w + B / (s * 1.0D));
        }
        /** prize pi. */
        this.prizesPi = new ArrayList<Double>();
        if (z == null || z.length == 0) {
            System.out.println("Error: vector z is null ...");
            System.exit(0);
        }
        for (int i = 0; i < z.length; i++) {
            /** the prizes of pi equals to <z,z> */
            prizesPi.add(z[i] * z[i]);
        }
        this.s = s;
        this.g = g;
        /** C = 2.0*B in algorithm 1 */
        this.C = 2.0D * B;
        /** we set nu as constant value */
        this.nu = 2.5D;
        /** delta = min(1/2,1/v) */
        this.delta = Math.min(1 / 2.0D, 1 / this.nu);
        /** to check the parameters are all correct. */
        if (!isParametersValid()) {
            System.out.println("some parameters are not valid.");
            System.exit(0);
        }
        bestForest = run();
        if (verboseLevel >= 1) {
            System.out.println("prize : " + Arrays.toString(z));
            System.out.println("g: " + g);
            System.out.println("s: " + s);
            System.out.println("nu: " + nu);
            System.out.println("C: " + C);
            System.out.println(" delta : " + delta);
        }
        /** cT of equation 8 is > 1 and constant. */
        double cT = Math.sqrt(1.0D + 3.0D / (this.nu - 2.0D));
        double[] b = new double[z.length];
        for (int i = 0; i < z.length; i++) {
            b[i] = z[i];
        }
        double[] bS = new double[b.length];
        double[] bSPrime = new double[b.length];
        for (int i = 0; i < b.length; i++) {
            if (this.bestForest.nodesInF.contains(i)) {
                bS[i] = b[i];
            } else {
                bS[i] = 0.0D;
            }
            if (ArrayUtils.contains(this.trueSubGraph, i)) {
                bSPrime[i] = b[i];
            } else {
                bSPrime[i] = 0.0D;
            }
        }
        /** check equation 8 */
        if (this.checkEqu8Valid(cT, b, bS, bSPrime)) {
            this.valid = true;
            if (verboseLevel > 0) {
                System.out.println("Tail approximation is valid");
            }
        } else {
            this.valid = false;
            String errorMes = "result of Tail approximation is invalid.";
            System.out.println(errorMes);
            System.exit(0);
        }
    }

    private boolean isParametersValid() {

        /** to check the prize is valid */
        for (double p : prizesPi) {
            if (p < 0.0D) {
                String errorMes = "the prize should not be negative.";
                System.out.println(errorMes);
                return false;
            }
        }
        /** to check the graph is a connected graph. */
        if (!isConnected(edges)) {
            String errorMes = "the graph is not connected.";
            System.out.println(errorMes);
            return false;
        }
        return true;
    }

    /**
     * This is the tail algorithm part.
     *
     * @return F a forest is returned by Tail algorithm
     */
    private F run() {

        /** the minimum positive prize entry in pi vector. */
        double minPi = getMinPi();
        /** lambda is trade off parameter */
        double lambda0 = minPi / (2.0D * C);
        /** to get a forest F using the fast pcsf algorithm. */
        F forest = PCSF_GW(edges, getCostsLambda(lambda0), prizesPi, g);
        /** the special case satisfied, then return the forest. */
        if ((forest.costF <= 2.0D * C) && (forest.piFBar <= 0.0D)) {
            return forest;
        }
        /** setting these three parameters \lambda_r, \lambda_l, and epsilon */
        double lambdaR = 0.0D;
        double lambdaL = 3.0D * getSumPrizePi();
        double lambdaUnit = lambdaL;
        double epsilon = (minPi * delta) / C;
        /** binary search to find a forest: (G,2 nu s + g, g,2 nu B)-WGM */
        while (true) {
            boolean changed_R = false;
            boolean changed_L = false;
            while ((lambdaL - lambdaR) > epsilon) {
                double lambdaM = (lambdaL + lambdaR) / 2.0D;
                /** get forest for the lambda_m */
                forest = PCSF_GW(edges, getCostsLambda(lambdaM), prizesPi, g);
                if ((forest.costF >= 2.0D * C) && (forest.costF <= nu * C)) {
                    /** the forest always has c(F) \leq \nu \cdot C */
                    return forest;
                }
                if (forest.costF > nu * C) {
                    /** the cost of F_r always has c(F_r) \geq \nu \cdot C */
                    lambdaR = lambdaM;
                    changed_R = true;
                } else {
                    /** the cost of F_l always has c(F_l) \leq 2C */
                    lambdaL = lambdaM;
                    changed_L = true;
                }
            }
            if (changed_R && changed_L) {
                break;
            } else {
                F tmp_f = PCSF_GW(edges, getCostsLambda(lambdaL), prizesPi, g);
                if (tmp_f.nodesInF.size() < (2.0 * this.s * this.nu + this.g)) {
                    return tmp_f;
                }
                lambdaR = 0.0D;
                lambdaL = 2.0D * lambdaUnit;
                lambdaUnit = lambdaL;
            }
        }
        /** this forest always satisfies the WGM model. */
        return PCSF_GW(edges, getCostsLambda(lambdaL), prizesPi, g);
    }

    /**
     * we need to check the validation of our result. We assume the norm
     * ||b-bSPrime|| is the minimum one. (This means you need to know which
     * subset is optimal)
     *
     * @return true the equation 8 is satisfied. ; if it returns false, the
     * algorithm must have error(s).
     */
    private boolean checkEqu8Valid(double cT, double[] b, double[] bS, double[] bSPrime) {
        /** if there is no true subgraph the equation will not be checked. */
        if (this.trueSubGraph == null) {
            return true;
        }
        boolean flag = false;
        double leftNorm = 0.0D;
        boolean isNull = (b == null || bS == null || bSPrime == null);
        boolean isZeroLen = (b.length == 0 || bS.length == 0 || bSPrime.length == 0);
        if (isNull || isZeroLen) {
            return flag;
        }
        for (int i = 0; i < b.length; i++) {
            leftNorm += (b[i] - bS[i]) * (b[i] - bS[i]);
        }
        leftNorm = Math.sqrt(leftNorm);
        double rightNorm = 0.0D;
        for (int i = 0; i < b.length; i++) {
            rightNorm += (b[i] - bSPrime[i]) * (b[i] - bSPrime[i]);
        }
        rightNorm = Math.sqrt(rightNorm);
        if (leftNorm <= cT * rightNorm) {
            flag = true;
        } else {
            flag = false;
        }
        return flag;
    }

    /**
     * @return the minimum positive prize entry in pi vector.
     */
    private double getMinPi() {
        double minPi = Double.MAX_VALUE;
        for (int i = 0; i < prizesPi.size(); i++) {
            if (prizesPi.get(i) < minPi && prizesPi.get(i) > 0.0D) {
                minPi = prizesPi.get(i);
            }
        }
        return minPi;
    }

    /**
     * @return the summation of pi, that is \pi(G) = \sum_{i \in G} pi(i)
     */
    private double getSumPrizePi() {
        double sumPrizePi = 0.0D;
        for (int i = 0; i < prizesPi.size(); i++) {
            sumPrizePi += prizesPi.get(i);
        }
        return sumPrizePi;
    }

    /**
     * @param lambda the parameter lambda
     * @return c_{\lambda}(e) = \lambda*c(e) \forall e \in G.
     */
    private ArrayList<Double> getCostsLambda(double lambda) {
        ArrayList<Double> cLambda = new ArrayList<Double>();
        for (int i = 0; i < edgeCostsc.size(); i++) {
            cLambda.add(lambda * edgeCostsc.get(i));
        }
        return cLambda;
    }

    /**
     * The PCSF algorithm of GW pruning version
     *
     * @param edges   the edges contain in the graph.
     * @param cLambda the edge costs for current PCSF.
     * @param pi      the prize for the FCSF.
     * @param g       the number of trees(components) returned fast FCSF algorithm.
     * @return the corresponding forest F
     */
    private F PCSF_GW(ArrayList<Integer[]> edges, ArrayList<Double> cLambda, ArrayList<Double> pi, int g) {
        FastPCST pcstFast = new FastPCST(edges, pi, cLambda, FastPCST.kNoRoot, g, -1);
        ArrayList<Integer> nodesInF = null;
        ArrayList<Integer> resultEdges = null;
        if (!pcstFast.run()) {
            String errorMes = "Error: There must be an error in PCSF. \n";
            System.out.println(errorMes);
            System.exit(0);
        } else {
            /** indices of edges. */
            resultEdges = pcstFast.resultEdges;
            nodesInF = pcstFast.resultNodes;
        }
        ArrayList<Integer[]> edgesInF = new ArrayList<Integer[]>();
        ArrayList<Double> costsInF = new ArrayList<Double>();
        if (resultEdges != null) {
            for (int i : resultEdges) {
                edgesInF.add(edges.get(i));
                costsInF.add(edgeCostsc.get(i));
            }
        }
        ArrayList<Double> piInF = new ArrayList<Double>();
        for (int i : nodesInF) {
            piInF.add(pi.get(i));
        }
        double totalPrizesInG = 0.0D;
        for (int i = 0; i < pi.size(); i++) {
            totalPrizesInG += pi.get(i);
        }
        F forest = new F(nodesInF, edgesInF, costsInF, piInF, totalPrizesInG);
        if (forest.gamma != g) {
            String errorMes = "Error: # of trees in forest is not equal to g !!";
            System.out.println(errorMes);
            System.exit(0);
        }
        if (verboseLevel >= 1) {
            System.out.println("c(F): " + forest.costF);
            System.out.println("2*C: " + 2.0D * C);
            System.out.println("pi(\bar{F}): " + forest.piFBar);
        }
        return forest;
    }

    /**
     * Use disjoint set to check whether this graph is connected.
     *
     * @param edges
     * @return true is the graph is connected , false : the graph is not
     * connected
     */
    private boolean isConnected(ArrayList<Integer[]> edges) {

        DisjointSet<Integer> dis = new DisjointSet<Integer>();
        Set<Integer> nodes = new HashSet<Integer>();
        for (Integer[] edge : edges) {
            nodes.add(edge[0]);
            nodes.add(edge[1]);
        }
        for (Integer node : nodes) {
            dis.makeSet(node);
        }
        for (Integer[] edge : edges) {
            dis.union(edge[0], edge[1]);
        }
        if (dis.numConnectedComponents == 1) {
            return true;
        } else {
            return false;
        }
    }

    /**
     * Forest represents the PCSF result.
     *
     * @author baojian
     */
    public class F {

        public final ArrayList<Integer> nodesInF;
        public final ArrayList<Integer[]> edgesInF;
        public final ArrayList<Double> edgesCostInF;
        public final ArrayList<Double> prizePiInF;
        public final double costF;
        public final double piFBar;
        public final double prizeF;
        public final ArrayList<Tree> trees;
        public ArrayList<ArrayList<Integer[]>> pruningTrees;
        public ArrayList<Integer> pruningNodes;
        /**
         * the number of trees in the forest.
         */
        public int gamma;

        public F(ArrayList<Integer> nodesInF, ArrayList<Integer[]> edgesInF, ArrayList<Double> edgesCostInF,
                 ArrayList<Double> prizePiInF, double totalPrizesInG) {
            this.nodesInF = nodesInF;
            this.edgesInF = edgesInF;
            this.prizePiInF = prizePiInF;
            this.edgesCostInF = edgesCostInF;
            this.costF = getCostF();
            double result = 0.0D;
            for (double d : prizePiInF) {
                result += d;
            }
            this.prizeF = result;
            this.piFBar = totalPrizesInG - result;
            this.gamma = this.getConnectedComponents();
            this.trees = this.constructTrees();
        }

        /**
         * construct a forest from trees
         *
         * @param trees the trees are combined to be a new forest
         */
        public F(ArrayList<Tree> trees) {
            this.nodesInF = new ArrayList<Integer>();
            this.edgesInF = new ArrayList<Integer[]>();
            this.edgesCostInF = new ArrayList<Double>();
            this.prizePiInF = new ArrayList<Double>();

            if (trees == null) {
                new IllegalArgumentException("Input trees are null ...");
                System.exit(0);
            }
            /** check duplicated nodes */
            HashSet<Integer> allnodes = new HashSet<Integer>();
            for (Tree tree : trees) {
                for (Integer node : tree.nodesInT) {
                    /** check duplicated nodes. */
                    if (allnodes.add(node)) {
                    } else {
                        String errorMes = "Error: duplicated nodes found in pruneForest.";
                        System.out.println(errorMes);
                        System.exit(0);
                    }
                }
            }
            for (Tree tree : trees) {
                this.nodesInF.addAll(tree.nodesInT);
                this.prizePiInF.addAll(tree.prizePiInT);
                /** check the null object. */
                if (tree.edgesInT != null) {
                    this.edgesInF.addAll(tree.edgesInT);
                    this.edgesCostInF.addAll(tree.edgesCostInT);
                }
            }
            this.costF = getCostF();
            this.prizeF = getPrizeF();
            /** be careful : this variable will not be used. just for test. */
            this.piFBar = 0.0D;
            this.trees = constructTrees();
            this.gamma = this.getNumConnectedComponents();
            /** gamma is the number of trees in this forest. */
            if (this.gamma != this.trees.size()) {
                System.out.println("gamma is " + gamma + ", but it is not equal to " + trees.size());
                System.out.println("Error : the number of trees is not equal to gamma function ...");
                System.out.println("trees size : " + trees.size());
                System.exit(0);
            }
        }

        /**
         * @return the cost of this forest
         */
        private double getCostF() {
            double result = 0.0D;
            for (double d : this.edgesCostInF) {
                result += d;
            }
            return result;
        }

        /**
         * using disjoint set data structure to find number of connected
         * component in a forest
         *
         * @return number of connected components in F
         */
        private int getNumConnectedComponents() {
            DisjointSet<Integer> dis = new DisjointSet<Integer>();
            for (Integer node : nodesInF) {
                dis.makeSet(node);
            }
            for (Integer[] edge : edgesInF) {
                dis.union(edge[0], edge[1]);
            }
            return dis.numConnectedComponents;
        }

        /**
         * @return the number of connected components in the graph
         */
        private int getConnectedComponents() {
            DisjointSet<Integer> dis = new DisjointSet<Integer>();
            for (Integer node : nodesInF) {
                dis.makeSet(node);
            }
            for (Integer[] edge : edgesInF) {
                dis.union(edge[0], edge[1]);
            }
            return dis.numConnectedComponents;
        }

        /**
         * construct trees for forest F
         *
         * @return
         */
        public ArrayList<Tree> constructTrees() {

            DisjointSet<Integer> dis = new DisjointSet<Integer>();
            for (Integer node : nodesInF) {
                dis.makeSet(node);
            }
            for (Integer[] edge : edgesInF) {
                dis.union(edge[0], edge[1]);
            }
            /** all of components in F */
            HashMap<Integer, Set<Integer>> componentsMap = dis.getConnectedComponents();

            ArrayList<Tree> trees = new ArrayList<Tree>();
            /** for each component, create a new tree. */
            for (Integer key : componentsMap.keySet()) {
                Tree tree;
                /** nodes in this component. */
                Set<Integer> nodes = componentsMap.get(key);
                /** get all nodes in a tree. */
                ArrayList<Integer> nodesInT = new ArrayList<Integer>(nodes);
                ArrayList<Integer[]> edgesInT = new ArrayList<Integer[]>();
                ArrayList<Double> edgesCostInT = new ArrayList<Double>();
                ArrayList<Double> prizePiInT = new ArrayList<Double>();
                /** get all edges in a tree. */
                for (Integer[] edge : this.edgesInF) {
                    if (nodes.contains(edge[0]) || nodes.contains(edge[1])) {
                        edgesInT.add(edge);
                        int index = this.edgesInF.indexOf(edge);
                        /** get all edge cost in a tree. */
                        edgesCostInT.add(this.edgesCostInF.get(index));
                    }
                }
                for (Integer node : nodesInT) {
                    int index = this.nodesInF.indexOf(node);
                    prizePiInT.add(this.prizePiInF.get(index));
                }
                tree = new Tree(nodesInT, edgesInT, edgesCostInT, prizePiInT);
                trees.add(tree);
            }
            return trees;
        }

        /**
         * @return total prizes of this forest
         */
        private double getPrizeF() {
            double result = 0.0D;
            if (prizePiInF != null) {
                for (double d : prizePiInF) {
                    result += d;
                }
            }
            return result;
        }

        public int[] pruningTrees(double[] y, ArrayList<Integer> tailNodes) {
            ArrayList<Integer> allNodes = new ArrayList<Integer>();
            ArrayList<ArrayList<Integer[]>> allEdges = new ArrayList<ArrayList<Integer[]>>();
            for (Tree tree : trees) {
                HashMap<Integer, HashMap<Integer, Double>> adj = tree.adj;
                while (true) {
                    boolean flag = false;
                    HashMap<Integer, HashMap<Integer, Double>> updatedAdj = new HashMap<Integer, HashMap<Integer, Double>>();
                    for (int k : adj.keySet()) {
                        updatedAdj.put(k, new HashMap<Integer, Double>());
                        for (int kk : adj.get(k).keySet()) {
                            updatedAdj.get(k).put(kk, 0.0D);
                        }
                    }
                    for (int i : adj.keySet()) {
                        int nei = tree.isLeafNode(i);
                        /** i is a leaf node and i is normal node */
                        if (nei != -1) {
                            System.out.println("leaf : " + i + "y[" + i + "]=" + y[i]);
                            /** i is a leaf node and it is normal node. */
                            if (y[i] == 0.0D) {
                                updatedAdj.get(nei).remove(i);
                                updatedAdj.remove(i);
                                flag = true;
                                /** i is a leaf node and it is abnormal node */
                            } else {
                                int backNode = i;
                                int next = nei;
                                ArrayList<Integer> path = new ArrayList<Integer>();
                                path.add(backNode);

                                while (adj.get(next).size() == 2) {
                                    ArrayList<Integer> neis = new ArrayList<Integer>(adj.get(next).keySet());
                                    int nei0 = neis.get(0);
                                    int nei1 = neis.get(1);
                                    if (nei0 == backNode) {
                                        backNode = next;
                                        next = nei1;
                                    } else {
                                        backNode = next;
                                        next = nei0;
                                    }
                                    path.add(backNode);
                                }
                                /** update adj */
                                boolean flag1 = false;
                                for (int kk : adj.keySet()) {
                                    if (adj.get(kk).size() > 2) {
                                        flag1 = true;
                                        break;
                                    }
                                }
                                if (path.size() >= 3 && flag1 == true) {
                                    double ratio = 0.0D;
                                    for (int k : path) {
                                        if (y[k] == 0.0D) {
                                            ratio += 1.0D;
                                        }
                                    }
                                    ratio = ratio / (path.size() + 0.0D);
                                    System.out.println("ratio : " + ratio + " path is : " + path.toString());
                                    if (ratio > 0.5) {
                                        int lastNode = path.get(path.size() - 1);
                                        int nextLast = path.get(path.size() - 2);
                                        updatedAdj.get(lastNode).remove(nextLast);
                                        path.remove(path.size() - 1);
                                        for (int pathNode : path) {
                                            updatedAdj.remove(pathNode);
                                        }
                                        flag = true;
                                    }
                                }
                            }
                        }
                    }
                    if (flag == false) {
                        break;
                    }
                    adj = updatedAdj;
                    System.out.println("end of iteration ...");
                }
                ArrayList<Integer[]> currentTree = new ArrayList<Integer[]>();
                for (int key : adj.keySet()) {
                    allNodes.add(key);
                    for (int k : adj.get(key).keySet()) {
                        currentTree.add(new Integer[]{key, k});
                    }
                }
                allEdges.add(currentTree);
            } /** next tree */
            int[] result = new int[allNodes.size()];
            for (int i = 0; i < allNodes.size(); i++) {
                result[i] = allNodes.get(i);
            }
            this.pruningTrees = allEdges;
            this.pruningNodes = allNodes;
            return result;
        }
    }

    /**
     * The tree represented by edges and nodes.
     *
     * @author baojian
     */
    public class Tree implements Comparable<Tree> {
        public final ArrayList<Integer> nodesInT;
        /**
         * Note: edges in T are undirected. Therefore, each of edge only save
         * once.
         */
        public final ArrayList<Integer[]> edgesInT;
        public final ArrayList<Double> edgesCostInT;
        public final ArrayList<Double> prizePiInT;
        public final double costTree;
        public final double prizeTree;
        public final double ratio;
        /**
         * adj, the adjacency list for path computing
         */
        public final HashMap<Integer, HashMap<Integer, Double>> adj;

        /**
         * Construct a new tree with more than one node
         *
         * @param nodesInT
         * @param edges
         * @param edgesCost
         * @param prizePi
         */
        public Tree(ArrayList<Integer> nodesInT, ArrayList<Integer[]> edges, ArrayList<Double> edgesCost,
                    ArrayList<Double> prizePi) {
            this.nodesInT = nodesInT;
            this.edgesInT = edges;
            this.edgesCostInT = edgesCost;
            this.prizePiInT = prizePi;
            this.costTree = getCostTree();
            this.prizeTree = getPrizeTree();
            this.ratio = this.prizeTree / this.costTree;
            this.adj = contructAdj();
        }

        /**
         * construct a new tree with only one node
         *
         * @param i
         * @param pi
         */
        public Tree(int i, double pi) {
            ArrayList<Integer> nodes = new ArrayList<Integer>();
            nodes.add(i);
            ArrayList<Integer[]> edges = new ArrayList<Integer[]>();
            ArrayList<Double> edgesCost = new ArrayList<Double>();
            ArrayList<Double> prizePi = new ArrayList<Double>();
            prizePi.add(pi);
            this.nodesInT = nodes;
            this.edgesInT = edges;
            this.edgesCostInT = edgesCost;
            this.prizePiInT = prizePi;
            this.costTree = getCostTree();
            this.prizeTree = getPrizeTree();
            this.ratio = this.prizeTree / this.costTree;
            this.adj = contructAdj();
        }

        public int isLeafNode(int node) {
            if (this.adj.get(node).size() == 1) {
                for (int nei : this.adj.get(node).keySet()) {
                    return nei;
                }
                return -1;
            } else {
                return -1;
            }
        }

        /**
         * @return adjacency list of this tree <node_i,nei(node_i)>
         */
        private HashMap<Integer, HashMap<Integer, Double>> contructAdj() {
            HashMap<Integer, HashMap<Integer, Double>> adj = new HashMap<Integer, HashMap<Integer, Double>>();
            for (int node : this.nodesInT) {
                HashMap<Integer, Double> nei = new HashMap<Integer, Double>();
                adj.put(node, nei);
            }
            if (this.edgesCostInT != null) {
                for (int i = 0; i < this.edgesCostInT.size(); i++) {
                    Integer[] edge = this.edgesInT.get(i);
                    adj.get(edge[0]).put(edge[1], this.edgesCostInT.get(i));
                    adj.get(edge[1]).put(edge[0], this.edgesCostInT.get(i));
                }
            }
            return adj;
        }

        /**
         * @return the total cost of this tree
         */
        private double getCostTree() {
            double result = 0.0D;
            if (this.edgesCostInT != null) {
                for (double d : this.edgesCostInT) {
                    result += d;
                }
            }
            return result;
        }

        /**
         * @return the total prize of this tree
         */
        private double getPrizeTree() {
            double result = 0.0D;
            if (prizePiInT != null) {
                for (double d : prizePiInT) {
                    result += d;
                }
            }
            return result;
        }

        /**
         * @return a one node tree with this tree's maximum prize
         */
        public Tree getOneNodeMaxTree() {
            double maximumPi = -Double.MAX_VALUE;
            int nodeID = -1;
            for (int i = 0; i < prizePiInT.size(); i++) {
                double d = prizePiInT.get(i);
                if (d > maximumPi) {
                    maximumPi = d;
                    nodeID = nodesInT.get(i);
                }
            }
            ArrayList<Integer> nodeInT = new ArrayList<Integer>();
            nodeInT.add(nodeID);
            ArrayList<Double> prizePi = new ArrayList<Double>();
            prizePi.add(maximumPi);
            return new Tree(nodeInT, null, null, prizePi);
        }

        /**
         * we simulate a tree as directed graph, and each edge in the tree are
         * simulated as two edges with different edge direction. For example, an
         * edge [i,j] in this tree will be simulated as [i->j] and [j->i].
         * Therefore, we can construct a connected Directed Graph. After that,
         * we can find a tour trough depth first search algorithm.
         *
         * @return a Euler tour with (2|T| -1) nodes.
         */
        public ArrayList<Integer> getEulerTour() {
            ArrayList<Integer> tour = new ArrayList<Integer>();
            Digraph digraph = new Digraph(nodesInT.size());
            ArrayList<Integer[]> directedEdges = new ArrayList<Integer[]>();
            for (Integer[] edge : edgesInT) {
                int index0 = nodesInT.indexOf(edge[0]);
                int index1 = nodesInT.indexOf(edge[1]);
                Integer[] edge1 = new Integer[]{index0, index1};
                Integer[] edge2 = new Integer[]{index1, index0};
                directedEdges.add(edge1);
                directedEdges.add(edge2);
            }
            for (Integer[] directedEdge : directedEdges) {
                digraph.addEdge(directedEdge[0], directedEdge[1]);
            }
            DirectedEulerianCycle di = new DirectedEulerianCycle(digraph);
            for (Iterator<Integer> iterator = di.cycle().iterator(); iterator.hasNext(); ) {
                int indexNode = iterator.next();
                tour.add(nodesInT.get(indexNode));
            }
            if (tour.size() != (2 * this.nodesInT.size() - 1)) {
                System.out.println("Error : the length of the tour should be (2*|V_T| -1 )");
                System.exit(0);
            }
            return tour;
        }

        // @Override
        public int compareTo(Tree o) {
            if ((this.ratio - o.ratio) < 0) {
                return -1;
            } else if ((this.ratio - o.ratio) == 0.0D) {
                return 0;
            } else {
                return 1;
            }
        }
    }// Tree class
}
