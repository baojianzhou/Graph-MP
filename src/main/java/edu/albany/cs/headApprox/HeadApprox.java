package edu.albany.cs.headApprox;

import edu.albany.cs.base.DisjointSet;
import edu.albany.cs.fastPCST.FastPCST;
import org.apache.commons.lang3.ArrayUtils;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

/**
 * Head approximation
 * 
 * Aglorithm 3 of
 * "http://people.csail.mit.edu/ludwigs/papers/icml15_graphsparsity.pdf" Title :
 * "A Nearly-Linear Time Framework for Graph-Structured Sparsity" Authors :
 * Ludwig Schmidt, Chinmay Hegde, and Piotr Indyk
 *
 * @author baojian bzhou6@albany.edu
 */
public class HeadApprox {
	/** Graph G is split into 3 parts edges, c, and pi. */
	/** the edges in graph G */
	private ArrayList<Integer[]> edges;
	/** the edge costs c corresponding to edges */
	private ArrayList<Double> c;
	/** the node prizes pi corresponding to nodes {0,...,n-1} */
	private ArrayList<Double> pi;
	/** the number of connected components */
	private int g;
	/** the cost budget */
	private double C;
	/** the delta value, which is a constant value 1.0/169.0 */
	private double delta;

	/** additional parameters, may not be useful. */
	private double cH;
	private int[] trueSubGraph;
	private int verboseLevel = 0;

	/** the result forest that is returned by Tail algorithm. */
	public F bestForest;
	/** This parameter "valid" to show the result is valid. */
	/** But, it may not be useful. */
	public boolean valid;

	/** This constructor is just for test. */
	public HeadApprox() {
	}

	/**
	 * A general input constructor.
	 * 
	 * @param edges
	 *            the edges of graph G.
	 * @param c
	 *            the edge costs corresponding to the edges in G.
	 * @param pi
	 *            the node prize in pi
	 * @param g
	 *            the number of connected components in the result forest.
	 * @param C
	 *            the cost budget
	 */
	public HeadApprox(ArrayList<Integer[]> edges, ArrayList<Double> c, ArrayList<Double> pi, int g, double C) {

		this.edges = edges;
		this.c = c;
		this.pi = pi;
		this.C = C;
		this.g = g;
		/** By Theorem 11 in this paper (p.20). */
		this.delta = 1.0D / 169.0D; //
		this.cH = Math.sqrt(1.0D / 14.0D);
		this.checkInputValidation(edges, pi, c);
		this.bestForest = run();
		double[] b = new double[pi.size()];
		for (int i = 0; i < pi.size(); i++) {
			b[i] = pi.get(i);
		}
		if (this.checkEqu9Valid(b)) {
			this.valid = true;
		} else {
			String errorMes = "Head approximation is invalid ...";
			System.out.println(errorMes);
			System.exit(0);
			this.valid = false;
		}
	}

	/**
	 * This constructor is for Algorithm 3 only
	 * 
	 * @param edges
	 *            the edges corresponding to WGM.
	 * @param edgeCostsW
	 *            the edge costs corresponding to WGM.
	 * @param b
	 *            the real vector that we want to estimate.
	 * @param s
	 *            the sparsity parameter.
	 * @param g
	 *            the number of connected components.
	 * @param B
	 *            the budget constraint, it is a real number.
	 * @param trueSubGraph
	 *            let it be null if you do not need it.
	 */
	public HeadApprox(ArrayList<Integer[]> edges, ArrayList<Double> edgeCostsW, double[] b, int s, int g, double B,
			int[] trueSubGraph) {
		this.edges = edges;
		/** edge costs c, let c(e) = w(e) + (B / s) */
		c = new ArrayList<Double>();
		for (double w : edgeCostsW) {
			c.add(w + B / (s * 1.0D));
		}
		/** let pi(i) = bi*bi */
		pi = new ArrayList<Double>();
		if (b == null || b.length == 0) {
			System.out.println("Error: vector z is null ...");
			System.exit(0);
		}
		for (int i = 0; i < b.length; i++) {
			pi.add(b[i] * b[i]);
		}
		/** let C = 2*B */
		C = 2.0D * B;
		/** we set this constant value by Theorem 11 of that paper. */
		delta = 1.0D / 169.0D;
		this.g = g;
		/** They may not be useful. */
		cH = Math.sqrt(1.0D / 14.0D);
		this.trueSubGraph = trueSubGraph;
		/** we skip the validation of equation 9 */
		this.trueSubGraph = null;
		if (verboseLevel >= 1) {
			System.out.println("pi : ");
			for (int k = 0; k < pi.size(); k++) {
				System.out.format(",%.1f", pi.get(k));
			}
			System.out.println();
			System.out.println("c(e): " + c.toString());
			System.out.println("B: " + B);
			System.out.println("C: " + C);
			System.out.println("delta: " + delta);
			System.out.println("cH: " + cH);
		}
		if (!isParametersValid()) {
			System.out.println("some parameters are not valid.");
			System.exit(0);
		}
		bestForest = run();
		if (checkEqu9Valid(b)) {
			valid = true;
		} else {
			System.out.println("Head approximation is invalid ...");
			System.exit(0);
			valid = false;
		}
	}

	/**
	 * check input validation, to make sure the following properties : 1. the
	 * graph is connected. 2. the graph does not have self cyclic edge
	 *
	 * @param edges
	 * @param pi
	 */
	private void checkInputValidation(ArrayList<Integer[]> edges, ArrayList<Double> pi, ArrayList<Double> edgeCosts) {
		/** check input validation */
		DisjointSet<Integer> dis = new DisjointSet<Integer>();
		HashSet<Integer> nodes = new HashSet<Integer>();
		for (Integer[] edge : edges) {
			nodes.add(edge[0]);
			nodes.add(edge[1]);
			dis.makeSet(edge[0]);
			dis.makeSet(edge[1]);
			dis.union(edge[0], edge[1]);
			if (edge[0].intValue() == edge[1].intValue()) {
				new IllegalArgumentException("bad edge : [" + edge[0] + "," + edge[1] + "]");
			}
		}
		if (dis.numConnectedComponents != 1) {
			new IllegalArgumentException("Error : the graph is not connected ...");
			System.exit(0);
		}
		if (nodes.size() != pi.size()) {
			new IllegalArgumentException("Error : edges of graph do not have whole nodes  ...");
			System.exit(0);
		}
		/** check validation of every weight of edge in the graph */
		for (double edgecost : edgeCosts) {
			if (edgecost <= 1e-9D) {
				new IllegalArgumentException("bad edge : " + edgecost);
				System.exit(0);
			}
		}
	}

	private boolean isParametersValid() {

		/** to check the graph is a connected graph. */
		if (!isConnected(edges)) {
			String errorMes = "the graph is not connected.";
			System.out.println(errorMes);
			return false;
		}

		/** check edgeCosts are valid */
		for (double edge : c) {
			if (edge <= 0.0D) {
				new IllegalAccessException("the edge cost should not be non-positive ...");
				System.exit(0);
			}
		}
		return true;
	}

	/**
	 * This is the head algorithm part.
	 *
	 * @return F a forest is returned by Head algorithm.
	 */
	private F run() {
		/** the minimum positive prize entry in pi vector. */
		double minPi = getMinPi();
		/** lambda_r = the minimum binary search value. */
		BigDecimal lambdaR = new BigDecimal((2.0D * C) / minPi);
		F forest = PCSF_GW(edges, c, getLambdaPi(lambdaR.doubleValue()), g);
		/** make sure we have the invariant c(Fr) > 2*C */
		if (forest.costF <= (2.0D * C)) {
			return forest;
		}
		BigDecimal epsilon = new BigDecimal((delta * C) / (2.0D * getSigmaPi()));
		BigDecimal lambdaL = new BigDecimal(1.0D / (4.0D * getSigmaPi()));
		if (verboseLevel >= 1) {
			System.out.println("costF : " + forest.costF);
			System.out.println("2*C: " + 2.0D * C);
			System.out.println("min pi : " + minPi);
			System.out.println("lambdaR : " + lambdaR.doubleValue());
			System.out.println("epsilon : " + epsilon.doubleValue());
			System.out.println("lambdaL : " + lambdaL.doubleValue());
		}
		int iter = 0;
		/** Binary search over the Lagrange parameter lambda */
		/** (lambda_r -lambda_l > epsilon) */
		while (lambdaR.subtract(lambdaL).compareTo(epsilon) == 1) {
			BigDecimal lambdaM = (lambdaL.add(lambdaR)).divide(new BigDecimal(2.0D));
			forest = PCSF_GW(edges, c, getLambdaPi(lambdaM.doubleValue()), g);
			if (forest.costF > (2.0D * C)) {
				lambdaR = lambdaM;
			} else {
				lambdaL = lambdaM;
			}
			if (verboseLevel >= 1) {
				iter++;
				System.out.println("iteration: " + iter);
				System.out.println("lambdaM: " + lambdaM);
				System.out.println("epsilon: " + epsilon);
				System.out.println("lambdaR: " + lambdaR);
				System.out.println("lambdaL: " + lambdaL);
				System.out.println("lambdaM: " + lambdaM);
			}
		}
		if (verboseLevel >= 1) {
			System.out.println("final lambdaL: " + lambdaL.doubleValue());
			System.out.println("final lambdaR: " + lambdaR.doubleValue());
			System.out.println("total iterations: " + iter);
		}
		F forestL = PCSF_GW(edges, c, getLambdaPi(lambdaL.doubleValue()), g);
		F forestR = PCSF_GW(edges, c, getLambdaPi(lambdaR.doubleValue()), g);
		if (verboseLevel >= 1) {
			System.out.println("nodes in L : " + forestL.nodesInF.size());
			System.out.println("edges in L :" + forestL.edgesInF.size());
			System.out.println("CC : " + forestL.gamma);
			System.out.println("nodes in R : " + forestR.nodesInF.size());
			System.out.println("edges in R :" + forestR.edgesInF.size());
			System.out.println("CC : " + forestR.gamma);
		}
		/** Prune the potentially large solution Fr */
		F forestRPrime = pruneForest(forestR, c, pi, C);
		if (forestL.prizeF >= forestRPrime.prizeF) {
			return forestL;
		} else {
			return forestRPrime;
		}
	}

	/**
	 * get total prize of this graph
	 *
	 * @return the summation of pi, that is \pi(G) = \sum_{i \in G} pi(i)
	 */
	private double getSigmaPi() {
		double result = 0.0D;
		for (double p : pi) {
			result += p;
		}
		return result;
	}

	/**
	 * we need to check the validation of our result. We assume the norm
	 * ||b-bSPrime|| is the minimum one. (This means you need to know which
	 * subset is optimal)
	 *
	 * @return true the equation 9 is satisfied. ; if it returns false, the
	 *         algorithm must have error(s).
	 */
	private boolean checkEqu9Valid(double[] b) {
		/** do not need to check equation 9 */
		if (this.trueSubGraph == null) {
			return true;
		}
		double[] bS = new double[b.length];
		double[] bSPrime = new double[b.length];
		for (int i = 0; i < b.length; i++) {
			if (ArrayUtils.contains(trueSubGraph, i)) {
				bSPrime[i] = b[i];
			} else {
				bSPrime[i] = 0.0D;
			}
			if (this.bestForest.nodesInF.contains(i)) {
				bS[i] = b[i];
			} else {
				bS[i] = 0.0D;
			}
		}
		boolean flag = false;
		if (bS == null || bSPrime == null || bS.length == 0 || bSPrime.length == 0) {
			return flag;
		}
		double leftNorm = 0.0D;
		for (int i = 0; i < bS.length; i++) {
			leftNorm += bS[i] * bS[i];
		}
		leftNorm = Math.sqrt(leftNorm);
		double rightNorm = 0.0D;
		for (int i = 0; i < bSPrime.length; i++) {
			rightNorm += bSPrime[i] * bSPrime[i];
		}
		rightNorm = Math.sqrt(rightNorm);
		if (leftNorm >= this.cH * rightNorm) {
			flag = true;
		} else {
			flag = false;
		}
		return flag;
	}

	/**
	 * algorithm 4 Head Approximation for the WGM : subroutine PruneForest
	 *
	 * @param forest
	 * @param edgeCosts
	 * @param prizes
	 * @param costBudget
	 * @return
	 */
	private F pruneForest(F forest, ArrayList<Double> edgeCosts, ArrayList<Double> prizes, double costBudget) {
		ArrayList<Tree> trees = forest.getDescendingSortedTrees();
		double costBudgetR = costBudget;
		ArrayList<Tree> treePrimes = new ArrayList<Tree>();
		for (int i = 0; i < forest.size(); i++) {
			Tree treePrimei;
			Tree treei = trees.get(i);
			if (costBudgetR >= treei.costTree) {
				treePrimei = treei;
				treePrimes.add(treePrimei);
				/** cost budget C^i = c(T_i) */
				costBudgetR = costBudgetR - treei.costTree;
			} else if (costBudgetR > 0) {
				treePrimei = this.pruneTree(treei, edgeCosts, prizes, costBudgetR);
				treePrimes.add(treePrimei);
				/** cost budget C^i = Cr */
				costBudgetR = 0.0D;
			} else {
				treePrimei = treei.getOneNodeMaxTree();
				/** cost budget C^i = 0 */
				treePrimes.add(treePrimei);
			}
		}
		return new F(treePrimes);
	}

	private Tree pruneTree(Tree tree, ArrayList<Double> edgeCosts, ArrayList<Double> pi, double CPrime) {

		/** T = (Vt,Et), the method of getTour has been tested */
		ArrayList<Integer> tourL = tree.getEulerTour();
		ArrayList<Double> piPrime = new ArrayList<Double>();
		HashSet<Integer> hashSet = new HashSet<Integer>();
		for (int j = 0; j < tourL.size(); j++) {
			int nodeI = tourL.get(j);
			int indexNodeI = tree.nodesInT.indexOf(nodeI);
			/** if position j is the first appearance of vj in L */
			if (hashSet.add(nodeI)) {
				if (pi.get(nodeI) != tree.prizePiInT.get(indexNodeI)) {
					System.out.println("the prize of this node is inconsistent ...");
					System.out.println(pi.get(nodeI) + "is not equal to " + tree.prizePiInT.get(indexNodeI));
					System.exit(0);
				}
				piPrime.add(pi.get(nodeI));
			} else {
				piPrime.add(0.0D);
			}
		}
		double phi = tree.prizeTree / tree.costTree;
		for (int i = 0; i < tree.nodesInT.size(); i++) {
			int nodeI = tree.nodesInT.get(i);
			if (pi.get(nodeI) != tree.prizePiInT.get(i)) {
				System.out.println("the prize of this node is inconsistent ...");
				System.out.println(pi.get(i) + "is not equal to " + tree.prizePiInT.get(tree.nodesInT.indexOf(i)));
				System.exit(0);
			}
			if (pi.get(nodeI) >= (CPrime * phi / 6.0D)) {
				return new Tree(nodeI, pi.get(i));
			}
		}

		int l = 1;
		ArrayList<Integer> Pl = new ArrayList<Integer>();
		ArrayList<Integer> PlPlus1 = new ArrayList<Integer>();
		for (int i = 0; i < tourL.size(); i++) {
			Pl.add(i);
			double CPrimeInPl = 0.0D;
			double piPrimeInPl = 0.0D;
			/** Pl has at least two nodes. It has at least one edge */
			if (Pl.size() >= 2) {
				for (int ii = 0; ii < Pl.size() - 1; ii++) {
					int Pi = tourL.get(Pl.get(ii));
					int PiPlus1 = tourL.get(Pl.get(ii + 1));
					CPrimeInPl += tree.adj.get(Pi).get(PiPlus1);
				}
			}
			/** Pl has at least one node, any path has at least one node */
			if (Pl.size() >= 1) {
				for (int ii = 0; ii < Pl.size(); ii++) {
					piPrimeInPl += piPrime.get(Pl.get(ii));
				}
			}
			if (CPrimeInPl > CPrime) {
				l++;
				PlPlus1 = Pl;
				Pl = new ArrayList<Integer>();
			} else if (piPrimeInPl >= (CPrime * phi / 6.0D)) {
				/** Pl has only one node */
				if (Pl.size() == 1) {
					int currentNode = tourL.get(Pl.get(0));
					int index = tree.nodesInT.indexOf(currentNode);
					if (pi.get(currentNode) != tree.prizePiInT.get(index)) {
						System.out.println("the prize of this node is inconsistent ...");
						System.out.println(
								pi.get(i) + "is not equal to " + tree.prizePiInT.get(tree.nodesInT.indexOf(i)));
						System.exit(0);
					}
					/** create single node tree (a tree without any edge) */
					return new Tree(currentNode, pi.get(currentNode));
				} else if (Pl.size() >= 2) {
					ArrayList<Integer> nodesInT = new ArrayList<Integer>();
					ArrayList<Integer[]> edgesInT = new ArrayList<Integer[]>();
					ArrayList<Double> edgesCostInT = new ArrayList<Double>();
					ArrayList<Double> prizePiInT = new ArrayList<Double>();
					for (int ii = 0; ii < Pl.size() - 1; ii++) {
						int Pi = Pl.get(ii);
						int PiPlus1 = Pl.get(ii + 1);
						boolean flag0 = false;
						boolean flag1 = false;
						int nodeI = tourL.get(Pi);
						int nodePlusI = tourL.get(PiPlus1);
						if (!nodesInT.contains(nodeI)) {
							flag0 = true;
							nodesInT.add(nodeI);
							prizePiInT.add(pi.get(nodeI));
							/** check pi validation */
							if (pi.get(nodeI) != tree.prizePiInT.get(tree.nodesInT.indexOf(nodeI))) {
								System.out.println("the prize of this node is inconsistent ...");
								System.out.println(
										pi.get(i) + "is not equal to " + tree.prizePiInT.get(tree.nodesInT.indexOf(i)));
								System.exit(0);
							}
						}
						if (!nodesInT.contains(nodePlusI)) {
							flag1 = true;
							nodesInT.add(nodePlusI);
							prizePiInT.add(pi.get(nodePlusI));
							/** check pi validation */
							if (pi.get(nodePlusI) != tree.prizePiInT.get(tree.nodesInT.indexOf(nodePlusI))) {
								System.out.println("the prize of this node is inconsistent ...");
								System.out.println(
										pi.get(i) + "is not equal to " + tree.prizePiInT.get(tree.nodesInT.indexOf(i)));
								System.exit(0);
							}
						}
						/** add new edge if not true */
						if ((flag0 == false) && (flag1 == false)) {
							/** do nothing as the edge already exists */
						} else {
							Integer[] edge = new Integer[] { nodeI, nodePlusI };
							double cost = tree.adj.get(nodeI).get(nodePlusI);
							edgesInT.add(edge);
							edgesCostInT.add(cost);
						}
					}
					/** return the subtree of tree T on the nodes in Pl. */
					return new Tree(nodesInT, edgesInT, edgesCostInT, prizePiInT);
				} else {
					System.out.println("Error : path Pl has at leat one node ...");
					System.exit(0);
				}
			}
		} // end for
		/** algorithm will never reach this point */
		this.mergePNodes(Pl, PlPlus1, l);
		return null;
	}

	/**
	 * This method could not be called.
	 *
	 * @param P1
	 * @param P2
	 * @param l
	 */
	private void mergePNodes(ArrayList<Integer> P1, ArrayList<Integer> P2, int l) {
		new IllegalAccessException("Error : the algorithm will never reach this point .... ");
		System.exit(0);
	}

	/**
	 * get minimum positive entry of pi vector
	 *
	 * @return the minimum positive prize entry in pi vector.
	 */
	private double getMinPi() {
		double minPi = Double.MAX_VALUE;
		for (int i = 0; i < pi.size(); i++) {
			if (pi.get(i) < minPi && pi.get(i) > 0.0D) {
				minPi = pi.get(i);
			}
		}
		return minPi;
	}

	/**
	 * get lambda*pi vector
	 *
	 * @param lambda
	 * @return pi*lambda
	 */
	private ArrayList<Double> getLambdaPi(double lambda) {
		ArrayList<Double> piLambda = new ArrayList<Double>();
		for (int i = 0; i < pi.size(); i++) {
			piLambda.add(lambda * pi.get(i));
		}
		return piLambda;
	}

	/**
	 * The PCSF algorithm of GW pruning version
	 *
	 * @param edges
	 *            the edges of the graph G
	 * @param edgeCosts
	 *            the edge costs in the graph
	 * @param pi
	 *            the prize for the FCSF.
	 * @param g
	 *            the number of connected components
	 * @return the corresponding forest F of pcsf algorithm
	 */
	public F PCSF_GW(ArrayList<Integer[]> edges, ArrayList<Double> edgeCosts, ArrayList<Double> pi, int g) {

		/** check prizes are valid */
		for (double p : pi) {
			if (p < 0.0D) {
				new IllegalAccessException("the prize should not be negative ...");
				System.exit(0);
			}
		}
		FastPCST pcstFast = new FastPCST(edges, pi, edgeCosts, FastPCST.kNoRoot, g, -1);
		ArrayList<Integer> nodesInF = null;
		ArrayList<Integer> resultEdges = null;
		if (!pcstFast.run()) {
			String errorMes = "Error : There must be an error. \n";
			System.out.println(errorMes);
			System.exit(0);
		} else {
			/** indexes of edges */
			resultEdges = pcstFast.resultEdges;
			nodesInF = pcstFast.resultNodes;
		}
		ArrayList<Integer[]> edgesInF = new ArrayList<Integer[]>();
		ArrayList<Double> costsInF = new ArrayList<Double>();
		if (resultEdges != null) {
			for (int i : resultEdges) {
				edgesInF.add(edges.get(i));
				costsInF.add(edgeCosts.get(i));
			}
		}
		/** the prizes should be the original prize */
		ArrayList<Double> piInF = new ArrayList<Double>();
		for (int i : nodesInF) {
			piInF.add(pi.get(i));
		}
		double totalPrizes = 0.0D;
		for (int i = 0; i < pi.size(); i++) {
			totalPrizes += pi.get(i);
		}
		if (verboseLevel >= 2) {
			System.out.println("number of nodes in F: " + nodesInF.size());
			System.out.println("nodesInF: " + nodesInF.toString());
			System.out.println("edgesInF: ");
			for (Integer[] edge : edgesInF) {
				System.out.println("," + "[" + edge[0] + "," + edge[1] + "]");
			}
			System.out.println();
		}
		/** the nodes, prizes, edges and costs must be consistent */
		F forest = new F(nodesInF, edgesInF, costsInF, piInF, totalPrizes);
		if (verboseLevel >= 2) {
			System.out.println("forest.gamma: " + forest.gamma);
			System.out.println("g: " + g);
		}
		if (forest.gamma != g) {
			String errorMes = "Error: # of trees in forest is not equal to g !!";
			System.out.println(errorMes);
			System.exit(0);
		}
		return forest;
	}

	/**
	 * Use disjoint set to check whether this graph is connected.
	 *
	 * @param edges
	 * @return true is the graph is connected , false : the graph is not
	 *         connected
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
	 * The forest represents the subset of edges found by algorithm, The forest
	 * can be seen as a combination of several trees
	 *
	 * @author baojian
	 */
	public class F {
		public final ArrayList<Integer> nodesInF;
		public final ArrayList<Integer[]> edgesInF;
		public final ArrayList<Double> edgesCostInF;
		public final ArrayList<Double> prizePiInF;
		public final double costF;
		public final double prizeF;
		public final double piFBar;
		public final int gamma;
		/** we construct corresponding trees of forest */
		public final ArrayList<Tree> trees;

		/**
		 * make sure the nodes and prizes, edges and costs are correspondingly
		 * consistent, which means the edges and nodes do not start from 0 ...
		 *
		 * @param nodesInF
		 * @param edgesInF
		 * @param edgesCostInF
		 * @param prizePiInF
		 * @param totalPrizes
		 */
		public F(ArrayList<Integer> nodesInF, ArrayList<Integer[]> edgesInF, ArrayList<Double> edgesCostInF,
				ArrayList<Double> prizePiInF, double totalPrizes) {
			this.nodesInF = nodesInF;
			this.edgesInF = edgesInF;
			this.edgesCostInF = edgesCostInF;
			this.prizePiInF = prizePiInF;
			this.costF = getCostF();
			this.prizeF = getPrizeF();
			this.trees = constructTrees();
			this.piFBar = totalPrizes - this.prizeF;
			this.gamma = this.getNumConnectedComponents();
			/** make sure gamma is the number of trees in this forest */
			if (this.gamma != this.trees.size()) {
				System.out.println("Error : the number of trees is not equal ...");
				System.exit(0);
			}
		}

		/**
		 * construct a forest from trees
		 *
		 * @param trees
		 *            the trees are combined to be a new forest
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
					if (!allnodes.add(node)) {
						System.out.println("Error : duplicated nodes found in pruneForest ...");
						System.exit(0);
					}
				}
			}
			for (Tree tree : trees) {
				this.nodesInF.addAll(tree.nodesInT);
				this.prizePiInF.addAll(tree.prizePiInT);
				/** check the null value */
				if (tree.edgesInT != null) {
					this.edgesInF.addAll(tree.edgesInT);
					this.edgesCostInF.addAll(tree.edgesCostInT);
				}
			}
			this.costF = getCostF();
			this.prizeF = getPrizeF();
			this.piFBar = 0.0D; // be careful this will not be used.
			this.trees = constructTrees();
			this.gamma = this.getNumConnectedComponents();
			/** check gamma is the number of trees in this forest */
			if (this.gamma != this.trees.size()) {
				System.out.println("Error : the number of trees is not equal to gamma function ...");
				System.out.println("gamma is " + this.gamma + " is not equal to " + this.trees.size());
				System.out.println("trees size : " + trees.size());
				System.exit(0);
			}
		}

		/**
		 * construct trees for forest F
		 *
		 * @return
		 */
		private ArrayList<Tree> constructTrees() {

			DisjointSet<Integer> dis = new DisjointSet<Integer>();
			for (Integer node : nodesInF) {
				dis.makeSet(node);
			}
			for (Integer[] edge : edgesInF) {
				dis.union(edge[0], edge[1]);
			}
			/** all of components in forest F */
			HashMap<Integer, Set<Integer>> componentsMap = dis.getConnectedComponents();
			ArrayList<Tree> trees = new ArrayList<Tree>();
			/** for each component create a new tree */
			for (Integer key : componentsMap.keySet()) {
				Tree tree;
				/** nodes in this component */
				Set<Integer> nodes = componentsMap.get(key);
				/** all nodes in a tree */
				ArrayList<Integer> nodesInT = new ArrayList<Integer>(nodes);
				ArrayList<Integer[]> edgesInT = new ArrayList<Integer[]>();
				ArrayList<Double> edgesCostInT = new ArrayList<Double>();
				ArrayList<Double> prizePiInT = new ArrayList<Double>();
				/** all edges in a tree */
				for (Integer[] edge : this.edgesInF) {
					if (nodes.contains(edge[0]) || nodes.contains(edge[1])) {
						edgesInT.add(edge);
						int index = this.edgesInF.indexOf(edge);
						/** all edge cost in a tree */
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
		 * sort the trees in this forest by descending ratio r1 > r2 > r3 > ...
		 *
		 * @return the sorted trees
		 */
		public ArrayList<Tree> getDescendingSortedTrees() {
			/** Note : with descending order */
			Collections.sort(this.trees, new Comparator<Tree>() {
				public int compare(Tree o1, Tree o2) {
					return o2.compareTo(o1);
				}
			});
			/** check duplicated nodes */
			HashSet<Integer> allnodes = new HashSet<Integer>();
			for (Tree tree : trees) {
				for (Integer node : tree.nodesInT) {
					/** do nothing */
					if (allnodes.add(node)) {
					} else {
						System.out.println("Error : duplicated nodes found in pruneForest ...");
						System.exit(0);
					}
				}
			}
			return this.trees;
		}

		/**
		 * gamma function in that paper
		 * 
		 * @return the size of the trees in Forest
		 */
		public int size() {
			return this.trees.size();
		}

		/**
		 * @return total cost of this forest
		 */
		private double getCostF() {
			double result = 0.0D;
			if (this.edgesCostInF != null) {
				for (double d : this.edgesCostInF) {
					result += d;
				}
			}
			return result;
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
	}// class F

	/**
	 * The tree represented by edges and nodes.
	 *
	 * @author baojian
	 */
	public class Tree implements Comparable<Tree> {
		public final ArrayList<Integer> nodesInT;
		/** Note: edges in T are undirected. So, each of edge only save once. */
		public final ArrayList<Integer[]> edgesInT;
		public final ArrayList<Double> edgesCostInT;
		public final ArrayList<Double> prizePiInT;
		public final double costTree;
		public final double prizeTree;
		public final double ratio;
		/** adj for path computing */
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
				Integer[] edge1 = new Integer[] { index0, index1 };
				Integer[] edge2 = new Integer[] { index1, index0 };
				directedEdges.add(edge1);
				directedEdges.add(edge2);
			}
			for (Integer[] directedEdge : directedEdges) {
				digraph.addEdge(directedEdge[0], directedEdge[1]);
			}
			DirectedEulerianCycle di = new DirectedEulerianCycle(digraph);
			for (Iterator<Integer> iterator = di.cycle().iterator(); iterator.hasNext();) {
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

	/**
	 * test for gamma and disjoint set
	 */
	public void testGammaDisjointSet() {
		ArrayList<Integer> nodes = new ArrayList<Integer>();
		nodes.add(1);
		nodes.add(2);
		nodes.add(3);
		nodes.add(4);
		nodes.add(5);
		ArrayList<Integer[]> edges = new ArrayList<Integer[]>();
		edges.add(new Integer[] { 1, 2 });
		edges.add(new Integer[] { 1, 3 });
		edges.add(new Integer[] { 4, 5 });
		ArrayList<Double> edgeCosts = new ArrayList<Double>();
		ArrayList<Double> prizePi = new ArrayList<Double>();
		edgeCosts.add(1.0D);
		edgeCosts.add(2.0D);
		edgeCosts.add(3.0D);
		prizePi.add(1.1D);
		prizePi.add(1.2D);
		prizePi.add(1.3D);
		prizePi.add(1.4D);
		prizePi.add(1.5D);
		F f = new F(nodes, edges, edgeCosts, prizePi, 0.0D);
		System.out.println(f.gamma);
	}

	/**
	 * test the tour algorihtm
	 */
	public void testTOUR() {
		ArrayList<Integer> nodesInT = new ArrayList<Integer>();
		nodesInT.add(12);
		nodesInT.add(13);
		nodesInT.add(14);
		nodesInT.add(21);
		nodesInT.add(129);
		nodesInT.add(30);
		ArrayList<Integer[]> edges = new ArrayList<Integer[]>();
		edges.add(new Integer[] { 12, 13 });
		edges.add(new Integer[] { 12, 14 });
		edges.add(new Integer[] { 30, 14 });
		edges.add(new Integer[] { 21, 14 });
		edges.add(new Integer[] { 129, 14 });
		ArrayList<Double> edgesCost = new ArrayList<Double>();
		edgesCost.add(1.0);
		edgesCost.add(2.0);
		edgesCost.add(3.0);
		edgesCost.add(4.0);
		edgesCost.add(5.0);
		ArrayList<Double> prizePi = new ArrayList<Double>();
		prizePi.add(1.0);
		prizePi.add(2.0);
		prizePi.add(3.0);
		prizePi.add(4.0);
		prizePi.add(5.0);
		prizePi.add(6.0);
		Tree tree = new Tree(nodesInT, edges, edgesCost, prizePi);
		for (int i : tree.getEulerTour()) {
			System.out.print(i + " ");
		}
		System.out.println();
	}

	public static void main(String args[]) {
		HeadApprox head = new HeadApprox();
		head.testTOUR();
		head.testGammaDisjointSet();
	}
}
