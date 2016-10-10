package edu.albany.cs.graphMP;

import edu.albany.cs.base.Utils;
import edu.albany.cs.headApprox.HeadApprox;
import edu.albany.cs.scoreFuncs.Function;
import edu.albany.cs.tailApprox.TailApprox;

import org.apache.commons.lang3.ArrayUtils;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Random;

/**
 * Algorithm 1 : Graph-Mp Aglorithm in our IJCAI paper.
 *
 * @author Baojian bzhou6@albany.edu
 */
public class GraphMP {

	/** 1Xn dimension, input data */
	private final double[] c;
	private final int graphSize;
	/** graph info */
	private final HashSet<Integer> nodes;
	private final ArrayList<Integer[]> edges;
	private final ArrayList<Double> edgeCosts;
	/** the total sparsity of S */
	private final int s;
	/** the maximum number of connected components formed by F */
	private final int g;
	/** bound on the total weight w(F) of edges in the forest F */
	private final double B;
	/** number of iterations */
	private final int t;
	private final int[] trueSubGraph;
	private final Function function;

	/** results */
	public double[] x;
	public int[] resultNodes_supportX;
	public int[] resultNodes_Tail = null;
	public double funcValue = -1.0D;
	public double runTime;

	private int verboseLevel = 0;

	public GraphMP(ArrayList<Integer[]> edges, ArrayList<Double> edgeCosts, double[] c, int s, int g, double B, int t,
			boolean singleNodeInitial, int[] trueSubGraph, Function func, String resultFileName) {
		this(edges, edgeCosts, c, s, g, B, t, trueSubGraph, func, resultFileName, null);
	}

	public GraphMP(ArrayList<Integer[]> edges, ArrayList<Double> edgeCosts, double[] c, int s, int g, double B, int t,
			int[] trueSubGraph, Function func, String resultFileName, String fileName) {
		this.edges = edges;
		this.nodes = new HashSet<>();
		for (Integer[] edge : edges) {
			nodes.add(edge[0]);
			nodes.add(edge[1]);
		}
		this.graphSize = nodes.size();
		this.edgeCosts = edgeCosts;
		this.c = c;
		this.s = s;
		this.g = g;
		this.B = B;
		this.t = t;
		this.trueSubGraph = trueSubGraph;
		this.function = func;
		/** run Graph-MP algorithm. */
		x = run();
	}

	private double[] run() {

		long startTime = System.nanoTime();
		double[] x = initializeRandom();
		ArrayList<Double> fValues = new ArrayList<>();
		for (int i = 0; i < this.t; i++) { // t iterations
			if (verboseLevel > 0) {
				System.out.println("------------iteration: " + i + "------------");
			}
			fValues.add(function.getFuncValue(x));
			double[] gradientF = function.getGradient(x);
			gradientF = normalizeGradient(x, gradientF);
			/** head approximation */
			HeadApprox head = new HeadApprox(edges, edgeCosts, gradientF, s, g, B, trueSubGraph);
			ArrayList<Integer> S = Utils.unionSets(head.bestForest.nodesInF, support(x));
			/** tail approximation */
			double[] b = function.getArgMinFx(S);
			TailApprox tail = new TailApprox(edges, edgeCosts, b, s, g, B, trueSubGraph);
			/** calculate x^{i+1} */
			for (int j = 0; j < b.length; j++) {
				x[j] = 0.0D;
			}
			for (int j : tail.bestForest.nodesInF) {
				x[j] = b[j];
			}
			if (verboseLevel > 0) {
				System.out.println("number of head nodes : " + head.bestForest.nodesInF.size());
				System.out.println("number of tail nodes : " + tail.bestForest.nodesInF.size());
			}
			resultNodes_Tail = Utils.getIntArrayFromIntegerList(tail.bestForest.nodesInF);
		}
		resultNodes_supportX = getSupportNodes(x);
		funcValue = function.getFuncValue(x);
		runTime = (System.nanoTime() - startTime) / 1e9;
		return x;
	}

	/**
	 * in order to fit the gradient, we need to normalize the gradient if the
	 * domain of function is within [0,1]^n
	 * 
	 * @param x
	 *            input vector x
	 * @param gradient
	 *            gradient vector
	 * @return normGradient the normalized gradient
	 */
	private double[] normalizeGradient(double[] x, double[] gradient) {
		double[] normalizedGradient = new double[graphSize];
		for (int i = 0; i < graphSize; i++) {
			if ((gradient[i] < 0.0D) && (x[i] == 0.0D)) {
				normalizedGradient[i] = 0.0D;
			} else if ((gradient[i] > 0.0D) && (x[i] == 1.0D)) {
				normalizedGradient[i] = 0.0D;
			} else {
				normalizedGradient[i] = gradient[i];
			}
		}
		return normalizedGradient;
	}

	private double[] initializeRandom() {
		double[] x0 = new double[c.length];
		Random rand = new Random();
		for (int i = 0; i < c.length; i++) {
			if (rand.nextDouble() < 0.5D) {
				x0[i] = 1.0D;
			} else {
				x0[i] = 0.0D;
			}
		}
		return x0;
	}

	/**
	 * get a support of a vector
	 *
	 * @param x
	 *            array x
	 * @return a subset of nodes corresponding the index of vector x with
	 *         entries not equal to zero
	 */
	public ArrayList<Integer> support(double[] x) {
		if (x == null) {
			return null;
		}
		ArrayList<Integer> nodes = new ArrayList<>();
		for (int i = 0; i < x.length; i++) {
			if (x[i] != 0.0D) {
				nodes.add(i);
			}
		}
		return nodes;
	}

	/**
	 * get the nodes returned by algorithm
	 *
	 * @param x
	 *            array x
	 * @return the result nodes
	 */
	private int[] getSupportNodes(double[] x) {
		int[] nodes = null;
		for (int i = 0; i < x.length; i++) {
			if (x[i] != 0.0D) {
				/** get nonzero nodes */
				nodes = ArrayUtils.add(nodes, i);
			}
		}
		return nodes;
	}
}
