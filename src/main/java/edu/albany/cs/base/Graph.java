package edu.albany.cs.base;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * @author baojian
 *
 */
public abstract class Graph {

	/** node size of this graph */
	public final int[] V;
	public final int numOfEdges;
	public final int numOfNodes;
	public final ArrayList<ArrayList<Integer>> arrayListAdj;
	public final int[][] integerAdj;
	public final int[] trueSubGraph;
	public final ArrayList<Integer[]> edges;
	public final ArrayList<Double> edgeCosts;

	public Graph(APDMInputFormat apdm) {
		integerAdj = apdm.data.graphAdj;
		arrayListAdj = apdm.data.graphAdjList;
		numOfEdges = apdm.data.numEdges;
		numOfNodes = apdm.data.numNodes;
		edges = apdm.data.intEdges;
		edgeCosts = apdm.data.identityEdgeCosts;

		int[] v = new int[numOfNodes];
		for (int i = 0; i < numOfNodes; i++) {
			v[i] = i;
		}
		V = v;
		trueSubGraph = apdm.data.trueSubGraphNodes;
	}

	public Graph(int n) {
		numOfNodes = n;
		int[] v = new int[numOfNodes];
		for (int i = 0; i < numOfNodes; i++) {
			v[i] = i;
		}
		V = v;
		integerAdj = null;
		arrayListAdj = null;
		numOfEdges = 0;
		edges = null;
		edgeCosts = null;
		trueSubGraph = null;
	}

	@Override
	public String toString() {
		return "Graph [V=" + Arrays.toString(V) + ", numOfEdges=" + numOfEdges + ", numOfNodes=" + numOfNodes
				+ ", arrayListAdj=" + arrayListAdj + ", integerAdj=" + Arrays.toString(integerAdj) + ", trueSubGraph="
				+ Arrays.toString(trueSubGraph) + ", edges=" + edges + ", edgeCosts=" + edgeCosts + "]";
	}
	
}
