package edu.albany.cs.base;

import org.apache.commons.lang3.ArrayUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * In this class, using APMDIOFormat class/object to get the data from STP-like
 * file. More details can access DIMACS
 * [http://dimacs11.cs.princeton.edu/downloads.html] and STP File Format Page
 * [http://steinlib.zib.de/format.php]
 *
 * @author Baojian Zhou
 * @see <a href="http://steinlib.zib.de/format.php">STP File Format</a>
 * @see <a href="http://dimacs11.cs.princeton.edu/downloads.html">DIMACS</a>
 */
public class APDMInputFormat {

	// header
	private final static String inputHeader = Utils.commentLine + "\n"
			+ "#APDM Input Graph, this input graph includes 3 sections:\n" + "#section1 : general information\n"
			+ "#section2 : nodes\n" + "#section3 : edges\n" + "#section4 : trueSubGraph (Optional)\n" + "#\n"
			+ "#if nodes haven't information set weight to null\n" + Utils.commentLine + "\n";
	public InputData data;

	/**
	 * read APDM input file
	 *
	 * @param APDMFile
	 *            APDM input file
	 */
	public APDMInputFormat(File APDMFile) {

		/** read inputFile */
		readAPDMFile(APDMFile);
		/** initialize true subgraph nodes */
		data.trueSubGraphNodes = null;
		if (data.trueSubGraphEdges != null) {
			int[] trueSubGraphNodes = null;
			for (int[] edge : data.trueSubGraphEdges.keySet()) {
				if (!ArrayUtils.contains(trueSubGraphNodes, edge[0])) {
					trueSubGraphNodes = ArrayUtils.add(trueSubGraphNodes, edge[0]);
				}
				if (!ArrayUtils.contains(trueSubGraphNodes, edge[1])) {
					trueSubGraphNodes = ArrayUtils.add(trueSubGraphNodes, edge[1]);
				}
			}
			data.trueSubGraphNodes = trueSubGraphNodes;
		}
		/** initialize adjacency information */
		ArrayList<ArrayList<Integer>> graphAdjList = new ArrayList<>();
		ArrayList<ArrayList<Double>> graphWeightedAdjList = new ArrayList<>();
		int[][] graphAdj = new int[data.numNodes][];
		double[][] graphWeightedAdj = new double[data.numNodes][];
		for (int i = 0; i < data.numNodes; i++) {
			graphAdjList.add(new ArrayList<>());
			graphWeightedAdjList.add(new ArrayList<>());
		}
		for (int[] edge : data.edges.keySet()) {
			if (!graphAdjList.get(edge[0]).contains(edge[1])) {
				graphAdjList.get(edge[0]).add(edge[1]);
				graphWeightedAdjList.get(edge[0]).add(data.edges.get(edge));
				graphAdj[edge[0]] = ArrayUtils.add(graphAdj[edge[0]], edge[1]);
				graphWeightedAdj[edge[0]] = ArrayUtils.add(graphWeightedAdj[edge[0]], data.edges.get(edge));
			}
			if (!graphAdjList.get(edge[1]).contains(edge[0])) {
				graphAdjList.get(edge[1]).add(edge[0]);
				graphWeightedAdjList.get(edge[1]).add(data.edges.get(edge));
				graphAdj[edge[1]] = ArrayUtils.add(graphAdj[edge[1]], edge[0]);
				graphWeightedAdj[edge[1]] = ArrayUtils.add(graphWeightedAdj[edge[1]], data.edges.get(edge));
			}
		}
		data.graphAdj = graphAdj;
		data.graphAdjList = graphAdjList;
		data.graphWeightedAdj = graphWeightedAdj;
		data.graphWeightedAdjList = graphWeightedAdjList;
		/** other parameters. */
		data.V = new int[data.numNodes];
		for (int i = 0; i < data.V.length; i++) {
			data.V[i] = i;
		}
		data.cc = new ConnectedComponents(data.graphAdjList);
	}

	public APDMInputFormat(String apdmFileName) {
		this(new File(apdmFileName));/** read inputFile */
	}

	/**
	 * According to the APDM input format, we generate APDM input file
	 *
	 * @param usedAlgorithm
	 *            this kind of input file is used for usedAlgorithm
	 * @param dataSource
	 *            dataSource [WaterPollutionDataset CivilUnrestDataset
	 *            GridDataset SNAPDataset TransportationDataset ]
	 * @param edges
	 *            the edges of the graph
	 * @param PValue
	 *            the nodes of the graph
	 * @param fileName
	 *            the input file name
	 */
	public static void generateAPDMFile(String usedAlgorithm, String dataSource, ArrayList<Edge> edges, double[] PValue,
			double[] counts, double[] averValue, HashMap<int[], Double> trueSubGraphEdges, String fileName) {
		DecimalFormat decimalFormat = new DecimalFormat("0.000000");
		try {
			FileWriter fw;
			fw = new FileWriter(fileName, false);
			fw.write(APDMInputFormat.inputHeader);
			// general information
			fw.write("SECTION1 (General Information)\n");
			if (PValue == null) {
				fw.write("numNodes = " + 0 + "\n");
			} else {
				fw.write("numNodes = " + PValue.length + "\n");
			}
			if (edges == null) {
				fw.write("numEdges = " + 0 + "\n");
			} else {
				fw.write("numEdges = " + edges.size() + "\n");
			}
			fw.write("usedAlgorithm = " + usedAlgorithm + "\n");
			fw.write("dataSource = " + dataSource + "\n");
			fw.write("END\n" + Utils.commentLine + "\n");

			// nodes information
			fw.write("SECTION2 (Nodes Information)\n");
			fw.write("NodeID counts mean std\n");
			if (PValue == null) {
				fw.write("null\n");
			} else {
				if (counts == null && averValue == null) {
					for (int i = 0; i < PValue.length; i++) {
						fw.write(i + " " + decimalFormat.format(PValue[i]) + "\n");
					}
				} else {
					for (int i = 0; i < PValue.length; i++) {
						fw.write(i + " " + decimalFormat.format(PValue[i]) + " " + counts[i] + " " + averValue[i]
								+ "\n");
					}
				}

			}
			fw.write("END\n" + Utils.commentLine + "\n");

			// edges information
			fw.write("SECTION3 (Edges Information)\n");
			fw.write("EndPoint0 EndPoint1 Weight\n");
			if (edges != null) {
				for (Edge e : edges) {
					fw.write(e.i + " " + e.j + " " + decimalFormat.format(e.cost) + "\n");
				}
			}
			fw.write("END\n" + Utils.commentLine + "\n");
			// edges information
			fw.write("SECTION4 (TrueSubGraph Information)\n");
			fw.write("EndPoint0 EndPoint1 Weight\n");
			if (trueSubGraphEdges != null) {
				for (int[] e : trueSubGraphEdges.keySet()) {
					fw.write(e[0] + " " + e[1] + " " + trueSubGraphEdges.get(e) + "\n");
				}
			}
			fw.write("END\n" + Utils.commentLine + "\n");
			fw.close();
		} catch (IOException e1) {
			e1.printStackTrace();
		}
	}

	/**
	 * According to the APDM input format, we generate APDM input file
	 *
	 * @param usedAlgorithm
	 *            this kind of input file is used for usedAlgorithm
	 * @param dataSource
	 *            dataSource [WaterPollutionDataset CivilUnrestDataset
	 *            GridDataset SNAPDataset TransportationDataset ]
	 * @param edges
	 *            the edges of the graph
	 * @param PValue
	 *            the nodes of the graph
	 * @param fileName
	 *            the input file name
	 */
	public static void generateAPDMFile(String usedAlgorithm, String dataSource, ArrayList<Edge> edges, double[] PValue,
			double[] counts, double[] averValue, ArrayList<Edge> trueSubGraphEdges, String fileName) {
		DecimalFormat decimalFormat = new DecimalFormat("0.000000");

		try {
			FileWriter fw = new FileWriter(fileName, false);
			fw.write(APDMInputFormat.inputHeader);
			// general information
			fw.write("SECTION1 (General Information)\n");
			if (PValue == null) {
				fw.write("numNodes = " + 0 + "\n");
			} else {
				fw.write("numNodes = " + PValue.length + "\n");
			}
			if (edges == null) {
				fw.write("numEdges = " + 0 + "\n");
			} else {
				fw.write("numEdges = " + edges.size() / 2 + "\n");
			}
			fw.write("usedAlgorithm = " + usedAlgorithm + "\n");
			fw.write("dataSource = " + dataSource + "\n");
			fw.write("END\n" + Utils.commentLine + "\n");

			/** section2 is nodes information */
			fw.write("SECTION2 (Nodes Information)\n");
			fw.write("NodeID Base Counts\n");
			if (PValue == null) {
				fw.write("null\n");
			} else {
				if (counts == null && averValue == null) {
					for (int i = 0; i < PValue.length; i++) {
						fw.write(i + " " + decimalFormat.format(PValue[i]) + "\n");
					}
				} else {
					if (averValue == null) {
						for (int i = 0; i < PValue.length; i++) {
							fw.write(i + " " + decimalFormat.format(PValue[i]) + " " + counts[i] + "\n");
						}
					} else {
						for (int i = 0; i < PValue.length; i++) {
							fw.write(i + " " + decimalFormat.format(PValue[i]) + " " + counts[i] + " " + averValue[i]
									+ "\n");
						}
					}

				}

			}
			fw.write("END\n" + Utils.commentLine + "\n");
			/** section3 is edges information */
			fw.write("SECTION3 (Edges Information)\n");
			fw.write("EndPoint0 EndPoint1 Weight\n");
			if (edges != null) {
				for (Edge e : edges) {
					fw.write(e.i + " " + e.j + " " + decimalFormat.format(e.cost) + "\n");
				}
			}
			fw.write("END\n" + Utils.commentLine + "\n");
			/** section4 is true subgraph information */
			fw.write("SECTION4 (TrueSubGraph Information)\n");
			fw.write("EndPoint0 EndPoint1 Weight\n");
			if (trueSubGraphEdges != null) {
				for (Edge e : trueSubGraphEdges) {
					fw.write(e.i + " " + e.j + " " + e.cost + "\n");
				}
			}
			fw.write("END\n" + Utils.commentLine + "\n");
			fw.close();
		} catch (IOException e1) {
			e1.printStackTrace();
		}

	}

	/**
	 * @param APDMFile
	 *            read APDM from file
	 */
	private boolean readAPDMFile(File APDMFile) {
		BufferedReader br = null;
		data = new InputData();
		try {
			String sCurrentLine;
			br = new BufferedReader(new FileReader(APDMFile));

			while ((sCurrentLine = br.readLine()) != null) {
				if (sCurrentLine.startsWith("#")) {
					continue;
				}
				/** general information of this graph. */
				if (sCurrentLine.startsWith("SECTION1")) {
					while (!(sCurrentLine = br.readLine()).equals("END")) {
						if (sCurrentLine.startsWith("numNodes")) {
							String[] str = sCurrentLine.split(" ");
							data.numNodes = Integer.parseInt(str[2]);
						}
						if (sCurrentLine.startsWith("numEdges")) {
							String[] str = sCurrentLine.split(" ");
							int numEdges = Integer.parseInt(str[2]);
							data.numEdges = numEdges;
						}
						if (sCurrentLine.startsWith("dataSource")) {
							String[] str = sCurrentLine.trim().split(" ");
							String dataSource = str[2];
							this.data.dataSource = dataSource;
						}
						if (sCurrentLine.startsWith("usedAlgorithm")) {
							String[] str = sCurrentLine.split(" ");
							String usedAlgorithm = str[2];
							this.data.usedAlgorithm = usedAlgorithm;
						}
					}
				}

				/** nodes information */
				if (sCurrentLine.startsWith("SECTION2")) {
					data.nodes = new HashMap<Integer, Double>();
					data.counts = new double[data.numNodes];
					data.base = new double[data.numNodes];
					int count = 0;
					while (!(sCurrentLine = br.readLine()).equals("END")) {
						if (sCurrentLine.startsWith("NodeID")) {
							continue;
						}
						String[] str = sCurrentLine.split(" ");
						data.nodes.put(Integer.parseInt(str[0]), Double.parseDouble(str[1]));
						if (data.dataSource.equals("GridData")) {
							data.base[count] = Double.parseDouble(str[1]);
							data.counts[count] = Double.parseDouble(str[2]);
						}
						count++;
					}
				}
				/** edges information */
				if (sCurrentLine.startsWith("SECTION3")) {
					data.edges = new HashMap<int[], Double>();
					data.newEdges = new ArrayList<Edge>();
					data.intEdges = new ArrayList<Integer[]>();
					data.edgeCosts = new ArrayList<Double>();
					int count = 0;
					while (!(sCurrentLine = br.readLine()).equals("END")) {
						if (sCurrentLine.startsWith("EndPoint0")) {
							continue;
						}
						String[] str = sCurrentLine.split(" ");
						int[] edge = new int[] { Integer.parseInt(str[0]), Integer.parseInt(str[1]) };
						data.edges.put(edge, Double.parseDouble(str[2]));
						data.newEdges.add(new Edge(edge[0], edge[1], count++, Double.parseDouble(str[2])));
						data.intEdges.add(new Integer[] { edge[0], edge[1] });
						data.edgeCosts.add(Double.parseDouble(str[2]));
					}
					ArrayList<Double> identityEdgeCosts = new ArrayList<Double>();
					for (int i = 0; i < data.edgeCosts.size(); i++) {
						identityEdgeCosts.add(1.0D);
					}
					data.identityEdgeCosts = identityEdgeCosts;
				}
				/** true subgraph edges information */
				if (sCurrentLine.startsWith("SECTION4")) {
					data.trueSubGraphEdges = new HashMap<int[], Double>();
					while (!(sCurrentLine = br.readLine()).equals("END")) {
						if (sCurrentLine.startsWith("EndPoint0") || sCurrentLine.startsWith("null")) {
							continue;
						}
						String[] str = sCurrentLine.split(" ");
						int[] edge = new int[] { Integer.parseInt(str[0]), Integer.parseInt(str[1]) };
						data.trueSubGraphEdges.put(edge, Double.parseDouble(str[2]));
					}
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			try {
				if (br != null)
					br.close();
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		}
		return true;
	}

	public class InputData {

		public String dataSource;
		public String usedAlgorithm;

		public ArrayList<ArrayList<Integer>> graphAdjList;
		public int[][] graphAdj;
		public ArrayList<ArrayList<Double>> graphWeightedAdjList;
		public double[][] graphWeightedAdj;
		public boolean connective;
		public ConnectedComponents cc;

		public int numNodes;
		public int numEdges;
		public HashMap<Integer, Double> nodes;
		public int[] V;
		public double[] base;
		public double[] counts;

		public HashMap<int[], Double> edges;
		public ArrayList<Integer[]> intEdges;
		public ArrayList<Double> edgeCosts;
		public ArrayList<Double> identityEdgeCosts;
		public ArrayList<Edge> newEdges;

		public int[] trueSubGraphNodes;
		public HashMap<int[], Double> trueSubGraphEdges = null;
	}

	public static void main(String args[]) throws IOException {
		/** your test code goes here. */
	}

}
