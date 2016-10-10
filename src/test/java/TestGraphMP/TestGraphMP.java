package TestGraphMP;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

import org.junit.Test;

import edu.albany.cs.base.APDMInputFormat;
import edu.albany.cs.base.PreRec;
import edu.albany.cs.graphMP.GraphMP;
import edu.albany.cs.scoreFuncs.EBPStat;
import edu.albany.cs.scoreFuncs.EMSStat;

public class TestGraphMP {

	private int verboseLevel = 0;

	public TestGraphMP() {
	}

	/**
	 * test EMSStat score function
	 * 
	 * @param inputFilePath
	 *            the file path of input data
	 */
	public void testEMSStatOnGrid(String inputFilePath) {

		System.out.println("\n------------------------------ test starts ------------------------------");
		System.out.println("testing file path: " + inputFilePath);
		/** step0: data file */
		APDMInputFormat apdm = new APDMInputFormat(inputFilePath);
		ArrayList<Integer[]> edges = apdm.data.intEdges;
		ArrayList<Double> edgeCosts = apdm.data.identityEdgeCosts;
		/** step1: score function */
		EMSStat func = new EMSStat(apdm.data.base, apdm.data.counts);
		/** step2: optimization */
		int[] candidateS = new int[] { 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20 };
		double optimalVal = -Double.MAX_VALUE;
		PreRec bestPreRec = new PreRec();
		GraphMP bestGraphMP = null;
		for (int s : candidateS) {
			double B = s - 1 + 0.0D;
			int t = 5;
			GraphMP graphMP = new GraphMP(edges, edgeCosts, apdm.data.base, s, 1, B, t, true/** maximumCC */
			, null, func, null);
			double[] yx = graphMP.x;
			if (func.getFuncValue(yx) > optimalVal) {
				optimalVal = func.getFuncValue(yx);
				bestPreRec = new PreRec(graphMP.resultNodes_supportX, apdm.data.trueSubGraphNodes);
				bestGraphMP = graphMP;
				if (verboseLevel == 0) {
					System.out.println("current best [pre,rec]: " + "[" + bestPreRec.pre + "," + bestPreRec.rec + "]");
				}
			}
		}
		System.out.println("precision : " + bestPreRec.pre + " ; recall : " + bestPreRec.rec);
		System.out.println("result subgraph is: " + Arrays.toString(bestGraphMP.resultNodes_Tail));
		System.out.println("true subgraph is: " + Arrays.toString(apdm.data.trueSubGraphNodes));
		System.out.println("------------------------------ test ends --------------------------------\n");
	}

	/**
	 * Test EBPStat score function
	 * 
	 * @param inputFilePath
	 *            the file path of input data
	 */
	public void testEBPStatOnGrid(String inputFilePath) {
		System.out.println("\n------------------------------ test starts ------------------------------");
		System.out.println("testing file path: " + inputFilePath);
		/** step0: data file */
		APDMInputFormat apdm = new APDMInputFormat(inputFilePath);
		ArrayList<Integer[]> edges = apdm.data.intEdges;
		ArrayList<Double> edgeCosts = apdm.data.identityEdgeCosts;
		/** step1: score function */
		EBPStat func = new EBPStat(apdm.data.base, apdm.data.counts);
		/** step2: optimization */
		int[] candidateS = new int[] { 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20 };
		double optimalVal = -Double.MAX_VALUE;
		PreRec bestPreRec = new PreRec();
		GraphMP bestGraphMP = null;
		for (int s : candidateS) {
			double B = s - 1 + 0.0D;
			int t = 5;
			GraphMP graphMP = new GraphMP(edges, edgeCosts, apdm.data.base, s, 1, B, t, true/** maximumCC */
			, null, func, null);
			double[] yx = graphMP.x;
			if (func.getFuncValue(yx) > optimalVal) {
				optimalVal = func.getFuncValue(yx);
				bestPreRec = new PreRec(graphMP.resultNodes_supportX, apdm.data.trueSubGraphNodes);
				bestGraphMP = graphMP;
				if (verboseLevel == 0) {
					System.out.println("current best [pre,rec]: " + "[" + bestPreRec.pre + "," + bestPreRec.rec + "]");
				}
			}
		}
		System.out.println("precision : " + bestPreRec.pre + " ; recall : " + bestPreRec.rec);
		System.out.println("result subgraph is: " + Arrays.toString(bestGraphMP.resultNodes_Tail));
		System.out.println("true subgraph is: " + Arrays.toString(apdm.data.trueSubGraphNodes));
		System.out.println("------------------------------ test ends --------------------------------\n");
	}

	/**
	 * To test all of the cases
	 */
	@Test
	public void testAllCases() {
		TestGraphMP testGraphMP = new TestGraphMP();
		for (File file : new File("data/GridDataEMS").listFiles()) {
			testGraphMP.testEMSStatOnGrid(file.getAbsolutePath());
		}
		for (File file : new File("data/GridDataEBP").listFiles()) {
			testGraphMP.testEBPStatOnGrid(file.getAbsolutePath());
		}
	}

	public static void main(String args[]) {
		/** just for debugging single test, your code goes here */
		new TestGraphMP().testEBPStatOnGrid("data/GridDataEBP/APDM-GridData-100_noise_0.0_trueSubSize_30_0.txt");
		new TestGraphMP().testEBPStatOnGrid("data/GridDataEBP/APDM-GridData-100_noise_0.0_trueSubSize_30_1.txt");
		new TestGraphMP().testEBPStatOnGrid("data/GridDataEBP/APDM-GridData-100_noise_0.0_trueSubSize_30_2.txt");
		new TestGraphMP().testEBPStatOnGrid("data/GridDataEBP/APDM-GridData-100_noise_0.0_trueSubSize_30_3.txt");
		new TestGraphMP().testEBPStatOnGrid("data/GridDataEBP/APDM-GridData-100_noise_0.0_trueSubSize_30_4.txt");
	}

}
