package edu.albany.cs.scoreFuncs;

import edu.albany.cs.base.ArrayIndexSort;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.stat.StatUtils;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * Expectation-based Poisson Statistic. More details can be found in our paper.
 * Please Table 1 of the paper, Fast subset scan for spatial pattern detection.
 * 
 * @author baojian bzhou6@albany.edu
 * 
 * @see <a href="http://www.ijcai.org/Proceedings/16/Papers/200.pdf">Our paper
 *      link</a>
 * @see <a href="http://www.cs.cmu.edu/~./neill/papers/jrssb2012.pdf"> Others
 *      </a>
 */
public class EBPStat implements Function {

	/** base */
	private final double[] b;
	/** counts of application data */
	private final double[] c;
	/** vector size */
	private final int n;

	private int verboseLevel = 0;
	private final FuncType funcID;

	public EBPStat(double[] b, double[] c) {

		this.funcID = FuncType.EBPStat;
		if (!checkInput(b, c)) {
			System.out.println(funcID + " input parameter is invalid.");
			System.exit(0);
		}
		this.b = b;
		this.c = c;
		this.n = b.length;

	}

	/**
	 * make sure that the inputs of EBP are valid
	 * 
	 * @return true if it is valid; otherwise, return false.
	 */
	private boolean checkInput(double[] b, double[] c) {

		if (verboseLevel > 0) {
			System.out.println("BAll: " + StatUtils.sum(b));
			System.out.println("CAll: " + StatUtils.sum(c));
		}
		if (b == null || c == null || b.length == 0 || c.length == 0) {
			return false;
		} else {
			if (StatUtils.sum(b) <= 0.0D || StatUtils.sum(c) <= 0.0D) {
				return false;
			}
			for (int i = 0; i < n; i++) {
				if (b[i] <= 0.0D) {
					return false;
				}
				if (c[i] < 0.0D) {
					return false;
				}
			}
			return true;
		}
	}

	/**
	 * get gradient of EBP
	 */
	@Override
	public double[] getGradient(double[] x) {

		checkIndictorVect(x);
		if (StatUtils.sum(x) <= 0.0D) {
			x[new Random().nextInt(n)] = 1.0D;
		}
		double[] gradient = new double[n];
		double B = new ArrayRealVector(x).dotProduct(new ArrayRealVector(b));
		double C = new ArrayRealVector(x).dotProduct(new ArrayRealVector(c));
		if (verboseLevel > 0) {
			System.out.println("B is : " + B + " ; C is : " + C);
		}
		if (B <= 0.0D || C <= 0.0D) {
			System.out.println("Poission Gradient Error : B or C is non-positive value ...");
			System.out.println("B is: " + B + " C is: " + C);
			System.exit(0);
		}
		if (C > B) {
			for (int i = 0; i < n; i++) {
				gradient[i] = c[i] * Math.log(C / B) + b[i] * (1 - C / B);
				if (!Double.isFinite(gradient[i])) {
					System.out.println("gradient error ...");
					System.exit(0);
				}
			}
		} else {
			Arrays.fill(gradient, 0.0D);
		}
		return gradient;
	}

	@Override
	public double getFuncValue(double[] x) {

		checkIndictorVect(x);
		double B = new ArrayRealVector(x).dotProduct(new ArrayRealVector(b));
		double C = new ArrayRealVector(x).dotProduct(new ArrayRealVector(c));
		double f = 0.0D;
		/** B <= C, the function value is 0.0 */
		if (C > B) {
			f = C * Math.log(C / B) + B - C;
			if (!Double.isFinite(f)) {
				System.out.println("EBP stat getFuncValue error. it is not finite. ");
				System.out.println("B: " + B + " C: " + C);
				System.exit(0);
			}
		}
		return f;
	}

	/**
	 * make sure x is an indicator vector
	 * 
	 * @param x
	 *            input vector
	 */
	private void checkIndictorVect(double[] x) {

		if (x == null || x.length != n) {
			System.out.println("Kulldorff gradient error : Invalid parameters ...");
			System.exit(0);
		}

		for (int i = 0; i < n; i++) {
			if (x[i] < 0.0D || x[i] > 1.0D) {
				System.out.println("x[i] should be in [0,1], but it is " + x[i]);
				System.exit(0);
			}
		}
	}

	@Override
	public FuncType getFuncID() {
		return funcID;
	}

	@Override
	public double[] getGradient(int[] S) {

		double[] x = new double[n];
		Arrays.fill(x, 0.0D);
		for (int i : S) {
			x[i] = 1.0D;
		}
		return getGradient(x);
	}

	@Override
	public double getFuncValue(int[] S) {

		double[] x = new double[n];
		Arrays.fill(x, 0.0D);
		for (int i : S) {
			x[i] = 1.0D;
		}
		return getFuncValue(x);
	}

	@Override
	public BigDecimal[] getGradientBigDecimal(BigDecimal[] x) {

		double[] xD = new double[n];
		for (int i = 0; i < n; i++) {
			xD[i] = x[i].doubleValue();
		}
		double[] gradient = getGradient(xD);
		BigDecimal[] grad = new BigDecimal[n];
		for (int i = 0; i < n; i++) {
			grad[i] = new BigDecimal(gradient[i]);
		}
		return grad;
	}

	@Override
	public double[] getArgMinFx(ArrayList<Integer> S) {
		/**
		 * as our objective function is f(S), we do min_{S} -f(S), this is
		 * equivalent to maximize f(S)
		 */
		return getArgMaxFx(S);
	}
	
	/**
	 * This is a heuristic method. It maximizes objective function f(S).
	 * 
	 * @param S
	 *            the constraint set S
	 * @return the maximizer of objective function
	 */
	public double[] getArgMaxFx(ArrayList<Integer> S) {

		double[] result = new double[n];
		Double[] vectorRatioCB = new Double[S.size()];
		for (int i = 0; i < S.size(); i++) {
			vectorRatioCB[i] = c[S.get(i)] / b[S.get(i)];
		}
		ArrayIndexSort arrayIndexComparator = new ArrayIndexSort(vectorRatioCB);
		Integer[] indexes = arrayIndexComparator.getIndices();
		Arrays.sort(indexes, arrayIndexComparator);
		/** v_1,v_2,...,v_m from large to small */
		ArrayList<Integer> sortedS = new ArrayList<Integer>();
		/** indexes from large to small */
		for (int index : indexes) {
			sortedS.add(S.get(index));
		}
		double maxF = -Double.MAX_VALUE;
		double[] argMaxX = null;
		for (int k = 1; k <= sortedS.size(); k++) {
			List<Integer> Rk = sortedS.subList(0, k);
			double[] x = new double[n];
			for (int i = 0; i < n; i++) {
				x[i] = 0.0D;
			}
			for (int index : Rk) {
				x[index] = 1.0D;
			}
			double fk = getFuncValue(x);
			if (fk > maxF) {
				maxF = fk;
				argMaxX = x;
			}
		}
		result = argMaxX;
		return result;
	}

	public static void main(String args[]) {
		/** TODO test cases goes here. */
	}
}
