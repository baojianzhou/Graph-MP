package edu.albany.cs.scoreFuncs;

import java.util.Arrays;

public final class Stat {

	private final double[] data;
	private final int size;

	public Stat(double[] data) {
		if (data == null || data.length == 0) {
			this.data = null;
			this.size = 0;
		} else {
			this.data = new double[data.length];
			for (int i = 0; i < data.length; i++) {
				this.data[i] = data[i];
			}
			size = data.length;
		}
	}

	public double mean() {
		double sum = 0.0;
		if (data == null || data.length == 0) {
			return sum;
		}
		for (double a : data) {
			sum += a;
		}
		return sum / size;
	}

	private double mean(double[] x) {
		double sum = 0.0;
		for (double a : x){
			sum += a;
		}
			
		return sum / (x.length + 0.0D);
	}

	double getVariance() {
		double mean = mean();
		double temp = 0;
		for (double a : data)
			temp += (mean - a) * (mean - a);
		return temp / size;
	}

	double getStdDev() {
		return Math.sqrt(getVariance());
	}

	public double median() {
		Arrays.sort(data);
		if (data.length % 2 == 0) {
			return (data[(data.length / 2) - 1] + data[data.length / 2]) / 2.0;
		} else {
			return data[data.length / 2];
		}
	}

	public double mad() {
		double[] abs = new double[data.length];
		for (int i = 0; i < data.length; i++) {
			abs[i] = Math.abs(data[i] - mean(data));
		}
		return mean(abs);
	}
}