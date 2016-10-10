package edu.albany.cs.base;

import org.apache.commons.lang3.ArrayUtils;

import java.util.Arrays;
import java.util.Comparator;

public class ArrayIndexSort implements Comparator<Integer> {

	/** we want to sort array by descending order */
	private final boolean isDescend;
	/** data */
	private final Double[] data;
	/** keep indexes in descending order */
	private final Integer[] indexes;

	public ArrayIndexSort(Double[] dataArr) {
		this.isDescend = true;
		this.data = dataArr;
		this.indexes = createIndexArray();
	}

	public ArrayIndexSort(double[] dataArr) {
		this.isDescend = true;
		this.data = ArrayUtils.toObject(dataArr);
		this.indexes = createIndexArray();
	}

	public Integer[] createIndexArray() {
		Integer[] indexes = new Integer[data.length];
		for (int i = 0; i < this.data.length; i++) {
			indexes[i] = i;
		}
		return indexes;
	}

	public boolean isDescend() {
		return isDescend;
	}

	public Integer[] getIndices() {
		return indexes;
	}

	public Double[] getSortedData() {
		Arrays.sort(indexes, this);
		return null;
	}

	@Override
	public int compare(Integer index1, Integer index2) {
		/** descending order i1 > i2 > ... > .... */
		return data[index2].compareTo(data[index1]);
	}

	public static void main(String args[]) {
		Double[] vectorRatioCB = new Double[] { 0.1, 0.03, 0.2, 0.04, 0.5 };
		ArrayIndexSort arrayIndexComparator = new ArrayIndexSort(vectorRatioCB);
		Integer[] indexes = arrayIndexComparator.indexes;
		System.out.println(Arrays.toString(indexes));
		Arrays.sort(indexes, arrayIndexComparator);
		System.out.println(Arrays.toString(indexes));
		for (int i : indexes) {
			System.out.print(vectorRatioCB[i] + " ");
		}
	}
}