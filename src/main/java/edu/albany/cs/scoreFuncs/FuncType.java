package edu.albany.cs.scoreFuncs;

/**
 * Types of different score functions.
 *
 * @author baojian bzhou6@albany.edu
 */
public enum FuncType {

	KulldorffStat, EBPStat, EMSStat, Unknown;

	public static FuncType defaultFuncType() {
		return FuncType.Unknown;
	}
}
