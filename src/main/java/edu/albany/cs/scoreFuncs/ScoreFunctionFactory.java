package edu.albany.cs.scoreFuncs;

public class ScoreFunctionFactory {

    public static Function getFunction(FuncType funcID, double[] b, double[] c, double[] sigma) {

        if (funcID == null) {
            return null;
        } else if (funcID.equals(FuncType.EMSStat)) {
            return new EMSStat(b, c);
        } else if (funcID.equals(FuncType.KulldorffStat)) {
            return new KulldorffStat(b, c);
        } else if (funcID.equals(FuncType.EBPStat)) {
            return new EBPStat(b, c);
        } else {
            System.out.println("Unknown Type ...");
            return null;
        }
    }

}
