package edu.albany.cs.fastPCST;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Utils used by Fast PCST aglorithm
 *
 * @param <Type>
 * @author baojian bzhou6@albany.edu
 */
public class Utils<Type> {

    //resize the array list

    public static ArrayList<EdgeInfo> resize(ArrayList<EdgeInfo> arr, int size, EdgeInfo val) {

        if (arr == null || arr.equals(null)) {
            arr = new ArrayList<EdgeInfo>();
        }

        if (size < 0) {
            new IllegalArgumentException("size should be larger than 0");
            System.exit(0);
            return null;
        } else if (size == 0) {
            return new ArrayList<EdgeInfo>();
        }

        int newSize = size - arr.size();
        if (newSize < 0) {
            for (int i = 0; i < Math.abs(newSize); i++) {
                arr.remove(arr.size() - 1);
            }
            return arr;
        } else if (newSize == 0) {
            return arr;
        } else if (newSize > 0) {
            for (int i = 0; i < newSize; i++) {
                arr.add(new EdgeInfo());
            }
            return arr;
        }
        return arr;
    }

    public static ArrayList<Integer> resize(ArrayList<Integer> arr, int size, Integer val) {

        if (arr == null || arr.equals(null)) {
            arr = new ArrayList<Integer>();
        }

        if (size < 0) {
            new IllegalArgumentException("size should be larger than 0");
            System.exit(0);
            return null;
        } else if (size == 0) {
            return new ArrayList<Integer>();
        }

        int newSize = size - arr.size();
        if (newSize < 0) {
            for (int i = 0; i < Math.abs(newSize); i++) {
                arr.remove(arr.size() - 1);
            }
            return arr;
        } else if (newSize == 0) {
            return arr;
        } else if (newSize > 0) {
            for (int i = 0; i < newSize; i++) {
                arr.add(val);
            }
            return arr;
        }
        return arr;
    }

    public static ArrayList<Boolean> resize(ArrayList<Boolean> arr, int size, Boolean val) {

        if (arr == null || arr.equals(null)) {
            arr = new ArrayList<Boolean>();
        }

        if (size < 0) {
            new IllegalArgumentException("size should be larger than 0");
            System.exit(0);
            return null;
        } else if (size == 0) {
            return new ArrayList<Boolean>();
        }

        int newSize = size - arr.size();
        if (newSize < 0) {
            for (int i = 0; i < Math.abs(newSize); i++) {
                arr.remove(arr.size() - 1);
            }
            return arr;
        } else if (newSize == 0) {
            return arr;
        } else if (newSize > 0) {
            for (int i = 0; i < newSize; i++) {
                arr.add(val);
            }
            return arr;
        }
        return arr;
    }


    public static ArrayList<Pair<Integer, Double>> resize(
            ArrayList<Pair<Integer, Double>> arr, int size, Pair<Integer, Double> pair) {

        if (arr == null || arr.equals(null)) {
            arr = new ArrayList<Pair<Integer, Double>>();
        }

        if (size < 0) {
            new IllegalArgumentException("size should be larger than 0");
            System.exit(0);
            return null;
        } else if (size == 0) {
            return new ArrayList<Pair<Integer, Double>>();
        }

        int newSize = size - arr.size();
        if (newSize < 0) {
            for (int i = 0; i < Math.abs(newSize); i++) {
                arr.remove(arr.size() - 1);
            }
            return arr;
        } else if (newSize == 0) {
            return arr;
        } else if (newSize > 0) {
            for (int i = 0; i < newSize; i++) {
                arr.add(pair);
            }
            return arr;
        }
        return arr;
    }

    public static ArrayList<ArrayList<Pair<Integer, Double>>> resize(
            ArrayList<ArrayList<Pair<Integer, Double>>> arr,
            int size, ArrayList<Pair<Integer, Double>> pair) {

        if (arr == null || arr.equals(null)) {
            arr = new ArrayList<ArrayList<Pair<Integer, Double>>>();
        }

        if (size < 0) {
            new IllegalArgumentException("size should be larger than 0");
            System.exit(0);
            return null;
        } else if (size == 0) {
            return new ArrayList<ArrayList<Pair<Integer, Double>>>();
        }

        int newSize = size - arr.size();
        if (newSize < 0) {
            for (int i = 0; i < Math.abs(newSize); i++) {
                arr.remove(arr.size() - 1);
            }
            return arr;
        } else if (newSize == 0) {
            return arr;
        } else if (newSize > 0) {
            for (int i = 0; i < newSize; i++) {
                arr.add(new ArrayList<Pair<Integer, Double>>());
            }
            return arr;
        }
        return arr;
    }

    public static ArrayList<Double> resize(ArrayList<Double> arr, int size, Double val) {

        if (arr == null || arr.equals(null)) {
            arr = new ArrayList<Double>();
        }

        if (size < 0) {
            new IllegalArgumentException("size should be larger than 0");
            System.exit(0);
            return null;
        } else if (size == 0) {
            return new ArrayList<Double>();
        }

        int newSize = size - arr.size();
        if (newSize < 0) {
            for (int i = 0; i < Math.abs(newSize); i++) {
                arr.remove(arr.size() - 1);
            }
            return arr;
        } else if (newSize == 0) {
            return arr;
        } else if (newSize > 0) {
            for (int i = 0; i < newSize; i++) {
                arr.add(val);
            }
            return arr;
        }
        return arr;
    }

    public static ArrayList<EdgePart> resize(ArrayList<EdgePart> arr, int size, EdgePart val) {

        if (arr == null || arr.equals(null)) {
            arr = new ArrayList<EdgePart>();
        }

        if (size < 0) {
            new IllegalArgumentException("size should be larger than 0");
            System.exit(0);
            return null;
        } else if (size == 0) {
            return new ArrayList<EdgePart>();
        }

        int newSize = size - arr.size();
        if (newSize < 0) {
            for (int i = 0; i < Math.abs(newSize); i++) {
                arr.remove(arr.size() - 1);
            }
            return arr;
        } else if (newSize == 0) {
            return arr;
        } else if (newSize > 0) {
            for (int i = 0; i < newSize; i++) {
                arr.add(new EdgePart());
            }
            return arr;
        }
        return arr;
    }

    public static ArrayList<EdgePart> resize(ArrayList<EdgePart> arr, int size, boolean flag) {

        if (arr == null || arr.equals(null)) {
            arr = new ArrayList<EdgePart>();
        }

        if (size < 0) {
            new IllegalArgumentException("size should be larger than 0");
            System.exit(0);
            return null;
        } else if (size == 0) {
            return new ArrayList<EdgePart>();
        }

        int newSize = size - arr.size();
        if (newSize < 0) {
            for (int i = 0; i < Math.abs(newSize); i++) {
                arr.remove(arr.size() - 1);
            }
            return arr;
        } else if (newSize == 0) {
            return arr;
        } else if (newSize > 0) {
            for (int i = 0; i < newSize; i++) {
                arr.add(new EdgePart());
            }
            return arr;
        }
        return arr;
    }

    public static void printIntegerArrayList(ArrayList<Integer[]> arr) {
        for (int i = 1; i <= arr.size(); i++) {
            System.out.print(Arrays.toString(arr.get(i - 1)) + "  ");
            if (i % 10 == 0) {
                System.out.println();
            }
        }
    }

    /**
     * stop for testing
     */
    public static void stop() {
        try {
            System.out.println("Press any key to continue ...");
            System.in.read();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    public static void main(String args[]) {
        ArrayList<Integer> arr = new ArrayList<Integer>(4);
        arr.add(1);
        arr.add(2);
        arr.add(4);
        arr.add(5);
        System.out.println(arr.toString());
        arr = Utils.resize(arr, 2, 10);
        System.out.println(arr.toString());
    }
}
