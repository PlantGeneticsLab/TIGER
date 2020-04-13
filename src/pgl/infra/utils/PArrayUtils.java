/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package pgl.infra.utils;

import java.util.*;

import org.apache.commons.lang.ArrayUtils;

/**
 *
 * @author Fei Lu 
 */
public class PArrayUtils {
    

    /**
     * Return the index bounds of (given-number) subsets of an array, bound[i][0] is inclusive, bound[i][1] is exclusive
     * @param arraySize
     * @param subsetNumber
     * @return 
     */
    public static int[][] getSubsetsIndicesBySubsetNumber (int arraySize, int subsetNumber) {
        if (arraySize < subsetNumber) return null;
        int remainder = arraySize % subsetNumber;
        int size = (arraySize-remainder)/subsetNumber;
        int[][] bound = new int[subsetNumber][2];
        int cnt = 0;
        for (int i = 0; i < remainder; i++) {
            bound[i][0] = cnt;
            cnt=cnt+size+1;
            bound[i][1] = cnt;
        }
        for (int i = remainder; i < subsetNumber; i++) {
            bound[i][0] = cnt;
            cnt+=size;
            bound[i][1] = cnt;
        }
        return bound;
    }
    
    /**
     * Return the index bounds of (given-size) subsets of an array, bound[i][0] is inclusive, bound[i][1] is exclusive
     * @param arraySize
     * @param subsetSize
     * @return 
     */
    public static int[][] getSubsetsIndicesBySubsetSize (int arraySize, int subsetSize) {
        int remainder = arraySize%subsetSize;
        int groupNumber;
        if (remainder == 0) groupNumber = arraySize/subsetSize;
        else groupNumber = arraySize/subsetSize + 1;
        int[][] bound = new int[groupNumber][2];
        for (int i = 0; i < groupNumber; i++) {
            bound[i][0] = i*subsetSize;
            bound[i][1] = i*subsetSize + subsetSize;
            if (bound[i][1] > arraySize) bound[i][1] = arraySize;
        }
        return bound;
    }
    
    /**
     * Return a random subset of an integer array, may have redundant value
     * @param ar
     * @param size
     * @return 
     */
    public static int[] getRandomSubset (int[] ar, int size) {
        if (size > ar.length) return null;
        int[] index = getRandomIntArray (ar.length, size);
        int[] nar = new int[size];
        for (int i = 0; i < nar.length; i++) nar[i] = ar[index[i]];
        return nar;
    }
    
    /**
     * Return a random subset of a double array, may have redundant value
     * @param ar
     * @param size
     * @return 
     */
    public static double[] getRandomSubset (double[] ar, int size) {
        if (size > ar.length) return null;
        int[] index = getRandomIntArray (ar.length, size);
        double[] nar = new double[size];
        for (int i = 0; i < nar.length; i++) nar[i] = ar[index[i]];
        return nar;
    }
    
    /**
     * Generate a random integer array, plus value, may have redundant value
     * Good when maxValue is large
     * Designed to generate subset indices
     * @param maxValue
     * @param size
     * @return 
     */
    public static int[] getRandomIntArray (int maxValue, int size) {
        int[] ar = new int[size];
        for (int i = 0; i < ar.length; i++) {
            ar[i] = (int)(maxValue*Math.random());
        }
        return ar;
    }
    
    /**
     * Return a random subset of an integer array, no redundant value
     * @param ar
     * @param size
     * @return 
     */
    public static int[] getNonredundantRandomSubset (int[] ar, int size) {
        if (size > ar.length) return null;
        int[] index = getNonredundantRandomIntArray (ar.length, size);
        int[] nar = new int[size];
        for (int i = 0; i < nar.length; i++) nar[i] = ar[index[i]];
        return nar;
    }
    
    /**
     * Return a random subset of a double array, no redundant value
     * @param ar
     * @param size
     * @return 
     */
    public static double[] getNonredundantRandomSubset (double[] ar, int size) {
        if (size > ar.length) return null;
        int[] index = getNonredundantRandomIntArray (ar.length, size);
        double[] nar = new double[size];
        for (int i = 0; i < nar.length; i++) nar[i] = ar[index[i]];
        return nar;
    }
    
    /**
     * Generate a random integer array, plus value
     * Good when maxValue is small
     * Designed to generate subset indices
     * @param maxValue
     * @param size
     * @return 
     */
    public static int[] getNonredundantRandomIntArray (int maxValue, int size) {
        int[] ar = new int[maxValue];
        for (int i = 0; i < maxValue; i++) ar[i] = i;
        shuffleArray(ar);
        int[] nar = new int[size];
        for (int i = 0; i < size; i++) {
            nar[i] = ar[i];
        }
        return nar;
    }
    
    /**
     * Shuffle an int array
     * Implementing Fisher–Yates shuffle
     * @param ar 
     */
    public static void shuffleArray(int[] ar) {
        Random rnd = new Random();
        for (int i = ar.length - 1; i > 0; i--) {
            int index = rnd.nextInt(i + 1);
            // Simple swap
            int a = ar[index];
            ar[index] = ar[i];
            ar[i] = a;
        }
    }
    
    /**
     * Shuffle an double array
     * Implementing Fisher–Yates shuffle
     * @param ar 
     */
    public static void shuffleArray(double[] ar) {
        Random rnd = new Random();
        for (int i = ar.length - 1; i > 0; i--) {
            int index = rnd.nextInt(i + 1);
            // Simple swap
            double a = ar[index];
            ar[index] = ar[i];
            ar[i] = a;
        }
    }

    /**
     * Return an index of an array by descending order of value
     * @param array
     * @return 
     */
    public static int[] getIndicesByDescendingValue(int[] array) {
        int[] index = getIndicesByAscendingValue(array);
        ArrayUtils.reverse(index);
        return index;
    }
    
    /**
     * Return an index of an array by ascending order of value
     * @param array
     * @return 
     */
    public static int[] getIndicesByAscendingValue(int[] array) {
        int[] inputArray = new int[array.length];
        System.arraycopy(array, 0, inputArray, 0, array.length);
        Integer[] idx = new Integer[inputArray.length];
        for( int i = 0 ; i < idx.length; i++ ) idx[i] = i;
        Arrays.sort(idx, new Comparator<Integer>() {
            public int compare(Integer i1, Integer i2) {
                return Integer.compare(inputArray[i1], inputArray[i2]);
            }
        });
        int[] index = new int[idx.length];
        for (int i = 0; i < index.length; i++) index[i] = idx[i];
        
        return index;
    }
    
    /**
     * Return an index of an array by descending order of value
     * @param array
     * @return 
     */
    public static int[] getIndicesByDescendingValue(double[] array) {
        int[] index = getIndicesByAscendingValue(array);
        ArrayUtils.reverse(index);
        return index;
    }
    
    /**
     * Return an index of an array by ascending order of value
     * @param array
     * @return 
     */
    public static int[] getIndicesByAscendingValue(double[] array) {
        double[] inputArray = new double[array.length];
        System.arraycopy(array, 0, inputArray, 0, array.length);
        Integer[] idx = new Integer[inputArray.length];
        for( int i = 0 ; i < idx.length; i++ ) idx[i] = i;
        Arrays.sort(idx, new Comparator<Integer>() {
            public int compare(Integer i1, Integer i2) {
                return Double.compare(inputArray[i1], inputArray[i2]);
            }
        });
        int[] index = new int[idx.length];
        for (int i = 0; i < index.length; i++) index[i] = idx[i];
        return index;
    }

    /**
     * Return an index of an array by descending order of value
     * @param array
     * @return
     */
    public static int[] getIndicesByDescendingValue(String[] array) {
        int[] index = getIndicesByAscendingValue(array);
        ArrayUtils.reverse(index);
        return index;
    }

    /**
     * Return an index of an array by ascending order of value
     * @param array
     * @return
     */
    public static int[] getIndicesByAscendingValue(String[] array) {
        String[] inputArray = new String[array.length];
        System.arraycopy(array, 0, inputArray, 0, array.length);
        Integer[] idx = new Integer[inputArray.length];
        for( int i = 0 ; i < idx.length; i++ ) idx[i] = i;
        Arrays.sort(idx, new Comparator<Integer>() {
            public int compare(Integer i1, Integer i2) {
                return inputArray[i1].compareTo(inputArray[i2]);
            }
        });
        int[] index = new int[idx.length];
        for (int i = 0; i < index.length; i++) index[i] = idx[i];
        return index;
    }

    /**
     * Return a list of index with a give size
     * @param size
     * @return
     */
    public static List<Integer> getIndexList (int size) {
        List<Integer> l = new ArrayList<>();
        for (int i = 0; i < size; i++) {
            l.add(i);
        }
        return l;
    }

    /**
     * Return a list of index from startIndex to the endIndex
     * @param startIndex inclusive
     * @param endIndex exclusive
     * @return
     */
    public static List<Integer> getIndexList (int startIndex, int endIndex) {
        List<Integer> l = new ArrayList<>(endIndex - startIndex);
        for (int i = startIndex; i < endIndex; i++) {
            l.add(i);
        }
        return l;
    }
}
