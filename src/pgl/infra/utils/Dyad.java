/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.utils;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;

/**
 *
 * @author feilu
 * @param <T>
 * @param <S>
 */
public class Dyad<T, S> {
    T first = null;
    S second = null;
    
    public Dyad (T first, S second) {
        this.first = first;
        this.second = second;
    }
    
    public T getFirstElement () {
        return first;
    }
    
    public S getSecondElement () {
        return second;
    }
    
    public boolean sortByFirstIntInt () {
        if (first.getClass().getSimpleName().equals("int[]")) {
            if (second.getClass().getSimpleName().equals("int[]")) {
                if (((int[])first).length != ((int[])second).length) return false;
                GenericSorting.quickSort(0, ((int[])first).length, compByFirstIntInt, swapperIntInt);
                return true;
            }
            return false;
        }
        return false;
    }
    
    Swapper swapperIntInt = new Swapper() {
        @Override
        public void swap(int index1, int index2) {
            int temp = ((int[])first)[index1];
            ((int[])first)[index1] = ((int[])first)[index2];
            ((int[])first)[index2] = temp;
            temp = ((int[])second)[index1];
            ((int[])second)[index1] = ((int[])second)[index2];
            ((int[])second)[index2] = temp;
        }
    };
    
    IntComparator compByFirstIntInt = new IntComparator() {
        @Override
        public int compare(int index1, int index2) {
            return ((int[])first)[index1]-((int[])first)[index2];
        }
    };
}
