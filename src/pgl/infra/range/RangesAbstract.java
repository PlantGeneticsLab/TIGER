/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.range;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import com.koloboke.collect.set.hash.HashIntSet;
import com.koloboke.collect.set.hash.HashIntSets;
import gnu.trove.list.array.TIntArrayList;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Basic method implementation of {@link RangesInterface}, providing functions of sorting, searching, and merging, etc.
 * @author feilu
 */
public abstract class RangesAbstract implements RangesInterface {

    //statistics of current range list, need to rebuild after write operations on ranges
    //0, unsorted; 1, by position; 2, by size; 3 by value
    protected int sortType = 0;
    protected int[] chrs = null;
    
    
    protected void resetStatistics () {
        this.chrs = null;
        this.sortType = 0;
    }
    
    @Override
    public void sortBySize() {
        GenericSorting.quickSort(0, this.getRangeNumber(), compBySize, swapper);
        this.sortType = 2;
    }

    @Override
    public void sortByStartPosition() {
        GenericSorting.quickSort(0, this.getRangeNumber(), compByStartPosition, swapper);
        this.sortType = 1;
    }


    protected Swapper swapper = new Swapper() {
        @Override
        public void swap(int a, int b) {
            Range temp = getRange(a);
            setRange(a, getRange(b));
            setRange(b, temp);
        }
    };
    
    protected IntComparator compBySize = new IntComparator() {
        @Override
        public int compare(int a, int b) {
            int sA = getRange(a).getRangeSize();
            int sB = getRange(b).getRangeSize();
            return sA-sB;
        }
    };
    
    protected IntComparator compByStartPosition = new IntComparator() {
        @Override
        public int compare(int a, int b) {
            int chrA = getRangeChromosome(a);
            int chrB = getRangeChromosome(b);
            if (chrA == chrB) {
                int sA = getRange(a).getRangeStart();
                int sB = getRange(b).getRangeStart();
                return sA-sB;
            }
            else {
                return chrA-chrB;
            }
        }
    };

    @Override
    public int getStartIndexOfChromosome(int chr) {
        if (sortType != 1) this.sortByStartPosition();
        RangeInterface query = new Range(chr, Integer.MIN_VALUE, Integer.MIN_VALUE);
        int hit  = Collections.binarySearch(getRangeList(), query);
        if (hit < 0) {
            int index = -hit-1;
            if (this.getRangeChromosome(index) == chr) return index;
            return hit;
        }
        return hit;
    }

    @Override
    public int getEndIndexOfChromosome(int chr) {
        if (sortType != 1) this.sortByStartPosition();
        RangeInterface query = new Range(chr+1, Integer.MIN_VALUE, Integer.MIN_VALUE);
        int hit  = Collections.binarySearch(getRangeList(), query);
        if (hit < 0) {
            int index = -hit-2;
            if (this.getRangeChromosome(index) == chr) return index;
            return hit;
        }
        return hit;
    }

    @Override
    public int[] getChromosomes() {
        if (this.chrs != null) return chrs;
        HashIntSet s = HashIntSets.getDefaultFactory().newMutableSet();
        for (int i = 0; i < this.getRangeNumber(); i++) {
            s.add(this.getRangeChromosome(i));
        }
        this.chrs = s.toIntArray();
        Arrays.sort(chrs);
        return chrs;
    }

    @Override
    public int getChromosomeNumber() {
        return this.getChromosomes().length;
    }

    @Override
    public int getRangeStart(int rangeIndex) {
        return this.getRange(rangeIndex).getRangeStart();
    }

    @Override
    public int getRangeEnd(int rangeIndex) {
        return this.getRange(rangeIndex).getRangeEnd();
    }

    @Override
    public int getRangeChromosome(int rangeIndex) {
        return this.getRange(rangeIndex).getRangeChromosome();
    }
    
    @Override
    public int[] getRangesIndicesContainsPosition (int chr, int pos) {
        if (sortType != 1) this.sortByStartPosition();
        TIntArrayList indexList = new TIntArrayList();
        Range query = new Range(chr, pos, pos+1);
        int hit = Collections.binarySearch(getRangeList(), query);
        if (hit < 0) hit = -hit-1;
        while (this.getRange(hit).isContain(chr, pos)) {
            indexList.add(hit);
            hit++;
            if (hit > this.getRangeNumber()-1) break;
        }
        return indexList.toArray();
    }
    
    @Override
    public Ranges getNonOverlapRanges() {
        int[] chromosomes = this.getChromosomes();
        int[] mins = new int[chromosomes.length];
        int[] maxs = new int[chromosomes.length];
        for (int i = 0; i < mins.length; i++) {
            mins[i] = Integer.MAX_VALUE;
            maxs[i] = Integer.MIN_VALUE;
        }
        for (int i = 0; i < this.getRangeNumber(); i++) {
            int index = Arrays.binarySearch(chromosomes, this.getRangeChromosome(i));
            int v = this.getRangeStart(i);
            if (v < mins[index]) mins[index] = v;
            v = this.getRangeEnd(i);
            if (v > maxs[index]) maxs[index] = v;
        }
        List<Range> rList = new ArrayList<>();
        for (int i = 0; i < chromosomes.length; i++) {
            System.out.println("Start collasping chromosome " + String.valueOf(chromosomes[i]));
            int base = mins[i];
            int length = maxs[i] - mins[i];
            byte[] status = new byte[length];
            int startIndex = this.getStartIndexOfChromosome(chromosomes[i]);
            int endIndex = this.getEndIndexOfChromosome(chromosomes[i]);
            for (int j = startIndex; j < endIndex; j++) {
                for (int k = this.getRangeStart(j); k < this.getRangeEnd(j); k++) {
                    status[k-base] = 1;
                }
            }
            int current = 0;
            while (current < length) {
                if (status[current] != 0) {
                    int start = current+base;
                    while (current < length && status[current] == 1) {
                        current++;
                    }
                    int end = current+base;
                    rList.add(new Range(chromosomes[i], start, end));
                }
                current++;
            }
        }
        return new Ranges(rList);
    }

    @Override
    public int getFirstRangeIndex(int chr, int pos) {
        if (sortType != 1) this.sortByStartPosition();
        Range query = new Range(chr, pos, pos+1);
        int hit = Collections.binarySearch(getRangeList(), query);
        int index;
        if (hit < 0) {
            index = -hit-1;
            if (this.getRange(index).isContain(chr, pos)) return index;
            return hit;
        }
        return hit;
    }

}
