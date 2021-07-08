/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.app.grt;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import gnu.trove.list.array.TByteArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.TByteHashSet;
import java.util.Arrays;

/**
 *
 * @author feilu
 */
public class AlleleDepth {
    TByteArrayList alleleList = null;
    TIntArrayList depthList = null;
    byte[] alleles = null;
    int[] depth = null;
    byte sortType = -1;
    
    public AlleleDepth () {
        alleleList = new TByteArrayList();
        depthList = new TIntArrayList();
    }
    
    public AlleleDepth (byte[] alleles, int[] depth) {
        this.alleles = alleles;
        this.depth = depth;
    }
    
    public void addAllele (byte allele) {
        alleleList.add(allele);
    }
    
    public void addDepth (int depth) {
        depthList.add(depth);
    }
    
    public void toArray () {
        TByteHashSet s = new TByteHashSet(this.alleleList);
        alleles = s.toArray();
        Arrays.sort(alleles);
        depth = new int[alleles.length];
        for (int i = 0; i < alleleList.size(); i++) {
            int index = Arrays.binarySearch(alleles, alleleList.get(i));
            depth[index]+=depthList.get(i);
        }
        sortType = 0;
        this.alleleList = null;
        this.depthList = null;
    }
    
    public byte[] getAlleles () {
        return this.alleles; 
    }
    
    public int[] getDepths () {
        return this.depth;
    }
    
    public byte getAllele (int index) {
        return alleles[index];
    }
    
    public int getDepth (int index) {
        return depth[index];
    }
    
    public int getTotalDepth () {
        int cnt = 0;
        for (int i = 0; i < this.getAlleleNumber(); i++) {
            cnt+=this.getDepth(i);
        }
        return cnt;
    }
    
    public int getAlleleNumber () {
        return alleles.length;
    }
    
    public int getDepth (byte allele) {
        if (sortType != 0) this.sortByAllele();
        int alleleIndex = this.binarySearch(allele);
        if (alleleIndex < 0) return 0;
        return this.getDepth(alleleIndex);
    }
    
    public int[] getDepth (byte allele1, byte allele2) {
        int[] dep = new int[2];
        dep[0] = this.getDepth(allele1);
        dep[1] = this.getDepth(allele2);
        return dep; 
    }
    
    public AlleleDepth getAltAlleleDepth (byte ref) {
        if (sortType != 0) this.sortByAllele();
        int refIndex = this.binarySearch(ref);
        TByteArrayList alleleList = new TByteArrayList(this.alleles);
        TIntArrayList depthList = new TIntArrayList(this.depth);
        if (!(refIndex < 0)) {
            alleleList.removeAt(refIndex);
            depthList.removeAt(refIndex);
        }
        return new AlleleDepth(alleleList.toArray(), depthList.toArray());
    }
    
    public int binarySearch (byte allele) {
        if (sortType != 0) this.sortByAllele();
        return Arrays.binarySearch(this.alleles, allele);
    }
    
    public void sortByAllele () {
        GenericSorting.quickSort(0, alleles.length, compByAllele, swapper);
        sortType = (byte)0;
    }
    
    public void sortByDepthAsending () {
        GenericSorting.quickSort(0, alleles.length, compByDepthAsending, swapper);
        sortType = (byte)1;
    }
    
    public void sortByDepthDesending () {
        GenericSorting.quickSort(0, alleles.length, compByDepthDesending, swapper);
        sortType = (byte)2;
    }
    
    Swapper swapper = new Swapper() {
        @Override
        public void swap(int index1, int index2) {
            byte tempByte = alleles[index1];
            alleles[index1] = alleles[index2];
            alleles[index2] = tempByte;
            int tempInt = depth[index1];
            depth[index1] = depth[index2];
            depth[index2] = tempInt;
        }
    };
    
    IntComparator compByAllele = new IntComparator() {
        @Override
        public int compare(int index1, int index2) {
            return alleles[index1] - alleles[index2];
        }
    };
    
    IntComparator compByDepthAsending = new IntComparator() {
        @Override
        public int compare(int index1, int index2) {
            return depth[index1] - depth[index2];
        }
    };
    
    IntComparator compByDepthDesending = new IntComparator() {
        @Override
        public int compare(int index1, int index2) {
            return depth[index2] - depth[index1];
        }
    };
}
