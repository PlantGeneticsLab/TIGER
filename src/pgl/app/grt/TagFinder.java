/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.app.grt;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import pgl.infra.dna.BaseEncoder;
import gnu.trove.list.array.TIntArrayList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import pgl.infra.utils.PArrayUtils;
import pgl.infra.utils.Dyad;

/**
 *
 * @author feilu
 */
public class TagFinder {
    int paraLevel = 32;
    TagAnnotations tas = null;
    GroupTagFinder[] gtfs = null;
    int ceaseSearchingSize = 3;
    int maxQueryExistance = 1000;
    
    public TagFinder (TagAnnotations tas) {
        this.tas = tas;
        this.initialize(tas);
    }
    
    public Dyad<int[], int[]> getMostSimilarTags (long[] tag, byte r1Length, byte r2Length, int groupIndex, int maxDivergence) {
        TIntArrayList r1IntSeqList = this.getR1IntSeq(tag, r1Length);
        TIntArrayList r2IntSeqList = this.getR2IntSeq(tag, r2Length);       
        int[] startEndIndex = null;
        int minSize = Integer.MAX_VALUE;
        boolean ifCease = false;
        for (int i = 0; i < r1IntSeqList.size(); i++) {
            int[] currentSE = gtfs[groupIndex].getStartEndIndex(r1IntSeqList.get(i));
            if (currentSE == null) continue;
            else {
                int currentSize = currentSE[1] - currentSE[0];
                if (currentSize < minSize) {
                    minSize = currentSize;
                    startEndIndex = currentSE;
                }
                if (minSize < ceaseSearchingSize) {
                    ifCease = true;
                    break;
                }
            }
        }
        if (!ifCease) {
            for (int i = 0; i < r2IntSeqList.size(); i++) {
                int[] currentSE = gtfs[groupIndex].getStartEndIndex(r2IntSeqList.get(i));
                if (currentSE == null) continue;
                else {
                    int currentSize = currentSE[1] - currentSE[0];
                    if (currentSize < minSize) {
                        minSize = currentSize;
                        startEndIndex = currentSE;
                    }
                }
                if (minSize < ceaseSearchingSize) {
                    ifCease = true;
                    break;
                }
            }
        }
        if (startEndIndex == null) return null;
        TIntArrayList mismatchList = new TIntArrayList();
        TIntArrayList tagIndicesList = new TIntArrayList();
        for (int i = startEndIndex[0]; i < startEndIndex[1]; i++) {
            int dbTagIndex = gtfs[groupIndex].getTagIndex(i);
            long[] dbTag = tas.getTag(groupIndex, dbTagIndex);
            
            int currentMismatch = 0;
            boolean toContinue = false;
            for (int j = 0; j < dbTag.length; j++) {
                byte diff = BaseEncoder.getSeqDifferences(dbTag[j], tag[j], maxDivergence);
                if (diff > maxDivergence) {
                    toContinue = true;
                    break;
                }
                currentMismatch+= diff;
            }
            if (toContinue) continue;
            if (currentMismatch > maxDivergence) continue;
            mismatchList.add(currentMismatch);
            tagIndicesList.add(dbTagIndex);
        }
        if (mismatchList.isEmpty()) return null;
        int[] mismatch = mismatchList.toArray();
        int[] tagIndices = tagIndicesList.toArray();
        Dyad<int[], int[]> t = new Dyad<>(mismatch, tagIndices);
        t.sortByFirstIntInt();
        mismatch = t.getFirstElement();
        tagIndices = t.getSecondElement();
        int minMismatch = mismatch[0];
        if (minMismatch > maxDivergence) return null;
        int size = 1;
        for (int i = 1; i < mismatch.length; i++) {
            if (mismatch[i] == minMismatch) size++;
            else break;
        }
        int[] mismatchArray = new int[size];
        int[] tagIndicesArray = new int[size];
        for (int i = 0; i < size; i++) {
            mismatchArray[i] = mismatch[i];
            tagIndicesArray[i] = tagIndices[i];
        }
        Dyad<int[], int[]> result = new Dyad<>(mismatchArray, tagIndicesArray);
        return result;
    }
    
    private void initialize (TagAnnotations tas) {
        System.out.println("Start building TagFinder from TagAnnotations");
        int groupNumber  = tas.getGroupNumber();
        gtfs = new GroupTagFinder[groupNumber];
        List<TagAnnotation> taList = tas.taList;
        int[][] indices = PArrayUtils.getSubsetsIndicesBySubsetSize(taList.size(), this.paraLevel);
        for (int i = 0; i < indices.length; i++) {
            List<Integer> subIndexList = new ArrayList();
            int size = indices[i][1]-indices[i][0];
            for (int j = 0; j < size; j++) {
                subIndexList.add(j+indices[i][0]);
            }
            subIndexList.parallelStream().forEach(currentIndex -> {
                TagAnnotation ta = taList.get(currentIndex);
                TIntArrayList intQueryList = new TIntArrayList();
                TIntArrayList indexList = new TIntArrayList();
                int tagNumber  = ta.getTagNumber();
                for (int j = 0; j < tagNumber; j++) {
                    long[] tag = ta.getTag(j);
                    byte r1Length = ta.getR1TagLength(j);
                    byte r2Length = ta.getR2TagLength(j);
                    TIntArrayList tagInt = this.getIntSeq(tag, r1Length, r2Length);                    
                    intQueryList.addAll(tagInt);
                    for (int k = 0; k < tagInt.size(); k++) indexList.add(j);
                }
                GroupTagFinder gsf = new GroupTagFinder(intQueryList, indexList, this.maxQueryExistance);
                gtfs[currentIndex] = gsf;
            });
        }
        System.out.println("Finish building TagFinder");
    }
    
    private TIntArrayList getR2IntSeq (long[] tag, byte r2Length) {
        TIntArrayList seqList = new TIntArrayList();
        int[] intTag = BaseEncoder.getIntsFromLongs(tag);
        int n = r2Length/BaseEncoder.intChunkSize;
        for (int i = 0; i < n; i++) {
            seqList.add(intTag[i+intTag.length/2]);
        }
        return seqList;
    }
    
    private TIntArrayList getR1IntSeq (long[] tag, byte r1Length) {
        TIntArrayList seqList = new TIntArrayList();
        int[] intTag = BaseEncoder.getIntsFromLongs(tag);
        int n = r1Length/BaseEncoder.intChunkSize;
        for (int i = 0; i < n; i++) {
            seqList.add(intTag[i]);
        }
        return seqList;
    }
    
    private TIntArrayList getIntSeq (long[] tag, byte r1Length, byte r2Length) {
        TIntArrayList seqList = new TIntArrayList();
        int[] intTag = BaseEncoder.getIntsFromLongs(tag);
        int n = r1Length/BaseEncoder.intChunkSize;
        for (int i = 0; i < n; i++) {
            seqList.add(intTag[i]);
        }
        n = r2Length/BaseEncoder.intChunkSize;
        for (int i = 0; i < n; i++) {
            seqList.add(intTag[i+intTag.length/2]);
        }
        return seqList;
    }
}

class GroupTagFinder {
    int[] intQuery;
    int[] tagIndices;
    
    public GroupTagFinder (TIntArrayList intQueryList, TIntArrayList indexList, int maxQueryExistence) {
        this.intQuery = intQueryList.toArray();
        this.tagIndices = indexList.toArray();
        this.sort();
        this.removeDuplicatedQuery(maxQueryExistence);
    }
    
    public void removeDuplicatedQuery (int maxQueryExistence) {
        if (intQuery.length == 0) return;
        TIntArrayList intQueryList = new TIntArrayList(intQuery);
        TIntArrayList indexList = new TIntArrayList(tagIndices);
        int cnt = 0;
        int currentIndex = 0;
        int query = intQueryList.get(0);
        for (int i = 0; i < intQueryList.size(); i++) {
            if (intQueryList.get(i) == query) {
                cnt++;
            }
            else {
                query = intQueryList.get(i);
                if (cnt > maxQueryExistence) {
                    intQueryList.remove(currentIndex, cnt);
                    indexList.remove(currentIndex, cnt);
                    i = currentIndex-1;
                    cnt = 0;
                }
                else {
                    currentIndex = i;
                    cnt = 1;
                }
            }
        }
        if (cnt > maxQueryExistence) {
            intQueryList.remove(currentIndex, cnt);
            indexList.remove(currentIndex, cnt);
        }
        intQuery = intQueryList.toArray();
        tagIndices = indexList.toArray();
    }
    
    public int getTagIndex (int intSeqIndex) {
        return tagIndices[intSeqIndex];
    }
    
    public int[] getStartEndIndex (int seq) {
        int index = this.binarySearch(seq);
        if (index < 0) return null;
        int[] startEnd = new int[2];
        startEnd[0] = index;//inclusive
        startEnd[1] = index+1;//exclusive
        for (int i = index-1; i > -1; i--) {
            if (intQuery[i] == seq) {
                startEnd[0] = i;
            }
            else break;
        }
        for (int i = index+1; i < intQuery.length; i++) {
            startEnd[1] = i;
            if (intQuery[i] != seq) {
                break;
            }
            else {
                if (i == intQuery.length-1) startEnd[1] = intQuery.length;
            }
        }
        return startEnd;
    }
    
    public int binarySearch (int seq) {
        return Arrays.binarySearch(intQuery, seq);
    }
    
    public void sort () {
        GenericSorting.quickSort(0, intQuery.length, comp, swapper);
    }
    
    Swapper swapper = new Swapper() {
        @Override
        public void swap(int index1, int index2) {
            int temp = intQuery[index1];
            intQuery[index1] = intQuery[index2];
            intQuery[index2] = temp;
            temp = tagIndices[index1];
            tagIndices[index1] = tagIndices[index2];
            tagIndices[index2] = temp;
        }
    };
    
    IntComparator comp = new IntComparator() {
        @Override
        public int compare(int index1, int index2) {
            if (intQuery[index1] < intQuery[index2]) return -1;
            else if (intQuery[index1] > intQuery[index2]) return 1;
            else {
                if (tagIndices[index1] < tagIndices[index2]) return -1;
                else if (tagIndices[index1] > tagIndices[index2]) return 1;
                return 0;
            }
        }
    };
}
