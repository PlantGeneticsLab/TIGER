/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.app.grt;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import pgl.infra.dna.snp.SNPOld;
import gnu.trove.list.array.TByteArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.list.array.TShortArrayList;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 *
 * @author feilu
 */
public class TagAnnotation implements Swapper, IntComparator {
    protected int groupIndex = -1;
    protected List<long[]> tagList = null;
    protected TByteArrayList r1LenList = null;
    protected TByteArrayList r2LenList = null;
    protected TIntArrayList readCountList = null;
    protected TShortArrayList r1ChrList = null;
    protected TShortArrayList r2ChrList = null;
    protected TIntArrayList r1StartPosList = null;
    protected TIntArrayList r2StartPosList = null;
    protected TByteArrayList r1StrandList = null;
    protected TByteArrayList r2StrandList = null;
    protected TByteArrayList r1MapQList = null;
    protected TByteArrayList r2MapQList = null;
    protected List<List<SNPOld>> SNPList = null;
    protected List<List<AlleleInfo>> alleleList = null;
    

    TagAnnotation (int groupIndex) {
        this.groupIndex = groupIndex;
        tagList = new ArrayList<>();
        r1LenList = new TByteArrayList();
        r2LenList = new TByteArrayList();
        readCountList = new TIntArrayList();
        r1ChrList = new TShortArrayList();
        r2ChrList = new TShortArrayList();
        r1StartPosList = new TIntArrayList();
        r2StartPosList = new TIntArrayList();
        r1StrandList = new TByteArrayList();
        r2StrandList = new TByteArrayList();
        r1MapQList = new TByteArrayList();
        r2MapQList = new TByteArrayList();
        SNPList = new ArrayList<>();
        alleleList = new ArrayList<>();
    }
    
    TagAnnotation (int tagLengthInLong, int groupIndex, int tagNumber, boolean ifSorted) {
        this.groupIndex = groupIndex;
        tagList = new ArrayList(tagNumber);
        r1LenList = new TByteArrayList(tagNumber);
        r2LenList = new TByteArrayList(tagNumber);
        readCountList = new TIntArrayList(tagNumber);
        r1ChrList = new TShortArrayList();
        r2ChrList = new TShortArrayList();
        r1StartPosList = new TIntArrayList();
        r2StartPosList = new TIntArrayList();
        r1StrandList = new TByteArrayList();
        r2StrandList = new TByteArrayList();
        r1MapQList = new TByteArrayList();
        r2MapQList = new TByteArrayList();
        SNPList = new ArrayList<>();
        alleleList = new ArrayList<>();
    }
    
    void appendTag (long[] tag, byte r1Len, byte r2Len, int readNumber) {
        tagList.add(tag);
        r1LenList.add(r1Len);
        r2LenList.add(r2Len);
        readCountList.add(readNumber);
        r1ChrList.add(Short.MIN_VALUE);
        r2ChrList.add(Short.MIN_VALUE);
        r1StartPosList.add(Integer.MIN_VALUE);
        r2StartPosList.add(Integer.MIN_VALUE);
        r1StrandList.add(Byte.MIN_VALUE);
        r2StrandList.add(Byte.MIN_VALUE);
        r1MapQList.add(Byte.MIN_VALUE);
        r2MapQList.add(Byte.MIN_VALUE);
        SNPList.add(new ArrayList<SNPOld>());
        alleleList.add(new ArrayList<AlleleInfo>());
    }
    
    void appendTag (long[] tag, byte r1Len, byte r2Len, int readNumber, short r1Chr, int r1Pos, byte r1Strand, byte r1MapQ, short r2Chr, int r2Pos, byte r2Strand, byte r2MapQ, List<SNPOld> tagSNPList, List<AlleleInfo> tagAlleleList) {
        tagList.add(tag);
        r1LenList.add(r1Len);
        r2LenList.add(r2Len);
        readCountList.add(readNumber);
        r1ChrList.add(r1Chr);
        r1StartPosList.add(r1Pos);
        r1StrandList.add(r1Strand);
        r1MapQList.add(r1MapQ);
        r2ChrList.add(r2Chr);
        r2StartPosList.add(r2Pos);
        r2StrandList.add(r2Strand);
        r2MapQList.add(r2MapQ);
        SNPList.add(tagSNPList);
        alleleList.add(tagAlleleList);
    }
    
    long[] getTag (int tagIndex) {
        return this.tagList.get(tagIndex);
    }
    
    byte getR1TagLength (int tagIndex) {
        return this.r1LenList.get(tagIndex);
    }
    
    byte getR2TagLength (int tagIndex) {
        return this.r2LenList.get(tagIndex);
    }
    
    short getR1Chromosome (int tagIndex) {
        return this.r1ChrList.get(tagIndex);
    }
    
    short getR2Chromosome (int tagIndex) {
        return this.r2ChrList.get(tagIndex);
    }
    
    int getR1StartPosition (int tagIndex) {
        return this.r1StartPosList.get(tagIndex);
    }
    
    int getR2StartPosition (int tagIndex) {
        return this.r2StartPosList.get(tagIndex);
    }
    
    byte getR1Strand (int tagIndex) {
        return this.r1StrandList.get(tagIndex);
    }
    
    byte getR2Strand (int tagIndex) {
        return this.r2StrandList.get(tagIndex);
    }
    
    byte getR1MapQ (int tagIndex) {
        return this.r1MapQList.get(tagIndex);
    }
    
    byte getR2MapQ (int tagIndex) {
        return this.r2MapQList.get(tagIndex);
    }
    
    List<SNPOld> getSNPOfTag (int tagIndex) {
        return this.SNPList.get(tagIndex);
    }
    
    
    List<AlleleInfo> getAlleleOfTag (int tagIndex) {
        return this.alleleList.get(tagIndex);
    }
    
    int getTotalReadNum () {
        int cnt = 0;
        for (int i = 0; i < this.getTagNumber(); i++) {
            cnt+=this.getReadNumber(i);
        }
        return cnt;
    }
    
    int getReadNumber (int tagIndex) {
        return this.readCountList.get(tagIndex);
    }
    
    int getTagNumber () {
        return this.tagList.size();
    }
    
    int getTagIndex (long[] tag) {
        return Collections.binarySearch(tagList, tag, TagUtils.tagCom);
    }  
    
    void setR1Chromosome (int tagIndex, short chr) {
        this.r1ChrList.set(tagIndex, chr);
    }
    
    void setR2Chromosome (int tagIndex, short chr) {
        this.r2ChrList.set(tagIndex, chr);
    }
    
    void setR1StartPosition (int tagIndex, int position) {
        this.r1StartPosList.set(tagIndex, position);
    }
    
    void setR2StartPosition (int tagIndex, int position) {
        this.r2StartPosList.set(tagIndex, position);
    }
    
    void setR1Strand (int tagIndex, byte strand) {
        this.r1StrandList.set(tagIndex, strand);
    }
    
    void setR2Strand (int tagIndex, byte strand) {
        this.r2StrandList.set(tagIndex, strand);
    }
    
    void setR1MapQ (int tagIndex, byte mapQ) {
        this.r1MapQList.set(tagIndex, mapQ);
    }
    
    void setR2MapQ (int tagIndex, byte mapQ) {
        this.r2MapQList.set(tagIndex, mapQ);
    }
    
    void setSNPOfTag (int tagIndex, List<SNPOld> tagSNPList) {
        this.SNPList.set(tagIndex, tagSNPList);
    }
    
    void setAlleleOfTag (int tagIndex, List<AlleleInfo> tagAlleleList) {
        this.alleleList.set(tagIndex, tagAlleleList);
    }
    
    void removeSNPOfTag (int tagIndex) {
        this.SNPList.get(tagIndex).clear();
    }
    
    void removeAlleleOfTag (int tagIndex) {
        this.alleleList.get(tagIndex).clear();
    }
    
    @Override
    public void swap (int index1, int index2) {
        long[] temp = tagList.get(index1);
        tagList.set(index1, tagList.get(index2));
        tagList.set(index2, temp);
        byte tb = r1LenList.get(index1);
        r1LenList.set(index1, r1LenList.get(index2));
        r1LenList.set(index2, tb);
        tb = r2LenList.get(index1);
        r2LenList.set(index1, r2LenList.get(index2));
        r2LenList.set(index2, tb);
        int tc = readCountList.get(index1);
        readCountList.set(index1, readCountList.get(index2));
        readCountList.set(index2, tc);
        short ts = r1ChrList.get(index1);
        r1ChrList.set(index1, r1ChrList.get(index2));
        r1ChrList.set(index2, ts);
        ts = r2ChrList.get(index1);
        r2ChrList.set(index1, r2ChrList.get(index2));
        r2ChrList.set(index2, ts);
        tc = r1StartPosList.get(index1);
        r1StartPosList.set(index1, r1StartPosList.get(index2));
        r1StartPosList.set(index2, tc);
        tc = r2StartPosList.get(index1);
        r2StartPosList.set(index1, r2StartPosList.get(index2));
        r2StartPosList.set(index2, tc);
        tb = r1StrandList.get(index1);
        r1StrandList.set(index1, r1StrandList.get(index2));
        r1StrandList.set(index2, tb);
        tb = r2StrandList.get(index1);
        r2StrandList.set(index1, r2StrandList.get(index2));
        r2StrandList.set(index2, tb);
        tb = r1MapQList.get(index1);
        r1MapQList.set(index1, r1MapQList.get(index2));
        r1MapQList.set(index2, tb);
        tb = r2MapQList.get(index1);
        r2MapQList.set(index1, r2MapQList.get(index2));
        r2MapQList.set(index2, tb);
        List<SNPOld> tempSNP = SNPList.get(index1);
        SNPList.set(index1, SNPList.get(index2));
        SNPList.set(index2, tempSNP);
        List<AlleleInfo> tempAllele = alleleList.get(index1);
        alleleList.set(index1, alleleList.get(index2));
        alleleList.set(index2, tempAllele);
    }

    @Override
    public int compare (int index1, int index2) {
        for (int i = 0; i < this.tagList.get(0).length; i++) {
            if (tagList.get(index1)[i] < tagList.get(index2)[i]) {
                return -1;
            }
            if (tagList.get(index1)[i] > tagList.get(index2)[i]) {
                return 1;
            }
        }
        return 0;
    }
    
    void sortAlleleListByPosition (int tagIndex) {
        Collections.sort(this.alleleList.get(tagIndex));
    }
  
    void sortSNPListByPosition (int tagIndex) {
        Collections.sort(this.SNPList.get(tagIndex));
    }
    
    void sort () {
        //System.out.println("TagCount sort begins");
        GenericSorting.quickSort(0, this.getTagNumber(), this, this);
        //System.out.println("TagCount sort ends");
    }
    
    protected int collapseCounts (int minReadCount) {
        int collapsedRows = 0;
        for (int i = 0; i < this.getTagNumber()-1; i++) {
            if (this.readCountList.get(i) == 0) continue;
            for (int j = i + 1; j < this.getTagNumber(); j++) {
                int index = this.compare(i, j);
                if (index < 0) break;
                else {
                    int sum = readCountList.get(i)+readCountList.get(j);
                    readCountList.set(i, sum);
                    collapsedRows++;
                    readCountList.set(j, 0);
                }
            }
        }
        List<long[]> aTagList = new ArrayList<>();
        TByteArrayList aR1LenList = new TByteArrayList();
        TByteArrayList aR2LenList = new TByteArrayList();
        TIntArrayList aReadCountList = new TIntArrayList();
        TShortArrayList aR1ChrList = new TShortArrayList();
        TShortArrayList aR2ChrList = new TShortArrayList();
        TIntArrayList aR1StartPosList = new TIntArrayList();
        TIntArrayList aR2StartPosList = new TIntArrayList();
        TByteArrayList aR1StrandList = new TByteArrayList();
        TByteArrayList aR2StrandList = new TByteArrayList();
        TByteArrayList aR1MapQList = new TByteArrayList();
        TByteArrayList aR2MapQList = new TByteArrayList();
        List<List<SNPOld>> aSNPList = new ArrayList<>();
        List<List<AlleleInfo>> aAlleleList = new ArrayList<>();
        for (int i = 0; i < this.getTagNumber(); i++) {
            if (readCountList.get(i) < minReadCount) continue;
            aTagList.add(tagList.get(i));
            aR1LenList.add(r1LenList.get(i));
            aR2LenList.add(r2LenList.get(i));
            aReadCountList.add(readCountList.get(i));
            aR1ChrList.add(r1ChrList.get(i));
            aR2ChrList.add(r2ChrList.get(i));
            aR1StartPosList.add(r1StartPosList.get(i));
            aR2StartPosList.add(r2StartPosList.get(i));
            aR1StrandList.add(r1StrandList.get(i));
            aR2StrandList.add(r2StrandList.get(i));
            aR1MapQList.add(r1MapQList.get(i));
            aR2MapQList.add(r2MapQList.get(i));
            aSNPList.add(SNPList.get(i));
            aAlleleList.add(alleleList.get(i));
        }
        tagList = aTagList;
        r1LenList = aR1LenList;
        r2LenList = aR2LenList;
        readCountList = aReadCountList;
        r1ChrList = aR1ChrList;
        r2ChrList = aR2ChrList;
        r1StartPosList = aR1StartPosList;
        r2StartPosList = aR2StartPosList;
        r1StrandList = aR1StrandList;
        r2StrandList = aR2StrandList;
        r1MapQList = aR1MapQList;
        r2MapQList = aR2MapQList;
        SNPList = aSNPList;
        alleleList = aAlleleList;
        aTagList = null;
        aR1LenList = null;
        aR2LenList = null;
        aReadCountList = null;
        aR1ChrList = null;
        aR2ChrList = null;
        aR1StartPosList = null;
        aR2StartPosList = null;
        aR1StrandList = null;
        aR2StrandList = null;
        aR1MapQList = null;
        aR2MapQList = null;
        aSNPList = null;
        aAlleleList = null;
        return collapsedRows;
    }
}
