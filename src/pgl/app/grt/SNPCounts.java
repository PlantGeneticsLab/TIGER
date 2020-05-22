/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.app.grt;

import pgl.infra.dna.snp.SNPOld;
import pgl.infra.pos.ChrPos;
import gnu.trove.list.array.TByteArrayList;
import gnu.trove.set.hash.TShortHashSet;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import pgl.infra.utils.IOUtils;

/**
 *
 * @author feilu
 */
public class SNPCounts {
    List<List<SNPCount>> chrSCLists = null;
    short[] chrs = null;
    
    public SNPCounts (TagAnnotations tas) {
        this.initilize(tas);
    }
    
    public SNPCounts (String infileS) {
        this.readBinaryFile(infileS);
    }
    
    public void selectSNPs (List<ChrPos> posList) {
        List<List<SNPCount>> nChrSCLists = new ArrayList<>();
        for (int i = 0; i < chrSCLists.size(); i++) {
            List<SNPCount> oList = chrSCLists.get(i);
            List<SNPCount> nList = new ArrayList<>();
            for (int j = 0; j < oList.size(); j++) {
                int index = Collections.binarySearch(posList, oList.get(j).getChrPos());
                if (index < 0) continue;
                nList.add(oList.get(j));
            }
            if (nList.size() == 0) continue;
            nChrSCLists.add(nList);
        }
        short[] nchrs = new short[nChrSCLists.size()];
        for (int i = 0; i <  nchrs.length; i++) {
            nchrs[i] = nChrSCLists.get(i).get(0).getChromosome();
        }
        chrSCLists = nChrSCLists;
        chrs = nchrs;
        nChrSCLists = null;
        nchrs = null;
    }
    
    public int getChromosomeNumber () {
        return chrs.length;
    }
    
    public short getChromosome (int chrIndex) {
        return chrs[chrIndex];
    }
    
    public int getSNPNumberOnChromosome (int chrIndex) {
        return chrSCLists.get(chrIndex).size();
    }
    
    public int getTotalSNPNumber () {
        int cnt = 0;
        for (int i = 0; i < this.getChromosomeNumber(); i++) {
            cnt+=this.getSNPNumberOnChromosome(i);
        }
        return cnt;
    }
    
    public int getReadNumberOfSNP (int chrIndex, int snpIndex) {
        return chrSCLists.get(chrIndex).get(snpIndex).getReadNumber();
    }
    
    public int getTagNumberOfSNP (int chrIndex, int snpIndex) {
        return chrSCLists.get(chrIndex).get(snpIndex).getTagNumber();
    }
    
    public int getPositionOfSNP (int chrIndex, int snpIndex) {
        return chrSCLists.get(chrIndex).get(snpIndex).getPosition();
    }
    
    public short getChromosomeOfSNP (int chrIndex, int snpIndex) {
        return chrSCLists.get(chrIndex).get(snpIndex).getChromosome();
    }
    
    public byte getRefAlleleByteOfSNP (int chrIndex, int snpIndex) {
        return chrSCLists.get(chrIndex).get(snpIndex).getRefAlleleByte();
    }
    
    public byte getAltAlleleByteOfSNP (int chrIndex, int snpIndex, int altIndex) {
        return chrSCLists.get(chrIndex).get(snpIndex).getAltAlleleByte(altIndex);
    }
    
    public int getAltAlleleIndex (int chrIndex, int snpIndex, byte alt) {
        return chrSCLists.get(chrIndex).get(snpIndex).getAltAlleleIndex(alt);
    }
    
    public int getAltAlleleIndex (int chrIndex, int snpIndex, char altAllele) {
        return chrSCLists.get(chrIndex).get(snpIndex).getAltAlleleIndex(altAllele);
    }
    
    public long getTotalReadNumber () {
        long cnt = 0;
        for (int i = 0; i < this.getChromosomeNumber(); i++) {
            for (int j = 0; j < this.getSNPNumberOnChromosome(i); j++) {
                cnt+=this.getReadNumberOfSNP(i, j);
            }
        }
        return cnt;
    }
    
    public void readBinaryFile (String infileS) {
        System.out.println("Reading SNPCounts file");
        try {
            DataInputStream dis  = IOUtils.getBinaryReader(infileS);
            short chrNumber = dis.readShort();
            chrs = new short[chrNumber];
            for (int i = 0; i < chrNumber; i++) {
                chrs[i] = dis.readShort();
            }
            chrSCLists = new ArrayList<>();
            long cnt = 0;
            for (int i = 0; i < chrNumber; i++) {
                List<SNPCount> cl = new ArrayList<>();
                int size = dis.readInt();
                for (int j = 0; j < size; j++) {
                    short chr = dis.readShort();
                    int pos = dis.readInt();
                    byte ref = dis.readByte();
                    byte altNumber = dis.readByte();
                    TByteArrayList alts = new TByteArrayList();
                    for (int k = 0; k < altNumber; k++) {
                        alts.add(dis.readByte());
                    }
                    SNPCount cs = new SNPCount(chr, pos, ref, alts, dis.readInt(), dis.readInt());
                    cl.add(cs);
                    cnt++;
                    if (cnt%1000000 == 0) System.out.println(String.valueOf(cnt) + " SNPs are read in");
                }
                chrSCLists.add(cl);
            }
            dis.close();
            System.out.println("A total of " + String.valueOf(cnt) + " SNPs are read in from " + infileS);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void writeBinaryFile (String outfileS) {
        System.out.println("Writing SNPCounts file");
        try {
            DataOutputStream dos = IOUtils.getBinaryWriter(outfileS);
            dos.writeShort(chrs.length);
            for (int i = 0; i < chrs.length; i++) {
                dos.writeShort(chrs[i]);
            }
            long cnt = 0;
            for (int i = 0; i < chrSCLists.size(); i++) {
                List<SNPCount> cl = chrSCLists.get(i);
                dos.writeInt(cl.size()); 
                for (int j = 0; j < cl.size(); j++) {
                    SNPCount s = cl.get(j);
                    dos.writeShort(s.getChromosome());
                    dos.writeInt(s.getPosition());
                    dos.writeByte(s.getRefAlleleByte());
                    dos.writeByte(s.getAltAlleleNumber());
                    for (int k = 0; k < s.getAltAlleleNumber(); k++) {
                        dos.writeByte(s.getAltAlleleByte(k));
                    }
                    dos.writeInt(s.getReadNumber());
                    dos.writeInt(s.getTagNumber());
                    cnt++;
                    if (cnt%1000000 == 0) System.out.println(String.valueOf(cnt) + " SNPs are written out");
                }
            }
            dos.flush();
            dos.close();
            System.out.println("A total of " + String.valueOf(cnt) + " SNPs are written to " + outfileS);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void writeBinaryFile (String outfileS, int minReadCount) {
        System.out.println("Writing SNPCounts file");
        try {
            DataOutputStream dos = IOUtils.getBinaryWriter(outfileS);
            dos.writeShort(chrs.length);
            for (int i = 0; i < chrs.length; i++) {
                dos.writeShort(chrs[i]);
            }
            long cnt = 0;
            for (int i = 0; i < chrSCLists.size(); i++) {
                List<SNPCount> cl = chrSCLists.get(i);
                int count = 0;
                for (int j = 0; j < cl.size(); j++) {
                    SNPCount s = cl.get(j);
                    if (s.getReadNumber() < minReadCount) continue;
                    count++;
                }
                dos.writeInt(count); 
                for (int j = 0; j < cl.size(); j++) {
                    SNPCount s = cl.get(j);
                    if (s.getReadNumber() < minReadCount) continue;
                    dos.writeShort(s.getChromosome());
                    dos.writeInt(s.getPosition());
                    dos.writeByte(s.getRefAlleleByte());
                    dos.writeByte(s.getAltAlleleNumber());
                    for (int k = 0; k < s.getAltAlleleNumber(); k++) {
                        dos.writeByte(s.getAltAlleleByte(k));
                    }
                    dos.writeInt(s.getReadNumber());
                    dos.writeInt(s.getTagNumber());
                    cnt++;
                    if (cnt%1000000 == 0) System.out.println(String.valueOf(cnt) + " SNPs are written out");
                }
            }
            dos.flush();
            dos.close();
            System.out.println("A total of " + String.valueOf(cnt) + " SNPs are written to " + outfileS);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public int getChrIndex (short chr) {
        return Arrays.binarySearch(chrs, chr);
    }
    
    public int getSNPIndex (int chrIndex, ChrPos query) {
        int snpIndex = Collections.binarySearch(this.chrSCLists.get(chrIndex), query);
        return snpIndex;
    }
    
    public int getSNPIndex (int chrIndex, int pos) {
        ChrPos query = new ChrPos (chrs[chrIndex], pos);
        int snpIndex = Collections.binarySearch(this.chrSCLists.get(chrIndex), query);
        return snpIndex;
    }
    
    public void sort () {
        chrSCLists.parallelStream().forEach(l -> {
            Collections.sort(l);
        });
    }
    
    public void mergeSNPs () {
        System.out.println("Merging identical or multi-allele SNPs in SNPCounts");
        AtomicInteger acnt = new AtomicInteger();
        chrSCLists.parallelStream().forEach(l -> {
            int collapsedRows = 0;     
            for (int i = 0; i < l.size()-1; i++) {
                if (l.get(i).getReadNumber() == 0) continue;
                for (int j = i + 1; j < l.size(); j++) {
                    int index = l.get(i).compareTo(l.get(j));//
                    if (index < 0) break;
                    else {
                        int readSum = l.get(i).getReadNumber()+l.get(j).getReadNumber();
                        int tagSum = l.get(i).getTagNumber()+l.get(j).getTagNumber();
                        l.get(i).setReadNumber(readSum);
                        l.get(i).setTagNumber(tagSum);
                        for (int k = 0; k < l.get(j).getAltAlleleNumber(); k++) {
                            l.get(i).addAltAlleleByte(l.get(j).getAltAlleleByte(k));
                        }
                        collapsedRows++;
                        l.get(j).setReadNumber(0);
                        l.get(j).setTagNumber(0);
                    }
                }
            }
            for (int i = 0; i < l.size(); i++) {
                if (l.get(i).getReadNumber() != 0) {
                    l.get(i).removeDuplicatedAltAlleles();
                    continue;
                }
                l.remove(i);
                i--;
            }
            acnt.addAndGet(collapsedRows);
        });
        System.out.println(String.valueOf(acnt) + " SNPs are merged");
        System.out.println(String.valueOf(this.getTotalSNPNumber()) + " unique SNPs remained in SNPCounts");
    }
    
    private void initilize (TagAnnotations tas) {
        System.out.println("Initializing SNPCounts from TagAnnotations");
        TShortHashSet chrSet = new TShortHashSet();
        for (int i = 0; i < tas.groupCount; i++) {
            for (int j = 0; j < tas.getTagNumber(i); j++) {
                List<SNPOld> tagSNPList = tas.getSNPOfTag(i, j);
                if (tagSNPList.size() == 0) continue;
                chrSet.add(tagSNPList.get(0).getChromosome());
            }
        }
        chrs = chrSet.toArray();
        Arrays.sort(chrs);
        System.out.println("Tags in TagAnnotations file are mapped to " + String.valueOf(chrs.length) + " chromosomes");
        chrSCLists = new ArrayList<>();
        for (int i = 0; i < chrs.length; i++) {
            List<SNPCount> al = new ArrayList<>();
            chrSCLists.add(al);
        }
        int index = -1;
        int cnt = 0;
        for (int i = 0; i < tas.getGroupNumber(); i++) {
            for (int j = 0; j < tas.getTagNumber(i); j++) {
                List<SNPOld> tagSNPList = tas.getSNPOfTag(i, j);
                if (tagSNPList.size() == 0) continue;
                index = Arrays.binarySearch(chrs, tagSNPList.get(0).getChromosome());
                for (int k = 0; k < tagSNPList.size(); k++) {
                    SNPCount sc = new SNPCount (tagSNPList.get(k), tas.getReadNumber(i, j), 1);
                    chrSCLists.get(index).add(sc);
                    cnt++;
                    if (cnt%100000 == 0) System.out.println(String.valueOf(cnt) + " SNPs found");
                }
            }
        }
        this.sort();
        System.out.println("A total of " + String.valueOf(cnt) + " SNPs are found from TagAnnotations");
        this.mergeSNPs();
    }
}

class SNPCount extends SNPOld {
    int readCount = 0;
    int tagCount = 0;
    
    public SNPCount(short chr, int pos, byte ref, byte alt, int readCount, int tagCount) {
        super(chr, pos, ref, alt);
        this.readCount = readCount;
        this.tagCount = tagCount;
    }
    
    public SNPCount(short chr, int pos, byte ref, TByteArrayList alts, int readCount, int tagCount) {
        super(chr, pos, ref, alts);
        this.readCount = readCount;
        this.tagCount = tagCount;
    }
    
    public SNPCount(SNPOld snp, int readCount, int tagCount) {
        super(snp.getChromosome(), snp.getPosition(), snp.getRefAlleleByte(), snp.getAltAlleleList());
        this.readCount = readCount;
        this.tagCount = tagCount;
    }
    
    public int getReadNumber () {
        return readCount;
    }
    
    public void setReadNumber (int readCount) {
        this.readCount = readCount;
    }
    
    public int getTagNumber () {
        return tagCount;
    }
    
    public void setTagNumber (int tagCount) {
        this.tagCount = tagCount;
    }
}
