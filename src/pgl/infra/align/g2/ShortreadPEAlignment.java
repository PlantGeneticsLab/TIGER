/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package pgl.infra.align.g2;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;


/**
 * @deprecated
 * Hold paired end alignment information from aligners (Bowtie2 and BWA-MEM)
 * @author Fei Lu
 */
public class ShortreadPEAlignment {
    PEAlignmentInfo[] pai = null;
    byte sortType = 0;
    
    public ShortreadPEAlignment () {}
    
    public ShortreadPEAlignment (String alignmentFileS, IOFileFormat f) {
        if (f == IOFileFormat.Text) this.readTextSimplePEAlignment(alignmentFileS);
        else if (f == IOFileFormat.Binary) this.readBinarySimplePEAlignment(alignmentFileS);
        else if (f == IOFileFormat.TextGzip) this.readTextGzipSimplePEAlignment(alignmentFileS);
        else throw new UnsupportedOperationException("Not supported yet.");
    }
    
    private void readTextGzipSimplePEAlignment (String alignmentFileS) {
        try {
            BufferedReader br = IOUtils.getTextGzipReader(alignmentFileS);
            String temp;
            br.readLine();
            int cnt = 0;
            ArrayList<PEAlignmentInfo> aList = new ArrayList();
            while ((temp = br.readLine()) != null) {
                PEAlignmentInfo info = new PEAlignmentInfo(temp);
                aList.add(info);
                if (cnt%100000 == 0) System.out.println(String.valueOf(cnt+1) + " is read from anlignmentFileS");
                cnt++;
            }
            pai = aList.toArray(new PEAlignmentInfo[aList.size()]);
            br.close();
            this.sortByQuery();
            System.out.println(alignmentFileS+" is load and sorted by query");
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    private void readBinarySimplePEAlignment (String alignmentFileS) {
        try {
            DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(alignmentFileS), 65536));
            this.sortType = dis.readByte();
            this.pai = new PEAlignmentInfo[dis.readInt()];
            for (int i = 0; i < pai.length; i++) {
                pai[i] = new PEAlignmentInfo(dis.readUTF(), dis.readUTF(), dis.readInt(), dis.readInt(), dis.readByte(), dis.readShort(), dis.readShort(), dis.readShort(), dis.readBoolean(), dis.readUTF(), dis.readInt(), dis.readInt(), dis.readByte(), dis.readShort(), dis.readShort(), dis.readShort(), dis.readBoolean());
                if (i%5000000 == 0) System.out.println("Read in " + String.valueOf(i) + " alignments");
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private void readTextSimplePEAlignment (String alignmentFileS) {
        try {
            BufferedReader br = new BufferedReader (new FileReader(alignmentFileS), 65536);
            String temp;
            br.readLine();
            int cnt = 0;
            ArrayList<PEAlignmentInfo> aList = new ArrayList();
            while ((temp = br.readLine()) != null) {
                PEAlignmentInfo info = new PEAlignmentInfo(temp);
                aList.add(info);
                if (cnt%100000 == 0) System.out.println(String.valueOf(cnt+1) + " is read from anlignmentFileS");
                cnt++;
            }
            pai = aList.toArray(new PEAlignmentInfo[aList.size()]);
            br.close();
            this.sortByQuery();
            System.out.println(alignmentFileS+" is load and sorted by query");
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    /**
     * Output SimplePEAlignment file
     * @param outfileS
     * @param f 
     */
    public void writeShortreadPEAlignment (String outfileS, IOFileFormat f) {
        if (f == IOFileFormat.Text) this.writeTextSimplePEAlignment(outfileS);
        else if (f == IOFileFormat.Binary) this.writeBinarySimplePEAlignment(outfileS);
        else throw new UnsupportedOperationException("Not supported yet.");
    }
    
    private void writeTextSimplePEAlignment (String outfileS) {
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(outfileS), 65536);
            bw.write("Query\tChromosomeF\tStartPosF\tEndPosF\tStrandF\tMappingQualityF\tIfPerfectMatchF\tChromosomeB\tStartPosB\tEndPosB\tStrandB\tMappingQualityB\tIfPerfectMatchB");
            bw.newLine();
            for (int i = 0; i < this.getAlignmentNumber(); i++) {
                bw.write(pai[i].toString());
                bw.newLine();
                if (i%1000000 == 0) System.out.println("wrote " + String.valueOf(i) + "PE alignments");
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    private void writeBinarySimplePEAlignment (String outfileS) {
        try {
            DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfileS), 65536));
            dos.writeByte(this.sortType);
            dos.writeInt(this.getAlignmentNumber());
            for (int i = 0; i < this.getAlignmentNumber(); i++) {
                dos.writeUTF(this.getQuery(i));
                dos.writeUTF(this.getHitF(i));
                dos.writeInt(this.getStartPosF(i));
                dos.writeInt(this.getEndPosF(i));
                dos.writeByte(this.getStrandF(i));
                dos.writeShort(this.getMappingQualityF(i));
                dos.writeShort(this.getMatchNumberF(i));
                dos.writeShort(this.getEditDistanceF(i));
                dos.writeBoolean(this.isPerfectMatchF(i));
                dos.writeUTF(this.getHitB(i));
                dos.writeInt(this.getStartPosB(i));
                dos.writeInt(this.getEndPosB(i));
                dos.writeByte(this.getStrandB(i));
                dos.writeShort(this.getMappingQualityB(i));
                dos.writeShort(this.getMatchNumberB(i));
                dos.writeShort(this.getEditDistanceB(i));
                dos.writeBoolean(this.isPerfectMatchB(i));
                if (i%1000000 == 0) System.out.println("wrote " + String.valueOf(i) + " alignments");
            }
            dos.flush();
            dos.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    /**
     * Make a SimplePEAlignment instance from BWAMEM sam file
     * Ignore those singletons and secondary alignment
     * @param samFileS 
     */
    public void readFromBWAMEM (String samFileS) {
        ArrayList<PEAlignmentInfo> aList = new ArrayList();
        System.out.println("Reading SAM format PE alignment (BWAMEM) from: " + samFileS);
        try {
            BufferedReader br;
            if (samFileS.endsWith(".gz")) {
                br = IOUtils.getTextGzipReader(samFileS);
            } else {
                br = IOUtils.getTextReader(samFileS);
            }
            while(br.readLine().startsWith("@PG")==false) {};
            String temp;
            String[] forwardArray = null;
            String[] backwardArray = null;
            ShortreadAlignment.AlignmentInfo fInfo = null;
            ShortreadAlignment.AlignmentInfo bInfo = null;
            ShortreadAlignment sa = new ShortreadAlignment();
            int cnt = 0;
            ArrayList<String[]> infoL = new ArrayList();
            infoL.add(br.readLine().split("\\s"));
            while((temp = br.readLine())!=null) {
                String[] tem = temp.split("\\s");
                if (tem[0].equals(infoL.get(infoL.size()-1)[0])) {
                    infoL.add(tem);
                }
                else {
                    Iterator<String[]> it = infoL.iterator();
                    while (it.hasNext()) {
                        String[] te = it.next();
                        int flag = Integer.valueOf(te[1]);
                        //the 11th(index) bit is essentially a bug of bwa-mem, the 11th bit is not defined in SAM protocol
                        if (sa.getBitFromInt(flag, 11) == 1 || sa.getBitFromInt(flag, 8) == 1 || sa.getBitFromInt(flag, 10) == 1) it.remove();
                    }
                    if (infoL.size() != 2) {
                        System.out.println(infoL.get(0)[0]);
                    }
                    if (sa.getBitFromInt(Integer.valueOf(infoL.get(0)[1]), 6) == 1) {
                        forwardArray = infoL.get(0);
                        backwardArray = infoL.get(1);
                    }
                    else {
                        forwardArray = infoL.get(1);
                        backwardArray = infoL.get(0);
                    }
                    fInfo = sa.getAlignmentInfoFromBWAMEM(forwardArray);
                    bInfo = sa.getAlignmentInfoFromBWAMEM(backwardArray);
                    aList.add(new PEAlignmentInfo(fInfo, bInfo));
                    infoL.clear();
                    infoL.add(tem);
                }
                cnt++;
                if (cnt%500000 == 0) System.out.println("Read in " + String.valueOf(cnt) + " lines");
            }
            if (!infoL.isEmpty()) {
                Iterator<String[]> it = infoL.iterator();
                while (it.hasNext()) {
                    String[] te = it.next();
                    int flag = Integer.valueOf(te[1]);
                    //the 11th(index) bit is essentially a bug of bwa-mem, the 11th bit is not defined in SAM protocol
                    if (sa.getBitFromInt(flag, 11) == 1 || sa.getBitFromInt(flag, 8) == 1 || sa.getBitFromInt(flag, 10) == 1) it.remove();
                }
                if (infoL.size() != 2) {
                    System.out.println(infoL.get(0)[0]);
                }
                forwardArray = infoL.get(0);
                backwardArray = infoL.get(1);
                fInfo = sa.getAlignmentInfoFromBWAMEM(forwardArray);
                bInfo = sa.getAlignmentInfoFromBWAMEM(backwardArray);
                aList.add(new PEAlignmentInfo(fInfo, bInfo));
            }
            pai = aList.toArray(new PEAlignmentInfo[aList.size()]);
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        this.sortByQuery();
        System.out.println("SimplePEAlignment has "+String.valueOf(pai.length) + " alignments");
    }
    
    public void readFromBlast (String r1AlignmentFileS, String r2AlignmentFileS, double eThresh, int readLength) {
        ArrayList<PEAlignmentInfo> aList = new ArrayList();
        System.out.println("Reading blast result (outfmt 7) from: " + r1AlignmentFileS + " and " + r2AlignmentFileS);
        ShortreadAlignment saR1 = new ShortreadAlignment();
        saR1.readFromBlast(r1AlignmentFileS, eThresh, readLength);
        ShortreadAlignment saR2 = new ShortreadAlignment();
        saR2.readFromBlast(r2AlignmentFileS, eThresh, readLength);
        String[] querys = saR1.getQuerys();
        for (int i = 0; i < querys.length; i++) {
            ShortreadAlignment.AlignmentInfo r1 = saR1.ai[saR1.getAlignmentStartIndexByQuery(querys[i])];
            ShortreadAlignment.AlignmentInfo r2 = saR2.ai[saR2.getAlignmentStartIndexByQuery(querys[i])];
            aList.add(new PEAlignmentInfo(r1, r2));
            if (i%500000 == 0) System.out.println("Processed " + String.valueOf(i) + " PE alignment");
        }
        pai = aList.toArray(new PEAlignmentInfo[aList.size()]);
        this.sortByQuery();
        System.out.println("SimplePEAlignment has "+String.valueOf(pai.length) + " alignments");
    }
    
    /**
     * Make a SimplePEAlignment instance from Bowtie2 sam file
     * @param samFileS 
     */
    public void readFromBowtie2 (String samFileS) {
        ArrayList<PEAlignmentInfo> aList = new ArrayList();
        System.out.println("Reading SAM format PE alignment (Bowtie2) from: " + samFileS);
        try {
            BufferedReader br;
            if (samFileS.endsWith(".gz")) {
                br = IOUtils.getTextGzipReader(samFileS);
            } else {
                br = IOUtils.getTextReader(samFileS);
            }
            while(br.readLine().startsWith("@PG")==false) {};
            String forwardStr = null;
            String backwardStr = null;
            ShortreadAlignment.AlignmentInfo fInfo = null;
            ShortreadAlignment.AlignmentInfo bInfo = null;
            ShortreadAlignment sa = new ShortreadAlignment();
            int cnt = 0;
            while((forwardStr = br.readLine())!=null) {
                backwardStr = br.readLine();
                fInfo = sa.getAlignmentInfoFromBowtie2(forwardStr);
                bInfo = sa.getAlignmentInfoFromBowtie2(backwardStr);
                aList.add(new PEAlignmentInfo(fInfo, bInfo));
                cnt++;
                if (cnt%500000 == 0) System.out.println("Read in " + String.valueOf(cnt) + " lines");
            }
            pai = aList.toArray(new PEAlignmentInfo[aList.size()]);
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        this.sortByQuery();
        System.out.println("SimplePEAlignment has "+String.valueOf(pai.length) + " alignments");
    }
    
    public double getPEMappedRatio () {
        int cnt = 0;
        for (int i = 0; i < this.getAlignmentNumber(); i++) {
            if (this.isMatchF(i) && this.isMatchB(i)) cnt++;
        }
        return (double)cnt/this.getAlignmentNumber();
    }
    
    public double getSEMapptedRatio () {
        int cnt = 0;
        for (int i = 0; i < this.getAlignmentNumber(); i++) {
            if (this.isMatchF(i))  {
                if (!this.isMatchB(i)) cnt++;
            }
            else {
                if (this.isMatchB(i)) cnt++;
            }
        }
        return (double)cnt/this.getAlignmentNumber();
    }
    
    public double getNoMapptedRatio () {
        int cnt = 0;
        for (int i = 0; i < this.getAlignmentNumber(); i++) {
            if (!this.isMatchF(i) && !this.isMatchB(i)) cnt++;
        }
        return (double)cnt/this.getAlignmentNumber();
    }
    
    /**
     * Return name of query
     * @param index
     * @return 
     */
    public String getQuery (int index) {
        return pai[index].query;
    }
    
    /**
     * Return name of hit of forward read, return empty String if the query does not match 
     * @param index
     * @return 
     */
    public String getHitF (int index) {
        return pai[index].hitF;
    }
    
    /**
     * Return start position of the alignment of forward read, return Integer.MIN_VALUE if the query does not match
     * @param index
     * @return 
     */
    public int getStartPosF (int index) {
        return pai[index].startPosF;
    }
    
    /**
     * Return end position of the alignment of forward read, return Integer.MIN_VALUE if the query does not match
     * @param index
     * @return 
     */
    public int getEndPosF (int index) {
        return pai[index].endPosF;
    }
    
    /**
     * Return strand of alignment of forward read, return Byte.MIN_VALUE if the query does not match
     * @param index
     * @return 
     */
    public byte getStrandF (int index) {
        return pai[index].strandF;
    }
    
    /**
     * Return mapping quality score of alignment of forward read, return Short.MIN_VALUE if the query does not match
     * @return 
     */
    public short getMappingQualityF (int index) {
        return pai[index].mappingQualityF;
    }
    
    /**
     * Return the length of matched range of forward read, return Short.MIN_VALUE if the query does not match
     * @param index
     * @return 
     */
    public short getMatchNumberF (int index) {
        return pai[index].matchNumberF;
    }
    
    /**
     * Return mismatch number in matched range of forward read, return Short.MIN_VALUE if the query does not match
     * @param index
     * @return 
     */
    public short getEditDistanceF (int index) {
        return pai[index].editDistanceF;
    }
    
    public float getEditDistanceRatioF (int index) {
        if (this.getEditDistanceF(index) < 0 || this.getMatchNumberF(index) < 0) return Float.NaN;
        return (float)this.getEditDistanceF(index)/(float)this.getMatchNumberF(index);
    }
    
    /**
     * Return if the alignment of forward read is a perfect match
     * @param index
     * @return 
     */
    public boolean isPerfectMatchF (int index) {
        return pai[index].ifPerfectMatchF;
    }
    
    /**
     * Return name of hit of backward read, return empty String if the query does not match 
     * @param index
     * @return 
     */
    public String getHitB (int index) {
        return pai[index].hitB;
    }
    
    /**
     * Return start position of the alignment of backward read, return Integer.MIN_VALUE if the query does not match
     * @param index
     * @return 
     */
    public int getStartPosB (int index) {
        return pai[index].startPosB;
    }
    
    /**
     * Return end position of the alignment of backward read, return Integer.MIN_VALUE if the query does not match
     * @param index
     * @return 
     */
    public int getEndPosB (int index) {
        return pai[index].endPosB;
    }
    
    /**
     * Return strand of alignment of backward read, return Byte.MIN_VALUE if the query does not match
     * @param index
     * @return 
     */
    public byte getStrandB (int index) {
        return pai[index].strandB;
    }
    
    /**
     * Return mapping quality score of alignment of backward read, return Short.MIN_VALUE if the query does not match
     * @return 
     */
    public short getMappingQualityB (int index) {
        return pai[index].mappingQualityB;
    }
    
    /**
     * Return the length of matched range of backward read, return Short.MIN_VALUE if the query does not match
     * @param index
     * @return 
     */
    public short getMatchNumberB (int index) {
        return pai[index].matchNumberB;
    }
    
    /**
     * Return mismatch number in matched range of backward read, return Short.MIN_VALUE if the query does not match
     * @param index
     * @return 
     */
    public short getEditDistanceB (int index) {
        return pai[index].editDistanceB;
    }
    
    public float getEditDistanceRatioB (int index) {
        if (this.getEditDistanceB(index) < 0 || this.getMatchNumberB(index) < 0) return Float.NaN;
        return (float)this.getEditDistanceB(index)/(float)this.getMatchNumberB(index);
    }
    
    public float getEditDistanceRatioBothEnd (int index) {
        int ed = 0;
        int matches = 0;
        if (this.getEditDistanceF(index) >= 0) ed+=this.getEditDistanceF(index);
        if (this.getEditDistanceB(index) >= 0) ed+=this.getEditDistanceB(index);
        if (this.getMatchNumberF(index) >= 0) matches+=this.getMatchNumberF(index);
        if (this.getMatchNumberB(index) >= 0) matches+=this.getMatchNumberB(index);
        if (matches < 0 || ed < 0) return Float.NaN;
        return (float)ed/matches;
    }
    /**
     * Return if the alignment of backward read is a perfect match
     * @param index
     * @return 
     */
    public boolean isPerfectMatchB (int index) {
        return pai[index].ifPerfectMatchB;
    }
    
    /**
     * Return if the forward read matches
     * @param index
     * @return 
     */
    public boolean isMatchF (int index) {
        if (this.getStartPosF(index)!=Integer.MIN_VALUE) return true;
        return false;
    }
    
    /**
     * Return if the backward read matches
     * @param index
     * @return 
     */
    public boolean isMatchB (int index) {
        if (this.getStartPosB(index)!=Integer.MIN_VALUE) return true;
        return false;
    }
    
    /**
     * Return if the forward and backward reads both match
     * @param index
     * @return 
     */
    public boolean isMatchBothEnds (int index) {
        if (this.isMatchF(index) && this.isMatchB(index)) return true;
        return false;
    }
    
    /**
     * Return if the forward and backward reads are aligned in reversed direction
     * @param index
     * @return 
     */
    public boolean isOnReverseDirectionFB (int index) {
        if (this.getStrandF(index)*this.getStrandB(index) == -1) return true;
        return false;
    }
    
    /**
     * Return best estimated fragment size of PE based on alignment. Return -1 if PE doesn't meet requirements, including:
     * Same chromosome, reversed direction of PE
     * @param index
     * @return 
     */
    public int getPEFragmentSize (int index) {
        if (this.isOnReverseDirectionFB(index)) {
            //if (!this.getHitF(index).equals(this.getHitB(index))) return -1;
            if (this.getStrandF(index) == 1) {
                if (this.getStartPosF(index) > this.getEndPosB(index)) return -1;
                return this.getEndPosB(index) - this.getStartPosF(index) + 1;
            }
            else {
                if (this.getEndPosF(index) < this.getStartPosB(index)) return -1;
                return this.getEndPosF(index) - this.getStartPosB(index) + 1;
            }
        }
        return -1;
    }
    
    /**
     * Return the number of alignment
     * @return 
     */
    public int getAlignmentNumber () {
        return pai.length;
    }
    
    public void sortByQuery () {
        sortType = 0;
        Arrays.sort(pai);
    }
    
    class PEAlignmentInfo implements Comparable <PEAlignmentInfo> {
        /**Query sequence*/
        String query = null;
        /**Chromosome forward*/
        String hitF = "";
        /**Chromosome backward*/
        String hitB = "";
        /**Reference starting position of forward sequence*/
        int startPosF = Integer.MIN_VALUE;
        /**Reference starting position of backward sequence*/
        int startPosB = Integer.MIN_VALUE;
        /**Reference ending position of forward sequence*/
        int endPosF = Integer.MIN_VALUE;
        /**Reference ending position of backward sequence*/
        int endPosB = Integer.MIN_VALUE;
        /**Strand of alignment of forward sequence, + = 1, - = -1*/
        byte strandF = Byte.MIN_VALUE;
        /**Strand of alignment of backward sequence, + = 1, - = -1*/
        byte strandB = Byte.MIN_VALUE;
        /**The length of matched range of forward sequence*/
        short matchNumberF = Short.MIN_VALUE;
        /**The length of matched range of backward sequence*/
        short matchNumberB = Short.MIN_VALUE;
        /**Mismatch number in the matched range of forward sequence*/
        short editDistanceF = Short.MIN_VALUE;
        /**Mismatch number in the matched range of backward sequence*/
        short editDistanceB = Short.MIN_VALUE;
        /**Mapping quality of forward sequence*/
        short mappingQualityF = Short.MIN_VALUE;
        /**Mapping quality of backward sequence*/
        short mappingQualityB = Short.MIN_VALUE;
        /**If forward sequence is a perfect match to reference*/
        boolean ifPerfectMatchF = false;
        /**If backward sequence is a perfect match to reference*/
        boolean ifPerfectMatchB = false;
        
        public PEAlignmentInfo (ShortreadAlignment.AlignmentInfo fInfo, ShortreadAlignment.AlignmentInfo bInfo) {
            this.query = fInfo.query;
            this.hitF = fInfo.hit;
            this.startPosF = fInfo.startPos;
            this.endPosF = fInfo.endPos;
            this.strandF = fInfo.strand;
            this.matchNumberF = fInfo.matchNumber;
            this.editDistanceF = fInfo.editDistance;
            this.mappingQualityF = fInfo.mappingQuality;
            this.ifPerfectMatchF = fInfo.ifPerfectMatch;
            this.hitB = bInfo.hit;
            this.startPosB = bInfo.startPos;
            this.endPosB = bInfo.endPos;
            this.strandB = bInfo.strand;
            this.matchNumberB = bInfo.matchNumber;
            this.editDistanceB = bInfo.editDistance;
            this.mappingQualityB = bInfo.mappingQuality;
            this.ifPerfectMatchB = bInfo.ifPerfectMatch;
        }
        
        public PEAlignmentInfo (String inputStr) {
            String[] temp = inputStr.split("\t");
            this.query = temp[0];
            this.hitF = temp[1];
            this.startPosF = Integer.valueOf(temp[2]);
            this.endPosF = Integer.valueOf(temp[3]);
            this.strandF = Byte.valueOf(temp[4]);
            this.mappingQualityF = Short.valueOf(temp[5]);
            this.matchNumberF = Short.valueOf(temp[6]);
            this.editDistanceF = Short.valueOf(temp[7]);
            this.ifPerfectMatchF = Boolean.valueOf(temp[8]);
            this.hitB = temp[9];
            this.startPosB = Integer.valueOf(temp[10]);;
            this.endPosB = Integer.valueOf(temp[11]);;
            this.strandB = Byte.valueOf(temp[12]);
            this.mappingQualityB = Short.valueOf(temp[13]);
            this.matchNumberB = Short.valueOf(temp[14]);
            this.editDistanceB = Short.valueOf(temp[15]);
            this.ifPerfectMatchB = Boolean.valueOf(temp[16]);
        }
        
        public PEAlignmentInfo (String query, String hitF, int startPosF, int endPosF, byte strandF, short mappingQualityF, short matchNumberF, short editDistanceF, boolean ifPerfectMatchF, String hitB, int startPosB, int endPosB, byte strandB, short mappingQualityB, short matchNumberB, short editDistanceB, boolean ifPerfectMatchB) {
            this.query = query;
            this.hitF = hitF;
            this.startPosF = startPosF;
            this.endPosF = endPosF;
            this.strandF = strandF;
            this.mappingQualityF = mappingQualityF;
            this.matchNumberF = matchNumberF;
            this.editDistanceF = editDistanceF;
            this.ifPerfectMatchF = ifPerfectMatchF;
            this.hitB = hitB;
            this.startPosB = startPosB;
            this.endPosB = endPosB;
            this.strandB = strandB;
            this.mappingQualityB = mappingQualityB;
            this.matchNumberB = matchNumberB;
            this.editDistanceB = editDistanceB;
            this.ifPerfectMatchB = ifPerfectMatchB;
        }
        
        @Override
        public String toString () {
            StringBuilder sb = new StringBuilder();
            sb.append(query).append("\t").append(hitF).append("\t").append(startPosF).append("\t").append(endPosF).append("\t").append(strandF).append("\t").append(mappingQualityF).append("\t");
            sb.append(matchNumberF).append("\t").append(editDistanceF).append("\t");
            if (ifPerfectMatchF) sb.append("1");
            else sb.append("0");
            sb.append(hitB).append("\t").append(startPosB).append("\t").append(endPosB).append("\t").append(strandB).append("\t").append(mappingQualityB).append("\t");
            sb.append(matchNumberB).append("\t").append(editDistanceB).append("\t");
            if (ifPerfectMatchB) sb.append("1");
            else sb.append("0");
            return sb.toString();
        }
        
        @Override
        public int compareTo(PEAlignmentInfo o) {
            if (sortType == 0) {
                return query.compareTo(o.query);
            }
            return 0;
        }
    }
}
