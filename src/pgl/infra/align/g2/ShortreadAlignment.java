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
import java.util.Comparator;
import java.util.TreeSet;

import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;

/**
 * Holding info from short read alignment result, including tools to convert results from aligners (Bowtie2 and BWA-MEM)
 * Alignments with no hits are included
 * Can be sort by query name, also (hit name and position)
 * @deprecated 
 * @author Fei Lu
 */
public class ShortreadAlignment {
    AlignmentInfo[] ai;
    /**0 by query (Default). 1 by hit and pos.*/
    byte sortType = 0; 
    
    /**
     * Construct a {@link ShortreadAlignment} object
     */
    public ShortreadAlignment() {}
    
    /**
     * Constructs the object from alignment files, supporting ".txt" or ".txt.gz"
     * @param alignmentFileS
     * @param f 
     */
    public ShortreadAlignment(String alignmentFileS, IOFileFormat f) {
        if (f == IOFileFormat.Text) this.readTextSimpleAlignment(alignmentFileS);
        else if (f == IOFileFormat.Binary) this.readBinarySimpleAlignment(alignmentFileS);
        else throw new UnsupportedOperationException("Not supported yet.");
    }
    
    /**
     * Return number of alignment
     * @return 
     */
    public int getAlignmentNumber () {
        return ai.length;
    }
    
    /**
     * Note: Need to be sorted by hit and pos first, otherwise program quit
     * @param hit
     * @return 
     */
    public int getAlignmentNumberByHit (String hit) {
        if (this.sortType != 1) {
            System.out.println("Alignment should be sorted by hit and pos first. Program quits");
            System.exit(0);
        }
        return getAlignmentEndIndexByHit(hit) - getAlignmentStartIndexByHit(hit);
    }
    
    /**
     * Note: Need to be sorted by query first, otherwise program quits
     * @param query
     * @return 
     */
    public int getAlignmentNumberByQuery (String query) {
        if (this.sortType != 0) {
            System.out.println("Alignment should be sorted by query first. Program quits");
            System.exit(0);
        }
        return getAlignmentEndIndexByQuery(query) - getAlignmentStartIndexByQuery(query);
    }
    
    /**
     * Return -1 if the hit doesn't exist
     * Note: Need to be sorted by hit and pos first, otherwise program quits
     * @param hit
     * @return 
     */
    public int getAlignmentStartIndexByHit (String hit) {
        if (this.sortType != 1) {
            System.out.println("Alignment should be sorted by hit and pos first. Program quits");
            System.exit(0);
        }
        int index = this.getAlignmentIndexByHitPos(hit, 0);
        if (index < 0) {
            index = -index -1;
            if (ai[index].hit.equals(hit)) return index;
            else return -1;
        }
        else {
            while (index > 0 && ai[index-1].hit.equals(hit)) {
                index--;
            }
            return index;
        }
    }
    
    /**
     * Return alignment index
     * Note: Need to be sorted by hit and pos first, otherwise program quits
     * @param hit
     * @param pos
     * @return 
     */
    public int getAlignmentIndexByHitPos (String hit, int pos) {
        if (this.sortType != 1) {
            System.out.println("Alignment should be sorted by hit and pos first. Program quits");
            System.exit(0);
        }
        return Arrays.binarySearch(ai, new AlignmentInfo(hit, pos));
    }
    
    /**
     * Note: Need to be sorted by query first, otherwise program quites
     * @param query
     * @return 
     */
    public int getAlignmentStartIndexByQuery (String query) {
        if (this.sortType != 0) {
            System.out.println("Alignment should be sorted by query first. Program quits");
            System.exit(0);
        }
        int index = Arrays.binarySearch(ai, new AlignmentInfo(query, null));
        while (index > 0 && ai[index-1].query.equals(query)) {
            index--;
        }
        return index;
    }
    
    /**
     * 
     * Note: Need to be sorted by query first
     * @param queryStartWith Partial query from the start
     * @return starting index of alignment whose query starts with queryStartWith, return Integer.MIN_VALUE when there is queryStartWith doesn't exist in any query
     */
    public int getAlignmentStartIndexByQueryStartWith (String queryStartWith) {
        if (this.sortType != 0) {
            System.out.println("Alignment should be sorted by query first. Program quits");
            System.exit(0);
        }
        int index = Arrays.binarySearch(ai, new AlignmentInfo(queryStartWith, null));
        if (index < 0) {
            index = - index - 1;
            if (this.getQuery(index).startsWith(queryStartWith)) return index;
            return Integer.MIN_VALUE;
        }
        else {
            while (index > 0 && ai[index-1].query.equals(queryStartWith)) {
                index--;
            }
        }
        return index;
    }
    
    /**
     * Return the index of the last alignment of a hit, exclusive.
     * Return -1 if the hit doesn't exist
     * Note: Need to be sorted by hit first
     * @param hit
     * @return 
     */
    public int getAlignmentEndIndexByHit (String hit) {
        if (this.sortType != 1) {
            System.out.println("Alignment should be sorted by hit and pos first. Program quits");
            System.exit(0);
        }
        int index = this.getAlignmentIndexByHitPos(hit, Integer.MAX_VALUE);
        if (index < 0) {
            index = - index - 2;
            if (ai[index].hit.equals(hit)) return index+1;
            else return -1;
        }
        else {
            while ((index+1) < this.getAlignmentNumber() && ai[index+1].hit.equals(hit)) {
                index++;
            }
            return index+1;
        }
    }
    
    /**
     * Return the index of the last alignment of a query, exclusive
     * Note: Need to be sorted by query first
     * @param query
     * @return 
     */
    public int getAlignmentEndIndexByQuery (String query) {
        if (this.sortType != 0) {
            System.out.println("Alignment should be sorted by query first. Program quits");
            System.exit(0);
        }
        int index = Arrays.binarySearch(ai, new AlignmentInfo(query, null));
        if (index < 0) return index;
        while ((index+1) < this.getAlignmentNumber() && ai[index+1].query.equals(query)) {
            index++;
        }
        return index+1;
    }
    
    /**
     * Return the index of the last alignment of a queryStartWith, exclusive
     * @param queryStartWith
     * @return Ending index of alignment whose query starts with queryStartWith, return Integer.MIN_VALUE when there is queryStartWith doesn't exist in any query
     */
    public int getAlignmentEndIndexByQueryStartWith (String queryStartWith) {
        if (this.sortType != 0) {
            System.out.println("Alignment should be sorted by query first. Program quits");
            System.exit(0);
        }
        int index = this.getAlignmentStartIndexByQueryStartWith(queryStartWith);
        if (index < 0) return index;
        while ((index+1) < this.getAlignmentNumber() && ai[index+1].query.startsWith(queryStartWith)) {
            index++;
        }
        return index+1;
    }
    
    /**
     * Return non-redundant names of queries
     * @return 
     */
    public String[] getQuerys () {
        TreeSet<String> querySet = new TreeSet();
        for (int i = 0; i < this.getAlignmentNumber(); i++) {
            querySet.add(ai[i].query);
        }
        return querySet.toArray(new String[querySet.size()]);
    }
    
    /**
     * Return non-redundant names of hits
     * @return 
     */
    public String[] getHits () {
        TreeSet<String> hitSet = new TreeSet();
        for (int i = 0; i < this.getAlignmentNumber(); i++) {
            if (ai[i].hit.equals("")) continue;
            hitSet.add(ai[i].hit);
        }
        return hitSet.toArray(new String[hitSet.size()]);
    }
    
    /**
     * Return name of query
     * @param index
     * @return 
     */
    public String getQuery (int index) {
        return ai[index].query;
    }
    
    /**
     * Return name of hit, return empty String if the query does not match 
     * @param index
     * @return 
     */
    public String getHit (int index) {
        return ai[index].hit;
    }
    
    /**
     * Return start position of the alignment, return Integer.MIN_VALUE if the query does not match
     * @param index
     * @return 
     */
    public int getStartPos (int index) {
        return ai[index].startPos;
    }
    
    /**
     * Return end position of the alignment, return Integer.MIN_VALUE if the query does not match
     * @param index
     * @return 
     */
    public int getEndPos (int index) {
        return ai[index].endPos;
    }
    
    /**
     * Return if a query does not match the database
     * @param query
     * @return 
     */
    public boolean isNoMatch (String query) {
        int index = this.getAlignmentStartIndexByQuery(query);
        if (this.getHit(index).equals("")) return true;
        return false;
    }
    
    /**
     * Return if a query is aligned
     * @param index
     * @return 
     */
    public boolean isMatch (int index) {
        if (ai[index].hit.isEmpty()) return false;
        return true;
    }
    
    /**
     * Return if this particular alignment is a perfect match
     * @param index
     * @return 
     */
    public boolean isPerfectMatch (int index) {
        return ai[index].ifPerfectMatch;
    }
    
    /**
     * Return if a query has only perfect match (there is no the second perfect match)
     * @param query
     * @return 
     */
    public boolean isOnlyPerfectMatch (String query) {
        if (this.sortType != 0) {
            System.out.println("Alignment should be sorted by query first. Program quits");
            System.exit(0);
        }
        int startIndex = this.getAlignmentStartIndexByQuery(query);
        int endIndex = this.getAlignmentEndIndexByQuery(query);
        int cnt = 0;
        for (int i = startIndex; i < endIndex; i++ ) {
            if (this.isPerfectMatch(i)) cnt++;
            if (cnt > 1) return false;
        }
        if (cnt == 0) return false;
        return true;
    }
    
    /**
     * Return if a query has unique perfect match (there is no the second match at all)
     * @param query
     * @return 
     */
    public boolean isUniquePerfectMatch (String query) {
        if (this.sortType != 0) {
            System.out.println("Alignment should be sorted by query first. Program quits");
            System.exit(0);
        }
        int startIndex = this.getAlignmentStartIndexByQuery(query);
        if (!this.isPerfectMatch(startIndex)) return false;
        int endIndex = this.getAlignmentEndIndexByQuery(query);
        if (endIndex == startIndex + 1) return true;
        return false;
    }
    
    /**
     * Return mapping quality score of alignment, return Short.MIN_VALUE if the query does not match
     * @return 
     */
    public short getMappingQuality (int index) {
        return ai[index].mappingQuality;
    }
    
    /**
     * Return length of matched range
     * @param index
     * @return 
     */
    public short getMatchNumber (int index) {
        return ai[index].matchNumber;
    }
    
    /**
     * Return number of mismatch within matched range
     * @param index
     * @return 
     */
    public short getEditDistance (int index) {
        return ai[index].editDistance;
    }
    
    /**
     * Return the ratio of edit distance vs length of matched range
     * @param index
     * @return NaN when there is no match
     */
    public float getEditDistanceRatio (int index) {
        if (this.getMatchNumber(index) > -1 && this.getEditDistance(index) > -1) return (float)this.getEditDistance(index)/(float)this.getMatchNumber(index);
        return Float.NaN;
    }
    
    /**
     * Return strand of alignment, return Byte.MIN_VALUE if the query does not match
     * @param index
     * @return 
     */
    public byte getStrand (int index) {
        return ai[index].strand;
    }
    
    /**
     * Return sortType
     * @return 
     */
    public byte getSortType () {
        return this.sortType;
    }
    
    /**
     * Generate a subset of SimpleAlignment by reference
     * @param keepList 
     */
    public void toSubset (boolean[] keepList) {
        int cnt = 0;
        for (int i = 0; i < keepList.length; i++) {
            if (keepList[i]) cnt++;
        }
        AlignmentInfo[] sub = new AlignmentInfo[cnt];
        cnt = 0;
        for (int i = 0; i < keepList.length; i++) {
            if (keepList[i] == false) continue;
            sub[cnt] = ai[i];
            cnt++;
        }
        ai = sub;
    }
    
    /**
     * Read from blast result, outfmt 7
     * @param inputFileS
     * @param eThresh
     * @param readLength read length of fasta or fastq
     */
    public void readFromBlast (String inputFileS, double eThresh, int readLength) {
        ArrayList<AlignmentInfo> aList = new ArrayList();
        System.out.println("Reading blast alignment (outfmt 7) from: " + inputFileS);
        try {
            BufferedReader br = IOUtils.getTextReader(inputFileS);
            String temp = null;
            String query = null;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("# Query:")) {
                    query = temp.split(" ")[2];
                    temp = br.readLine();
                    temp = br.readLine();
                    if (temp.startsWith("# 0")) {
                        AlignmentInfo ai = new AlignmentInfo (query, "");
                        aList.add(ai);
                    }
                    else {
                        temp = br.readLine();
                        int cnt = 0;
                        while (!(temp = br.readLine()).startsWith("#")) {
                            AlignmentInfo aInfo = this.getAlignmentInfoFromBlast(temp, eThresh, readLength);
                            if (aInfo == null) {
                                if (cnt == 0) aList.add(new AlignmentInfo(query, ""));
                            }
                            else {
                                aList.add(aInfo);
                            }
                            cnt++;
                        }
                    }
                }
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        ai = aList.toArray(new AlignmentInfo[aList.size()]);
        this.sortByQuery();
        System.out.println("ShortreadAlignment has "+String.valueOf(ai.length) + " alignments");
    }
    
    /**
     * Make a SimpleAlignment instance from Bowtie2 sam file, support text and textGzip format
     * @param inputFileS 
     */
    public void readFromBowtie2 (String inputFileS) {
        ArrayList<AlignmentInfo> aList = new ArrayList();
        System.out.println("Reading SAM format alignment (Bowtie2) from: " + inputFileS);
        try {
            BufferedReader br;
            if (inputFileS.endsWith(".gz")) {
                br = IOUtils.getTextGzipReader(inputFileS);
            } else {
                br = br = IOUtils.getTextReader(inputFileS);
            }
            while(br.readLine().startsWith("@PG")==false) {};
            String inputStr = null;
            int cnt = 0;
            while((inputStr = br.readLine())!=null) {
                AlignmentInfo aInfo = this.getAlignmentInfoFromBowtie2(inputStr);
                aList.add(aInfo);
                if (cnt%500000 == 0) System.out.println("Read in " + String.valueOf(cnt) + " lines");
                cnt++;
            }
            ai = aList.toArray(new AlignmentInfo[aList.size()]);
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        this.sortByQuery();
        System.out.println("ShortreadAlignment has "+String.valueOf(ai.length) + " alignments");
    }
    
    /**
     * Return bit of integer, bit by bit
     * @param n input integer
     * @param k index of bit, starting from 0, left to right
     * @return 
     */
    int getBitFromInt (int n, int k) {
        return (n >> k) & 1;
    }
    
    AlignmentInfo getAlignmentInfoFromBlast (String inputStr, double eTresh, int readLength) {
        String[] temp =inputStr.split("\\s+");
        String hit = "";
        byte strand = Byte.MIN_VALUE;
        int startPos = Integer.MIN_VALUE;
        int endPos = Integer.MIN_VALUE;
        short score = Short.MIN_VALUE;
        short matchNumber = Short.MIN_VALUE;
        short editDistance = Short.MIN_VALUE;
        boolean ifPerfectMatch = false;
        if (Double.valueOf(temp[10]) > eTresh) {
            return null;
        }
        else {
            hit = temp[1];
            startPos = Integer.valueOf(temp[8]);
            endPos = Integer.valueOf(temp[9]);
            if (startPos <= endPos) {
                strand = 1;
            }
            else {
                strand = -1;
                int tempInt = startPos;
                startPos = endPos;
                endPos = tempInt;
            }
            matchNumber = (short)(Integer.valueOf(temp[3]) - Integer.valueOf(temp[5]));
            editDistance = Short.valueOf(temp[4]);
            double d = Double.valueOf(temp[11]);
            score = (short)d;
            if (Integer.valueOf(temp[3]) == readLength && Double.valueOf(temp[2]) == 100) ifPerfectMatch = true;
            AlignmentInfo aInfo = new AlignmentInfo(temp[0], hit, startPos, endPos, strand, score, matchNumber, editDistance, ifPerfectMatch);
            return aInfo;
        }
    }
    
    AlignmentInfo getAlignmentInfoFromBowtie2 (String inputStr) {
        String[] temp =inputStr.split("\\s");
        int flag=Integer.parseInt(temp[1]);
        String hit = "";
        byte strand = Byte.MIN_VALUE;
        int startPos = Integer.MIN_VALUE;
        int endPos = Integer.MIN_VALUE;
        short score = Short.MIN_VALUE;
        short matchNumber = Short.MIN_VALUE;
        short editDistance = Short.MIN_VALUE;
        boolean ifPerfectMatch = false;
        if (temp[5].equals("*")) {
            //continue;
        }
        else {
            hit = temp[2];
            if (this.getBitFromInt(flag, 4) == 1)  strand = -1;
            else  strand = 1;
            startPos = Integer.parseInt(temp[3]);
            endPos = SAMUtils.getEndPos(temp[5], startPos);
            score = Short.valueOf(temp[4]);
            matchNumber = this.getMatchNumberFromCigar(temp[5]);
            if (temp[17].startsWith("NM")) editDistance = Short.parseShort(temp[17].split(":")[2]);
            else    editDistance = Short.parseShort(temp[16].split(":")[2]);
            if (temp[5].matches("\\d+M") && editDistance == 0) ifPerfectMatch = true;
        }
        AlignmentInfo aInfo = new AlignmentInfo(temp[0], hit, startPos, endPos, strand, score, matchNumber, editDistance, ifPerfectMatch);
        return aInfo;
    }
    
    /**
     * Make a SimpleAlignment instance from BWAMEM sam file, supporting ".txt" or ".txt.gz" format
     * @param inputFileS 
     */
    public void readFromBWAMEM (String inputFileS) {
        ArrayList<AlignmentInfo> aList = new ArrayList();
        System.out.println("Reading SAM format alignment (BWA-MEM) from: " + inputFileS);
        try {
            BufferedReader br;
            if (inputFileS.endsWith(".gz")) {
                br = IOUtils.getTextGzipReader(inputFileS);
            } else {
                br = br = IOUtils.getTextReader(inputFileS);
            }
            while(br.readLine().startsWith("@PG")==false) {};
            String inputStr = null;
            int cnt = 0;
            while((inputStr = br.readLine())!=null) {
                AlignmentInfo aInfo = this.getAlignmentInfoFromBWAMEM(inputStr);
                aList.add(aInfo);
                if (cnt%500000 == 0) System.out.println("Read in " + String.valueOf(cnt) + " lines");
                cnt++;
            }
            ai = aList.toArray(new AlignmentInfo[aList.size()]);
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        this.sortByQuery();
        System.out.println("SimpleAlignment has "+String.valueOf(ai.length) + " alignments");
    }
    
    AlignmentInfo getAlignmentInfoFromBWAMEM (String[] inputStrArray) {
        String[] temp = inputStrArray;
        int orientiation=Integer.parseInt(temp[1]);
        String hit = "";
        byte strand = Byte.MIN_VALUE;
        int startPos = Integer.MIN_VALUE;
        int endPos = Integer.MIN_VALUE;
        short score = Short.MIN_VALUE;
        short matchNumber = Short.MIN_VALUE;
        short editDistance = Short.MIN_VALUE;
        boolean ifPerfectMatch = false;
        if (temp[5].equals("*")) {
            //continue;
        }
        else {
            int NMIndex = this.getNMIndex(temp);
            hit = temp[2];
            if (this.getBitFromInt(orientiation, 4) == 1) strand = -1;
            else strand = 1;
            startPos = Integer.parseInt(temp[3]);
            endPos = SAMUtils.getEndPos(temp[5], startPos);
            score = Short.valueOf(temp[4]);
            matchNumber = this.getMatchNumberFromCigar(temp[5]);
            if (NMIndex > -1) editDistance = Short.valueOf(temp[NMIndex].split(":")[2]);
            if (temp[5].matches("\\d+M") && editDistance == 0) ifPerfectMatch = true;
        }
        AlignmentInfo aInfo = new AlignmentInfo(temp[0], hit, startPos, endPos, strand, score, matchNumber, editDistance, ifPerfectMatch);
        return aInfo;
    }
    
    AlignmentInfo getAlignmentInfoFromBWAMEM (String inputStr) {
        String[] temp =inputStr.split("\\s");
        return this.getAlignmentInfoFromBWAMEM(temp);
    }
    
    private int getNMIndex (String temp[]) {
        for (int i = 11; i < temp.length; i++) {
            if (temp[i].startsWith("NM")) {
                return i;
            }
        }
        return Integer.MIN_VALUE;
    }
    
    private short getMatchNumberFromCigar (String cigar) {
        short cnt = 0;
        if (cigar.contains("M")) {
            String[] temp = cigar.split("M");
            int n;
            if (cigar.endsWith("M")) n = temp.length;
            else n = temp.length-1;
            for (int i = 0; i < n; i++) {
                String[] tem = temp[i].split("\\D");
                cnt+=Short.valueOf(tem[tem.length-1]);
            }
            return cnt;
        }
        return Short.MIN_VALUE;
    }
    
    private void readBinarySimpleAlignment (String alignmentFileS) {
        try {
            DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(alignmentFileS), 65536));
            this.sortType = dis.readByte();
            this.ai = new AlignmentInfo[dis.readInt()];
            for (int i = 0; i < ai.length; i++) {
                ai[i] = new AlignmentInfo(dis.readUTF(), dis.readUTF(), dis.readInt(), dis.readInt(), dis.readByte(), dis.readShort(), dis.readShort(), dis.readShort(), dis.readBoolean());
                if (i%5000000 == 0) System.out.println("Read in " + String.valueOf(i) + " alignments");
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    /**
     * Read SimpleAlignment from HR format
     * @param alignmentFileS 
     */
    private void readTextSimpleAlignment (String alignmentFileS) {
        try {
            BufferedReader br = new BufferedReader (new FileReader(alignmentFileS), 65536);
            String temp;
            br.readLine();
            int cnt = 0;
            ArrayList<AlignmentInfo> aList = new ArrayList();
            while ((temp = br.readLine()) != null) {
                String[] tem = temp.split("\t");
                AlignmentInfo info = new AlignmentInfo(tem[0], tem[1], Integer.valueOf(tem[2]), Integer.valueOf(tem[3]), Byte.valueOf(tem[4]), Short.valueOf(tem[5]), Short.valueOf(tem[6]), Short.valueOf(tem[7]), tem[8].equals("1")?true:false);
                aList.add(info);
                if (cnt%100000 == 0) System.out.println(String.valueOf(cnt+1) + " is read from anlignmentFileS");
                cnt++;
            }
            ai = aList.toArray(new AlignmentInfo[aList.size()]);
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
     * Output SimpleAlignment file
     * @param outfileS
     * @param f 
     */
    public void writeSimpleAlignment (String outfileS, IOFileFormat f) {
        if (f == IOFileFormat.Text) this.writeTextSimpleAlignment(outfileS);
        else if (f == IOFileFormat.Binary) this.writeBinarySimpleAlignment(outfileS);
        else throw new UnsupportedOperationException("Not supported yet.");
        System.out.println(String.valueOf(this.getAlignmentNumber()) + " alignments output to " + outfileS);
    }
    
    private void writeBinarySimpleAlignment (String outfileS) {
        try {
            DataOutputStream dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfileS), 65536));
            dos.writeByte(this.sortType);
            dos.writeInt(this.getAlignmentNumber());
            for (int i = 0; i < this.getAlignmentNumber(); i++) {
                dos.writeUTF(this.getQuery(i));
                dos.writeUTF(this.getHit(i));
                dos.writeInt(this.getStartPos(i));
                dos.writeInt(this.getEndPos(i));
                dos.writeByte(this.getStrand(i));
                dos.writeShort(this.getMappingQuality(i));
                dos.writeShort(this.getMatchNumber(i));
                dos.writeShort(this.getEditDistance(i));
                dos.writeBoolean(this.isPerfectMatch(i));
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
     * Write alignment in HR format
     * @param outfileS 
     */
    private void writeTextSimpleAlignment (String outfileS) {
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(outfileS), 65536);
            bw.write("Query\tChromosome\tStartPos\tEndPos\tStrand\tMappingQuality\tMatchNumber\tEditDistance\tIfPerfectMatch");
            bw.newLine();
            for (int i = 0; i < this.getAlignmentNumber(); i++) {
                bw.write(ai[i].toString());
                bw.newLine();
                if (i%1000000 == 0) System.out.println("wrote " + String.valueOf(i) + " alignments");
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    /**
     * Sort alignments by query name
     */
    public void sortByQuery () {
        System.out.println("Start sorting by query");
        sortType = 0;
        Arrays.sort(ai);
        System.out.println("Finished sort");
    }
    
    /**
     * Sort alignments by hit and position on hit
     */
    public void sortByHitAndPos () {
        System.out.println("Start sorting by hit and pos");
        sortType = 1;
        Arrays.sort(ai);
        System.out.println("Finished sort");
    }
    
    private class SortByHitPos implements Comparator <AlignmentInfo> {
        @Override
        public int compare(AlignmentInfo o1, AlignmentInfo o2) {
            if (o1.hit.equals(o2.hit)) {
                return o1.startPos - o2.startPos;
            }
            else {
                return o1.hit.compareTo(o2.hit);
            }
        }
        
    }
    
    class AlignmentInfo implements Comparable <AlignmentInfo>{
        /**Query sequence*/
        String query = null;
        /**Chromosome*/
        String hit = "";
        /**Reference starting position of query sequence, inclusive*/
        int startPos = Integer.MIN_VALUE;
        /**Reference ending position of query sequence, inclusive*/
        int endPos = Integer.MIN_VALUE;
        /**Strand of alignment, + = 1, - = -1*/
        byte strand = Byte.MIN_VALUE;
        /**Mapping quality*/
        short mappingQuality = Short.MIN_VALUE;
         /**The length of matched range*/
        short matchNumber = Short.MIN_VALUE;
        /**Mismatch number in the matched range*/
        short editDistance = Short.MIN_VALUE;
        /**If query is a perfect match to reference*/
        boolean ifPerfectMatch = false;
        
        AlignmentInfo (String query, String hit) {
            this.query = query;
            this.hit = hit;
        }
        
        AlignmentInfo (String hit, int pos) {
            this.hit = hit;
            this.startPos = pos;
        }
        
        AlignmentInfo (String query, String hit, int startPos, int endPos, byte strand, short mappingQuality, short matchNumber, short editDistance, boolean ifPerfectMatch) {
            this.query = query;
            this.hit = hit;
            this.startPos = startPos;
            this.endPos = endPos;
            this.strand = strand;
            this.mappingQuality = mappingQuality;
            this.matchNumber = matchNumber;
            this.editDistance = editDistance;
            this.ifPerfectMatch = ifPerfectMatch;
        }
        
        @Override
        public int compareTo(AlignmentInfo o) {
            if (sortType == 0) {
                return query.compareTo(o.query);
            }
            else if (sortType == 1) {
                if (hit.equals(o.hit)) {
                    return startPos - o.startPos;
                }
                else {
                    return hit.compareTo(o.hit);
                }
            }
            else return 0;
        }
        
        @Override
        public String toString () {
            StringBuilder sb = new StringBuilder();
            sb.append(query).append("\t").append(hit).append("\t").append(startPos).append("\t").append(endPos).append("\t").append(strand).append("\t").append(mappingQuality).append("\t");
            sb.append(matchNumber).append("\t").append(editDistance).append("\t");
            if (ifPerfectMatch) sb.append("1");
            else sb.append("0");
            return sb.toString();
        }
    }
}
