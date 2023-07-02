/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.app.grt;

import pgl.infra.align.g2.SAMUtils;
import pgl.infra.dna.allele.AlleleEncoder;
import pgl.infra.dna.snp.SNPOld;
import pgl.infra.pos.ChrPos;
import gnu.trove.list.array.TByteArrayList;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

/**
 *
 * @author feilu
 */
public class TagAnnotations {
    int tagLengthInLong = Integer.MIN_VALUE;
    int offSet = 8;
    int groupIdentifierLength = 6;
    boolean ifSorted = false;
    int groupCount = -1;
    List<TagAnnotation> taList = null;
    
    public TagAnnotations (String infileS) {
        if (infileS.endsWith(".tp")) {
            this.readTPBinaryFile(infileS);
        }
        else if (infileS.endsWith(".tas")) {
            this.readBinaryFile(infileS);
        }
    }
    
    private void readTPBinaryFile (String infileS) {
        System.out.println("Reading TagAnnotations file from " + infileS);
        try {
            DataInputStream dis = IOUtils.getBinaryReader(infileS);
            this.tagLengthInLong = dis.readInt();
            int currentTagNum = dis.readInt();
            groupCount = (int)Math.pow(4, groupIdentifierLength);
            taList = new ArrayList<>();
            for (int i = 0; i < groupCount; i++) {
                TagAnnotation ta = new TagAnnotation(i);
                taList.add(ta);
            }
            if (currentTagNum == -1) currentTagNum = (int)((new File(infileS).length()-8)/(tagLengthInLong*2*8+2+4));
            for (int i = 0; i < currentTagNum; i++) {
                long[] tag = new long[2*this.tagLengthInLong];
                byte r1Len = 0;
                byte r2Len = 0;
                int readCount = 0;
                int groupIndex = 0;
                for (int j = 0; j < tag.length; j++) {
                    tag[j] = dis.readLong();
                }
                r1Len = dis.readByte();
                r2Len = dis.readByte();
                readCount = dis.readInt();
                groupIndex = TagUtils.getGroupIndexFromTag(tag, offSet, groupIdentifierLength);
                taList.get(groupIndex).appendTag(tag, r1Len, r2Len, readCount);
            }
            dis.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        this.sort();
    }
    
    public void writeFastqFile (String r1FastqFileS, String r2FastqFileS) {
        String polyQ = "????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????";
        try {
            System.out.println("Writing fastq file from TagAnnotations");
            BufferedWriter bw1 = IOUtils.getTextWriter(r1FastqFileS);
            BufferedWriter bw2 = IOUtils.getTextWriter(r2FastqFileS);
            StringBuilder sb = new StringBuilder();
            String identifier = null;
            String[] reads = null;
            long cnt = 0;
            for (int i = 0; i < this.getGroupNumber(); i++) {
                for (int j = 0; j < this.getTagNumber(i); j++) {
                    sb = new StringBuilder();
                    sb.append("@").append(i).append("_").append(j).append("_").append(this.getReadNumber(i, j));
                    identifier = sb.toString();
                    bw1.write(identifier); bw1.newLine();          
                    bw2.write(identifier); bw2.newLine();
                    reads = TagUtils.getReadsFromTag(this.getTag(i, j), this.getR1TagLength(i, j), this.getR2TagLength(i, j));
                    bw1.write(reads[0]);bw1.newLine();
                    bw2.write(reads[1]);bw2.newLine();
                    bw1.write("+");bw1.newLine();
                    bw2.write("+");bw2.newLine();
                    bw1.write(polyQ.substring(0, this.getR1TagLength(i, j)));bw1.newLine();
                    bw2.write(polyQ.substring(0, this.getR2TagLength(i, j)));bw2.newLine();
                    cnt++;
                    if (cnt%10000000 == 0) System.out.println(String.valueOf(cnt) + " tags have been converted to Fastq");
                }
            }
            bw1.flush();bw1.close();
            bw2.flush();bw2.close();
            System.out.println("Fastq files are written to " + String.valueOf(r1FastqFileS) + " " + String.valueOf(r2FastqFileS));
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
     
    public void readBinaryFile (String infileS) {
        try {
            System.out.println("Reading TagAnnotations file from " + infileS);
            DataInputStream dis = IOUtils.getBinaryReader(infileS);
            this.tagLengthInLong = dis.readInt();
            this.offSet = dis.readInt();
            this.groupIdentifierLength = dis.readInt();
            ifSorted = dis.readBoolean();
            this.groupCount = (int)Math.pow(4, groupIdentifierLength);
            taList = new ArrayList<>();
            long cnt = 0;
            for (int i = 0; i < groupCount; i++) {
                int tagNumber = dis.readInt();
                int groupIndex = dis.readInt();
                TagAnnotation ta = new TagAnnotation(this.tagLengthInLong, groupIndex, tagNumber, ifSorted);               
                for (int j = 0; j < tagNumber; j++) {
                    cnt++;
                    long[] tag = new long[this.tagLengthInLong*2];
                    for (int k = 0; k < tag.length; k++) {
                        tag[k] = dis.readLong();
                    }
                    byte r1Len = dis.readByte();
                    byte r2Len = dis.readByte();
                    int readNumber = dis.readInt();
                    short r1Chr = dis.readShort();
                    int r1Pos = dis.readInt();
                    byte r1Strand = dis.readByte();
                    byte r1MapQ = dis.readByte();
                    short r2Chr = dis.readShort();
                    int r2Pos = dis.readInt();
                    byte r2Strand = dis.readByte();
                    byte r2MapQ = dis.readByte();
                    byte snpNumber = dis.readByte();
                    List<SNPOld> tagSNPList = new ArrayList<>();
                    for (int k = 0; k < snpNumber; k++) {
                        short chr = dis.readShort();
                        int pos = dis.readInt();
                        byte ref = dis.readByte();
                        byte altNumber = dis.readByte();
                        TByteArrayList alts = new TByteArrayList();
                        for (int u = 0; u < altNumber; u++) {
                            alts.add(dis.readByte());
                        }
                        tagSNPList.add(new SNPOld(chr, pos, ref, alts));
                    }
                    List<AlleleInfo> tagAlleleList = new ArrayList();
                    byte alleleNumber = dis.readByte();
                    for (int k = 0; k < alleleNumber; k++) {
                        short chr = dis.readShort();
                        int pos = dis.readInt();
                        byte allele = dis.readByte();
                        byte base = dis.readByte();
                        byte end = dis.readByte();
                        byte relativePosition = dis.readByte();
                        AlleleInfo ai = new AlleleInfo (chr, pos, allele, base, end, relativePosition);
                        tagAlleleList.add(ai);
                    }
                    ta.appendTag(tag, r1Len, r2Len, readNumber, r1Chr, r1Pos, r1Strand, r1MapQ, r2Chr, r2Pos, r2Strand, r2MapQ, tagSNPList, tagAlleleList);
                    if (cnt%10000000 == 0) System.out.println("Reading in "+String.valueOf(cnt)+" tags");
                }
                taList.add(ta);
            }
            System.out.println(String.valueOf(cnt) + " tags are in " + infileS);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void writeBinaryFile (String outfileS) {
        System.out.println("Writing TagAnnotationss file to " + outfileS);
        try {
            DataOutputStream dos = IOUtils.getBinaryWriter(outfileS);
            dos.writeInt(this.getTagLengthInLong());
            dos.writeInt(this.getGroupIdentiferOffset());
            dos.writeInt(this.groupIdentifierLength);
            dos.writeBoolean(this.ifSorted);
            long cnt = 0;
            for (int i = 0; i < this.getGroupNumber(); i++) {
                dos.writeInt(this.getTagNumber(i));
                dos.writeInt(i);
                for (int j = 0; j < taList.get(i).getTagNumber(); j++) {
                    long[] tag = taList.get(i).getTag(j);
                    for (int k = 0; k < tag.length; k++) {
                        dos.writeLong(tag[k]);
                    }
                    dos.writeByte(this.getR1TagLength(i, j));
                    dos.writeByte(this.getR2TagLength(i, j));
                    dos.writeInt(this.getReadNumber(i, j));
                    dos.writeShort(this.getR1Chromosome(i, j));
                    dos.writeInt(this.getR1StartPosition(i, j));
                    dos.writeByte(this.getR1Strand(i, j));
                    dos.writeByte(this.getR1MapQ(i, j));
                    dos.writeShort(this.getR2Chromosome(i, j));
                    dos.writeInt(this.getR2StartPosition(i, j));
                    dos.writeByte(this.getR2Strand(i, j));
                    dos.writeByte(this.getR2MapQ(i, j));
                    byte snpNumber = this.getSNPNumberOfTag(i, j);
                    dos.writeByte(snpNumber);
                    List<SNPOld> snpList = this.getSNPOfTag(i, j);
                    for (int k = 0; k < snpNumber; k++) {
                        SNPOld s = snpList.get(k);
                        dos.writeShort(s.getChromosome());
                        dos.writeInt(s.getPosition());
                        dos.writeByte(s.getRefAlleleByte());  
                        dos.writeByte(s.getAltAlleleNumber());
                        for (int u = 0; u < s.getAltAlleleNumber(); u++) {
                            dos.writeByte(s.getAltAlleleByte(u));
                        }
                    }
                    List<AlleleInfo> alleleList = this.getAlleleOfTag(i, j);
                    byte alleleNumber = (byte)alleleList.size();
                    dos.writeByte(alleleNumber);
                    for (int k = 0; k < alleleNumber; k++) {
                        dos.writeShort(alleleList.get(k).getChromosome());
                        dos.writeInt(alleleList.get(k).getPosition());
                        dos.writeByte(alleleList.get(k).getAllele());
                        dos.writeByte(alleleList.get(k).getBase());
                        dos.writeByte(alleleList.get(k).getEnd());
                        dos.writeByte(alleleList.get(k).getRelativePosition());
                    }
                    cnt++;
                    if (cnt%10000000 == 0) System.out.println(String.valueOf(cnt) + " tags ouput to " + outfileS);
                }
            }
            dos.flush();
            dos.close();
            System.out.println(String.valueOf(cnt) + " tags are written to " + outfileS);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void writeTextFileOfGroup (String outfileS, int groupIndex) {
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("TagLengthInLong:\t" + String.valueOf(this.getTagLengthInLong()));
            bw.newLine();
            bw.write("GroupIdentifierOffset:\t" + String.valueOf(this.getGroupIdentiferOffset()));
            bw.newLine();
            bw.write("GroupIdentifierLength:\t" + String.valueOf(this.getGroupIdentiferLength()));
            bw.newLine();
            for (int i = 0; i < this.getGroupNumber(); i++) {
                if (i != groupIndex) continue;
                StringBuilder sb = new StringBuilder();
                sb.append("GroupIndex:\t").append(i).append("\nTagNumber:\t").append(this.getTagNumber(i));
                bw.write(sb.toString());
                bw.newLine();
                for (int j = 0; j < this.getTagNumber(i); j++) {
                    sb = new StringBuilder();
                    sb.append("Tag_ID\t:").append(i).append("_").append(j).append("\n");
                    sb.append(this.getR1TagLength(i, j)).append("\t").append(this.getR2TagLength(i, j)).append("\t").append(this.getReadNumber(i, j)).append("\n");
                    String[] reads = TagUtils.getReadsFromTag(this.getTag(i, j), this.getR1TagLength(i, j), this.getR2TagLength(i, j));
                    sb.append(reads[0]).append("\t").append(reads[1]).append("\n");
                    sb.append("R1 alignment:\t").append(this.getR1Chromosome(i, j)).append("\t").append(this.getR1StartPosition(i, j)).append("\t").append(this.getR1Strand(i, j)).append("\t").append(this.getR1MapQ(i, j)).append("\n");
                    sb.append("R2 alignment:\t").append(this.getR2Chromosome(i, j)).append("\t").append(this.getR2StartPosition(i, j)).append("\t").append(this.getR2Strand(i, j)).append("\t").append(this.getR2MapQ(i, j)).append("\n");
                    List<SNPOld> snpList = this.getSNPOfTag(i, j);
                    sb.append("SNPNumber:\t").append(snpList.size());
                    sb.append("\nSNPs:");
                    for (int k = 0; k < snpList.size(); k++) {
                        SNPOld s = snpList.get(k);
                        sb.append("\t|").append(k).append("->").append(s.getChromosome()).append("\t").append(s.getPosition()).append("\t").append(s.getRefAllele()).append("\t");
                        for (int u = 0; u < s.getAltAlleleNumber(); u++) {
                            sb.append(s.getAltAllele(u)).append("/");
                        }
                        sb.deleteCharAt(sb.length()-1);
                    }
                    sb.append("\n");
                    List<AlleleInfo> alleleList = this.getAlleleOfTag(i, j);
                    sb.append("AlleleNumber:\t").append(alleleList.size());
                    sb.append("\nAlleles:");
                    for (int k = 0; k < alleleList.size(); k++) {
                        sb.append("\t|").append(k).append("-").append(alleleList.get(k).getEnd()).append("-").append(alleleList.get(k).getRelativePosition())
                                .append("->").append(alleleList.get(k).getChromosome()).append("\t").append(alleleList.get(k).getPosition())
                                .append("\t").append(AlleleEncoder.alleleCodingToBaseMap.get(alleleList.get(k).getAllele()));
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                } 
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void writeTextFile (String outfileS) {
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("TagLengthInLong:\t" + String.valueOf(this.getTagLengthInLong()));
            bw.newLine();
            bw.write("GroupIdentifierOffset:\t" + String.valueOf(this.getGroupIdentiferOffset()));
            bw.newLine();
            bw.write("GroupIdentifierLength:\t" + String.valueOf(this.getGroupIdentiferLength()));
            bw.newLine();
            for (int i = 0; i < this.getGroupNumber(); i++) {
                StringBuilder sb = new StringBuilder();
                sb.append("GroupIndex:\t").append(i).append("\nTagNumber:\t").append(this.getTagNumber(i));
                bw.write(sb.toString());
                bw.newLine();
                for (int j = 0; j < this.getTagNumber(i); j++) {
                    sb = new StringBuilder();
                    sb.append("Tag_ID\t:").append(i).append("_").append(j).append("\n");
                    sb.append(this.getR1TagLength(i, j)).append("\t").append(this.getR2TagLength(i, j)).append("\t").append(this.getReadNumber(i, j)).append("\n");
                    String[] reads = TagUtils.getReadsFromTag(this.getTag(i, j), this.getR1TagLength(i, j), this.getR2TagLength(i, j));
                    sb.append(reads[0]).append("\t").append(reads[1]).append("\n");
                    sb.append("R1 alignment:\t").append(this.getR1Chromosome(i, j)).append("\t").append(this.getR1StartPosition(i, j)).append("\t").append(this.getR1Strand(i, j)).append("\t").append(this.getR1MapQ(i, j)).append("\n");
                    sb.append("R2 alignment:\t").append(this.getR2Chromosome(i, j)).append("\t").append(this.getR2StartPosition(i, j)).append("\t").append(this.getR2Strand(i, j)).append("\t").append(this.getR2MapQ(i, j)).append("\n");
                    List<SNPOld> snpList = this.getSNPOfTag(i, j);
                    sb.append("SNPNumber:\t").append(snpList.size());
                    sb.append("\nSNPs:");
                    for (int k = 0; k < snpList.size(); k++) {
                        SNPOld s = snpList.get(k);
                        sb.append("\t|").append(k).append("->").append(s.getChromosome()).append("\t").append(s.getPosition()).append("\t").append(s.getRefAllele()).append("\t");
                        for (int u = 0; u < s.getAltAlleleNumber(); u++) {
                            sb.append(s.getAltAllele(u)).append("/");
                        }
                        sb.deleteCharAt(sb.length()-1);
                    }
                    sb.append("\n");
                    List<AlleleInfo> alleleList = this.getAlleleOfTag(i, j);
                    sb.append("AlleleNumber:\t").append(alleleList.size());
                    sb.append("\nAlleles:");
                    for (int k = 0; k < alleleList.size(); k++) {
                        sb.append("\t|").append(k).append("-").append(alleleList.get(k).getEnd()).append("-").append(alleleList.get(k).getRelativePosition())
                                .append("->").append(alleleList.get(k).getChromosome()).append("\t").append(alleleList.get(k).getPosition())
                                .append("\t").append(AlleleEncoder.alleleCodingToBaseMap.get(alleleList.get(k).getAllele()));
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                } 
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public List<ChrPos> filterTagAnnotationsWithValidatedGenotype (String genotypeDirS) {
        File[] fs = new File(genotypeDirS).listFiles();
        File[] gzfs = IOUtils.listFilesEndsWith(fs, ".gz");
        File[] vcffs = IOUtils.listFilesEndsWith(fs, ".vcf");
        List<File> fList = new ArrayList();
        fList.addAll(Arrays.asList(gzfs));
        fList.addAll(Arrays.asList(vcffs));
        List<ChrPos> validatedSNPPosList = new ArrayList<>();
        for (int i = 0; i < fList.size(); i++) {
            try {
                BufferedReader br = null;
                if (fList.get(i).getName().endsWith(".gz")) {
                    br = IOUtils.getTextGzipReader(fList.get(i).getAbsolutePath());
                }
                else {
                    br = IOUtils.getTextReader(fList.get(i).getAbsolutePath());
                }
                String temp = null;
                while ((temp = br.readLine()).startsWith("##")) {}
                while ((temp = br.readLine()) != null) {
                    List<String> l = PStringUtils.fastSplit(temp.substring(0, 50));
                    validatedSNPPosList.add(new ChrPos(Short.parseShort(l.get(0)), Integer.parseInt(l.get(1))));
                }
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }
        Collections.sort(validatedSNPPosList);
        for (int i = 0; i < this.getGroupNumber(); i++) {
            for (int j = 0; j < this.getTagNumber(i); j++) {
                List<SNPOld> snpList = new ArrayList<>();
                List<SNPOld> currentSNPList = this.getSNPOfTag(i, j);
                for (int k = 0; k < currentSNPList.size(); k++) {
                    ChrPos query = currentSNPList.get(k).getChrPos();
                    int index = Collections.binarySearch(validatedSNPPosList, query);
                    if (index < 0) continue;
                    snpList.add(currentSNPList.get(k));
                }
                this.setSNPOfTag(i, j, snpList);
                List<AlleleInfo> alleleList = new ArrayList<>();
                List<AlleleInfo> currentAlleleList = this.getAlleleOfTag(i, j);
                for (int k = 0; k < currentAlleleList.size(); k++) {
                    ChrPos query = currentAlleleList.get(k).getChrPos();
                    int index = Collections.binarySearch(validatedSNPPosList, query);
                    if (index < 0) continue;
                    alleleList.add(currentAlleleList.get(k));
                }
                this.setAlleleOfTag(i, j, alleleList);
            }
        }
        return validatedSNPPosList;
    }
    
    public void callAllele (String samFileS, SNPCounts sc, int mapQThresh, int maxMappingIntervalThresh) {
        System.out.println("Start adding alleles to DB");
        String temp = null;
        try {
            BufferedReader br = IOUtils.getTextGzipReader(samFileS);
            while ((temp = br.readLine()).startsWith("@SQ")){}
            int queryCount = 0;
            List<AlleleInfo> tagAlleleList = new ArrayList();
            long cnt = 0;
            long snpCnt = 0;
            while ((temp = br.readLine()) != null) {
                List<String> l = SAMUtils.getAlignElements(temp);            
                if (Integer.parseInt(l.get(1)) > 2000) continue; //remove supplement alignment to have a pair of alignments for PE reads
                queryCount++;
                List<AlleleInfo> alleles = SAMUtils.getAlleles3(l, mapQThresh, sc, queryCount);
                if (queryCount == 1) {
                    if (alleles != null) {
                        tagAlleleList.addAll(alleles);
                    }
                }
                else if (queryCount == 2) {
                    if (l.get(6).equals("=")) {
                        double len = Math.abs(Double.valueOf(l.get(8)));
                        if (len < maxMappingIntervalThresh) {
                            if (alleles != null) {
                                tagAlleleList.addAll(alleles);
                            }
                            cnt++;
                            if (cnt%10000000 == 0) System.out.println(String.valueOf(cnt) + " tags are properly aligned for allele calling");                           
                            if (tagAlleleList.size() != 0) {
                                List<String> ll = PStringUtils.fastSplit(l.get(0), "_");                               
                                int groupIndex = Integer.parseInt(ll.get(0));
                                int tagIndex = Integer.parseInt(ll.get(1));
                                this.setAlleleOfTag(groupIndex, tagIndex, tagAlleleList);
                                this.sortAlleleListByPosition(groupIndex, tagIndex);
                                snpCnt++;
                            }
                        }
                    }
                    queryCount = 0;
                    tagAlleleList = new ArrayList();
                }
            }
            br.close();
            System.out.println("A total of "+String.valueOf(cnt) + " tags are properly aligned for allele calling");
            System.out.println("A total of "+String.valueOf(snpCnt) + " tags have allele calls");
        }
        catch (Exception e) {
            System.out.println(temp);
            e.printStackTrace();
        }
    }
    
    public void callSNP (String samFileS, int mapQThresh, int maxMappingIntervalThresh, int maxDivergence) {
        System.out.println("Start adding raw SNPs to DB");
        maxDivergence++;
        try {
            BufferedReader br = IOUtils.getTextGzipReader(samFileS);
            String temp = null;
            while ((temp = br.readLine()).startsWith("@SQ")){}
            int queryCount = 0;
            List<SNPOld> tagSNPList = new ArrayList();
            long cnt = 0;
            long snpCnt = 0;
            while ((temp = br.readLine()) != null) {
                List<String> l = SAMUtils.getAlignElements(temp);
                List<String> ll = PStringUtils.fastSplit(l.get(0), "_");                               
                int groupIndex = Integer.parseInt(ll.get(0));
                int tagIndex = Integer.parseInt(ll.get(1));
                if (Integer.parseInt(l.get(1)) > 2000) continue; //remove supplement alignment to have a pair of alignments for PE reads
                queryCount++;
                List<SNPOld> snpList = SAMUtils.getVariants(l, mapQThresh);
                if (queryCount == 1) {
                    if (!l.get(5).startsWith("*")) {
                        short chr = Short.parseShort(l.get(2));
                        int startPos = Integer.parseInt(l.get(3));
                        boolean ifMinus = SAMUtils.isReverseAligned(Integer.parseInt(l.get(1)));
                        byte strand = 1;
                        if (ifMinus) strand = 0;
                        byte mapQ = Byte.parseByte(l.get(4));
                        this.setR1Chromosome(groupIndex, tagIndex, chr);
                        this.setR1StartPosition(groupIndex, tagIndex, startPos);
                        this.setR1Strand(groupIndex, tagIndex, strand);
                        this.setR1MapQ(groupIndex, tagIndex, mapQ);
                    }
                    if (snpList != null) tagSNPList.addAll(snpList);
                }
                else if (queryCount == 2) {
                    if (!l.get(5).startsWith("*")) {
                        short chr = Short.parseShort(l.get(2));
                        int startPos = Integer.parseInt(l.get(3));
                        boolean ifMinus = SAMUtils.isReverseAligned(Integer.parseInt(l.get(1)));
                        byte strand = 1;
                        if (ifMinus) strand = 0;
                        byte mapQ = Byte.parseByte(l.get(4));
                        this.setR2Chromosome(groupIndex, tagIndex, chr);
                        this.setR2StartPosition(groupIndex, tagIndex, startPos);
                        this.setR2Strand(groupIndex, tagIndex, strand);
                        this.setR2MapQ(groupIndex, tagIndex, mapQ);
                    }
                    if (l.get(6).equals("=")) {
                        double len = Math.abs(Double.valueOf(l.get(8)));
                        if (len < maxMappingIntervalThresh) {
                            if (snpList != null) tagSNPList.addAll(snpList);
                            cnt++;
                            if (cnt%10000000 == 0) System.out.println(String.valueOf(cnt) + " tags are properly aligned for SNP calling");                           
                            if (tagSNPList.size() > 0 && tagSNPList.size() < maxDivergence) {                                
                                Collections.sort(tagSNPList);        
                                setSNPOfTag(groupIndex, tagIndex, tagSNPList);
                                snpCnt++;
                            }
                        }
                    }
                    queryCount = 0;
                    tagSNPList = new ArrayList();
                }
            }
            br.close();
            System.out.println("A total of "+String.valueOf(cnt) + " tags are properly aligned for SNP calling");
            System.out.println("A total of "+String.valueOf(snpCnt) + " tags have SNP calls");
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public SNPCounts getSNPCounts () {
        return new SNPCounts(this);
    }
    
    public boolean addTagAnnotations (TagAnnotations ata) {
        if (this.getTagLengthInLong() != ata.getTagLengthInLong()) return false;
        if (this.getGroupIdentiferOffset() != ata.getGroupIdentiferOffset()) return false;
        if (this.getGroupIdentiferLength() != ata.getGroupIdentiferLength()) return false;
        taList.parallelStream().forEach(ta -> {
            ta.r1LenList.addAll(ata.taList.get(ta.groupIndex).r1LenList);
            ta.r2LenList.addAll(ata.taList.get(ta.groupIndex).r2LenList);
            ta.readCountList.addAll(ata.taList.get(ta.groupIndex).readCountList);
            ta.r1ChrList.addAll(ata.taList.get(ta.groupIndex).r1ChrList);
            ta.r2ChrList.addAll(ata.taList.get(ta.groupIndex).r2ChrList);
            ta.r1StartPosList.addAll(ata.taList.get(ta.groupIndex).r1StartPosList);
            ta.r2StartPosList.addAll(ata.taList.get(ta.groupIndex).r2StartPosList);
            ta.r1StrandList.addAll(ata.taList.get(ta.groupIndex).r1StrandList);
            ta.r2StrandList.addAll(ata.taList.get(ta.groupIndex).r2StrandList);
            ta.r1MapQList.addAll(ata.taList.get(ta.groupIndex).r1MapQList);
            ta.r2MapQList.addAll(ata.taList.get(ta.groupIndex).r2MapQList);
            ta.tagList.addAll(ata.taList.get(ta.groupIndex).tagList);
            ta.SNPList.addAll(ata.taList.get(ta.groupIndex).SNPList);
            ta.alleleList.addAll(ata.taList.get(ta.groupIndex).alleleList);
        });
        this.ifSorted = false;
        return true;
    }
    
    public void sortSNPListByPosition (int groupIndex, int tagIndex) {
        this.taList.get(groupIndex).sortSNPListByPosition(tagIndex);
    }
    
    public void sortAlleleListByPosition (int groupIndex, int tagIndex) {
        this.taList.get(groupIndex).sortAlleleListByPosition(tagIndex);
    }
    
    public void sort () {
        System.out.println("TagAnnotations sort begins");
        taList.parallelStream().forEach(ta -> {
            ta.sort();
        });
        System.out.println("TagAnnotations sort ends");
        this.ifSorted = true;
    }
    
    public boolean isSorted () {
        return this.ifSorted;
    }
    
    public int getMaxTagNumberAcrossGroups () {
        int max = Integer.MIN_VALUE;
        for (int i = 0; i < this.getGroupNumber(); i++) {
            if (this.getTagNumber(i) > max) max = this.getTagNumber(i);
        }
        return max;
    }
    
    public int getReadNumber (int groupIndex, int tagIndex) {
        return taList.get(groupIndex).getReadNumber(tagIndex);
    }
    
    public byte getR1TagLength (int groupIndex, int tagIndex) {
        return taList.get(groupIndex).getR1TagLength(tagIndex);
    }
    
    public byte getR2TagLength (int groupIndex, int tagIndex) {
        return taList.get(groupIndex).getR2TagLength(tagIndex);
    }
    
    public short getR1Chromosome (int groupIndex, int tagIndex) {
        return taList.get(groupIndex).getR1Chromosome(tagIndex);
    }
    
    public short getR2Chromosome (int groupIndex, int tagIndex) {
        return taList.get(groupIndex).getR2Chromosome(tagIndex);
    }
    
    public int getR1StartPosition (int groupIndex, int tagIndex) {
        return taList.get(groupIndex).getR1StartPosition(tagIndex);
    }
    
    public int getR2StartPosition (int groupIndex, int tagIndex) {
        return taList.get(groupIndex).getR2StartPosition(tagIndex);
    }
    
    public byte getR1Strand (int groupIndex, int tagIndex) {
        return taList.get(groupIndex).getR1Strand(tagIndex);
    }
    
    public byte getR2Strand (int groupIndex, int tagIndex) {
        return taList.get(groupIndex).getR2Strand(tagIndex);
    }
    
    public byte getR1MapQ (int groupIndex, int tagIndex) {
        return taList.get(groupIndex).getR1MapQ(tagIndex);
    }
    
    public byte getR2MapQ (int groupIndex, int tagIndex) {
        return taList.get(groupIndex).getR2MapQ(tagIndex);
    }
    
    public int getTagIndex (long[] tag, int groupIndex) {
        return taList.get(groupIndex).getTagIndex(tag);
    }
    
    public int getGroupIndex (long[] tag) {
        return TagUtils.getGroupIndexFromTag(tag, this.getGroupIdentiferOffset(), this.getGroupIdentiferLength());
    }
    
    public long getTotalReadNumber () {
        long sum = 0; 
        for (int i = 0; i < this.getGroupNumber(); i++) {
            sum += taList.get(i).getTotalReadNum();
        }
        return sum;
    }
    
    public int getTagNumber (int groupIndex) {
        return taList.get(groupIndex).getTagNumber();
    }
    
    public long getTagNumber () {
        long sum = 0;
        for (int i = 0; i < this.getGroupNumber(); i++) {
            sum += taList.get(i).getTagNumber();
        }
        return sum;
    }
    
    public int getGroupIdentiferLength () {
        return this.groupIdentifierLength;
    }
    
    public int getGroupIdentiferOffset () {
        return this.offSet;
    }
    
    public int getGroupNumber () {
        return this.groupCount;
    }
    
    public int getTagLengthInLong () {
        return this.tagLengthInLong;
    }
    
    public long[] getTag (int groupIndex, int tagIndex) {
        return taList.get(groupIndex).getTag(tagIndex);
    }
    
    public byte getSNPNumberOfTag (int groupIndex, int tagIndex) {
        return (byte)this.getSNPOfTag(groupIndex, tagIndex).size();
    }
    
    public List<SNPOld> getSNPOfTag (int groupIndex, int tagIndex) {
        return this.taList.get(groupIndex).getSNPOfTag(tagIndex);
    } 
    
    public byte getAlleleNumberOfTag (int groupIndex, int tagIndex) {
        return (byte)this.getAlleleOfTag(groupIndex, tagIndex).size();
    }
        
    public List<AlleleInfo> getAlleleOfTag (int groupIndex, int tagIndex) {
        return this.taList.get(groupIndex).getAlleleOfTag(tagIndex);
    }
    
    public void setR1Chromosome (int groupIndex, int tagIndex, short chr) {
        taList.get(groupIndex).setR1Chromosome(tagIndex, chr);
    }
    
    public void setR2Chromosome (int groupIndex, int tagIndex, short chr) {
        taList.get(groupIndex).setR2Chromosome(tagIndex, chr);
    }
    
    public void setR1StartPosition (int groupIndex, int tagIndex, int pos) {
        taList.get(groupIndex).setR1StartPosition(tagIndex, pos);
    }
    
    public void setR2StartPosition (int groupIndex, int tagIndex, int pos) {
        taList.get(groupIndex).setR2StartPosition(tagIndex, pos);
    }
    
    public void setR1Strand (int groupIndex, int tagIndex, byte strand) {
        taList.get(groupIndex).setR1Strand(tagIndex, strand);
    }
    
    public void setR2Strand (int groupIndex, int tagIndex, byte strand) {
        taList.get(groupIndex).setR2Strand(tagIndex, strand);
    }
    
    public void setR1MapQ (int groupIndex, int tagIndex, byte mapQ) {
        taList.get(groupIndex).setR1MapQ(tagIndex, mapQ);
    }
    
    public void setR2MapQ (int groupIndex, int tagIndex, byte mapQ) {
        taList.get(groupIndex).setR2MapQ(tagIndex, mapQ);
    }
    
    public void setSNPOfTag (int groupIndex, int tagIndex, List<SNPOld> tagSNPList) {
        taList.get(groupIndex).setSNPOfTag(tagIndex, tagSNPList);
    }
    
    public void setAlleleOfTag (int groupIndex, int tagIndex, List<AlleleInfo> tagAlleleList) {
        taList.get(groupIndex).setAlleleOfTag(tagIndex, tagAlleleList);
    }
    
    public void removeSNPOfTag (int groupIndex, int tagIndex) {
        taList.get(groupIndex).removeSNPOfTag(tagIndex);
    }
    
    public void removeAlleleOfTag (int groupIndex, int tagIndex) {
        taList.get(groupIndex).removeAlleleOfTag(tagIndex);
    }
    
    public void removeAllSNP () {
        for (int i = 0; i < this.getGroupNumber(); i++) {
            for (int j = 0; j < this.getTagNumber(i); j++) {
                this.removeSNPOfTag(i, j);
            }
        }
    }
    
    public void removeAllAllele () {
        for (int i = 0; i < this.getGroupNumber(); i++) {
            for (int j = 0; j < this.getTagNumber(i); j++) {
                this.removeAlleleOfTag(i, j);
            }
        }
    }
    
    public void collapseCounts (int minReadCount) {
        System.out.println("Start collapsing read counts of TagAnnotations with "+this.getTagNumber()+" tags.");
        AtomicInteger acnt = new AtomicInteger();
        AtomicInteger gcnt = new AtomicInteger();
        AtomicInteger pcnt = new AtomicInteger();
        if (this.isSorted() == false) this.sort();
        NumberFormat defaultFormat = NumberFormat.getPercentInstance();
        defaultFormat.setMinimumFractionDigits(1);
	int step = (int)(this.getGroupNumber()*0.2);
        taList.parallelStream().forEach(ta -> {
            int cnt = ta.collapseCounts(minReadCount);
            acnt.addAndGet(cnt);
            int count = gcnt.addAndGet(1);            
            if (count%step == 0) {
                System.out.println("Colapsed " + defaultFormat.format(0.1*pcnt.addAndGet(2)));
            }
        });     
        System.out.println("Collapsing tags complected. Tag rows collapsed after sorting: " + acnt.get());
    }
}
