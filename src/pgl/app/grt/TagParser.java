/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.app.grt;

import com.koloboke.collect.map.hash.HashByteByteMap;
import pgl.infra.dna.BaseEncoder;
import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PArrayUtils;

/**
 *
 * @author feilu
 */
public class TagParser {
    LibraryInfo li = null;
    int tagLengthInLong = 3;
    String polyA = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    int setReadLength = tagLengthInLong*BaseEncoder.longChunkSize;
    int numThreads = 32;
    int minReadLength = 16;
    
    public TagParser (LibraryInfo li) {
        this.li = li;
    }
    
    public void setThreads (int numThreads) {
        this.numThreads = numThreads;
    }
    
    public void parseFastq(String tagBySampleDirS) {
        String[] libs = li.getLibArray();
        int[][] indices = PArrayUtils.getSubsetsIndicesBySubsetSize(libs.length, this.numThreads);
        for (int i = 0; i < indices.length; i++) {
            Integer[] subLibIndices = new Integer[indices[i][1]-indices[i][0]];
            for (int j = 0; j < subLibIndices.length; j++) {
                subLibIndices[j] = indices[i][0]+j;
            }
            List<Integer> indexList = Arrays.asList(subLibIndices);
            indexList.parallelStream().forEach(index -> {
                String fastqR1 = li.getFastqFileSR1(index);
                String fastqR2 = li.getFastqFileSR2(index);
                HashMap<String, Set<String>> barcodeR1TaxaMap = li.getbarcodeR1TaxaMap(index);
                HashMap<String, Set<String>> barcodeR2TaxaMap = li.getbarcodeR2TaxaMap(index);
                String[] taxaNames = li.getTaxaNames(index);
                String cutter1 = li.getCutter1();
                String cutter2 = li.getCutter2();
                this.splitFastq(fastqR1, fastqR2, barcodeR1TaxaMap, barcodeR2TaxaMap, taxaNames, tagBySampleDirS, cutter1, cutter2);
            });
        }
    }
    
    private void splitFastq (String fastqR1, String fastqR2, 
            HashMap<String, Set<String>> barcodeR1TaxaMap, HashMap<String, Set<String>> barcodeR2TaxaMap, String[] taxaNames, String tagBySampleDirS, String cutter1, String cutter2) {
        HashMap<String, DataOutputStream> taxaWriterMap = new HashMap<>();
        DataOutputStream[] doss = new DataOutputStream[taxaNames.length];
        for (int i = 0; i < taxaNames.length; i++) {
            String outfile = new File(tagBySampleDirS, taxaNames[i]+".tp").getAbsolutePath();
            DataOutputStream dos = IOUtils.getBinaryWriter(outfile);
            try {
                dos.writeInt(this.tagLengthInLong);
                dos.writeInt(-1);
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            taxaWriterMap.put(taxaNames[i], dos);
            doss[i] = dos;
        }          
        try {
            BufferedReader br1 = null;
            BufferedReader br2 = null;
            if (fastqR1.endsWith(".gz")) {
                br1 = IOUtils.getTextGzipReader(fastqR1);
            }
            else {
                br1 = IOUtils.getTextReader(fastqR1);
            }
            if (fastqR2.endsWith(".gz")) {
                br2 = IOUtils.getTextGzipReader(fastqR2);
            }
            else {
                br2 = IOUtils.getTextReader(fastqR2);
            }
            Set<String> bSetR1 = barcodeR1TaxaMap.keySet();
            String[] barcodeR1 = bSetR1.toArray(new String[bSetR1.size()]);
            Set<String> bSetR2 = barcodeR2TaxaMap.keySet();
            String[] barcodeR2 = bSetR2.toArray(new String[bSetR2.size()]);
            Arrays.sort(barcodeR1);
            Arrays.sort(barcodeR2);
            String temp1 = null;
            String temp2 = null;
            int index1 = -1;
            int index2 = -1;
            Set<String> taxaSR1 = null;
            Set<String> taxaSR2 = null;
            DataOutputStream dos = null;
            int totalCnt = 0;
            int processedCnt = 0;
            System.out.println("Parsing " + fastqR1 + "\t" + fastqR2);
            String readR1 = null;
            String readR2 = null;
            byte readR1Len = 0;
            byte readR2Len = 0;
            HashByteByteMap ascIIByteMap = BaseEncoder.getAscIIBaseCodingMap();
            while ((temp1 = br1.readLine()) != null) {
                temp2 = br2.readLine();
                totalCnt++;
                if (totalCnt%10000000 == 0) {
                    System.out.println("Total read count: "+String.valueOf(totalCnt)+"\tPassed read count: "+processedCnt);
                }
                temp1 = br1.readLine(); temp2 = br2.readLine();
                br1.readLine(); br2.readLine();
                br1.readLine(); br2.readLine();
                index1 = Arrays.binarySearch(barcodeR1, temp1);
                index2 = Arrays.binarySearch(barcodeR2, temp2);
                if (index1 == -1 || index2 == -1) {
                    continue;
                }
                index1 = -index1 - 2;
                index2 = -index2 - 2;
                taxaSR1 = barcodeR1TaxaMap.get(barcodeR1[index1]);
                taxaSR2 = barcodeR2TaxaMap.get(barcodeR2[index2]);
                Set<String> newSet = new HashSet<>(taxaSR1);
                newSet.retainAll(taxaSR2);
                if (newSet.size() != 1) {
                    continue;
                }      
                readR1 = this.getChimericRemovedRead(cutter1, cutter2, temp1, barcodeR1[index1].length());
                readR2 = this.getChimericRemovedRead(cutter1, cutter2, temp2, barcodeR2[index2].length());
                if (readR1.length()>this.setReadLength) readR1 = readR1.substring(0, this.setReadLength);
                if (readR2.length()>this.setReadLength) readR2 = readR2.substring(0, this.setReadLength);
                if (readR1.contains("N") || readR2.contains("N")) {
                    continue;
                }
                readR1Len = (byte)readR1.length();
                if (readR1Len < this.minReadLength) continue;
                readR2Len = (byte)readR2.length();
                if (readR2Len < this.minReadLength) continue;
                long[] tag = TagUtils.getTagFromReads(readR1, readR2, ascIIByteMap, tagLengthInLong);
                dos = taxaWriterMap.get(newSet.toArray(new String[newSet.size()])[0]);
                for (int i = 0; i < tag.length; i++) {
                    dos.writeLong(tag[i]);
                }
                dos.writeByte(readR1Len);
                dos.writeByte(readR2Len);
                dos.writeInt(1);               
                processedCnt++;
            }
            for (int i = 0; i < doss.length; i++) {
                doss[i].flush();
                doss[i].close();
            }
            br1.close();
            br2.close();
            System.out.println("Finished parsing " + fastqR1 + "\t" + fastqR2);
            System.out.println("Total read count: "+String.valueOf(totalCnt)+" \tPassed read count: "+processedCnt);
        }
        catch (Exception e) {
            e.printStackTrace();
            System.out.println("Something wrong while reading ");
            System.out.println(fastqR1);
            System.out.println(fastqR2);
            System.exit(1);
            System.out.println("Program quits");
        }
    }
    
    private String getChimericRemovedRead (String cutter1, String cutter2, String read, int barcodeLength) {
        read = read.substring(barcodeLength, read.length());
        int index1 = read.indexOf(cutter1);
        int index2 = read.indexOf(cutter2);
        if (index1 < 0) {
            if (index2 < 0) {
                return read;
            }
            else {
                return read.substring(0, index2);
            }
        }
        else {
            if (index2 < 0) {
                return read.substring(0, index1);
            }
            else {
                if (index1 < index2) return read.substring(0, index1);
                else return read.substring(0, index2);
            }
        }
    }
    
    public void compressTagsBySample (String tagBySampleDirS) {
        File[] fs = new File(tagBySampleDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".tp");
        Arrays.sort(fs);
        int[][] indices = PArrayUtils.getSubsetsIndicesBySubsetSize(fs.length, this.numThreads);
        for (int i = 0; i < indices.length; i++) {
            List<File> subFList = new ArrayList();
            for (int j = indices[i][0]; j < indices[i][1]; j++) {
                subFList.add(fs[j]);
            }
            subFList.parallelStream().forEach(f -> {
                String taxonName = f.getName().replaceFirst(".tp", "");
                String oufileS = new File (tagBySampleDirS, taxonName+".tas").getAbsolutePath();
                TagAnnotations ta = new TagAnnotations (f.getAbsolutePath());
                ta.collapseCounts(1);
                ta.writeBinaryFile(oufileS);
                //ta.writeTextFile(oufileS);
                f.delete();
            });
        }
    }
}
