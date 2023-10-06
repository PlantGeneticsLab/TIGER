/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.app.cpScore;

import com.koloboke.collect.map.hash.HashByteByteMap;
import com.koloboke.collect.map.hash.HashIntIntMap;
import com.koloboke.collect.map.hash.HashIntIntMaps;
import com.koloboke.collect.map.hash.HashLongIntMap;
import com.koloboke.collect.map.hash.HashLongIntMaps;
import gnu.trove.list.array.TIntArrayList;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import pgl.infra.dna.BaseEncoder;
import pgl.infra.dna.FastaByte;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.PArrayUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.IOUtils;

/**
 *
 * @author feilu
 */
public class GenomeProfiler {
    int kmerLength = -1;
    int barcodeLength = -1;
    HashLongIntMap[] longMaps = null; 
    HashIntIntMap[] intMaps = null;

    
    public GenomeProfiler (String libFileS, String referenceGenomeFileS, String anotherGenomeFileS, String outputDirS) {     
        this.initializeKmerCountMap(libFileS);
        this.countKmers(anotherGenomeFileS);
        this.writeCpScore(referenceGenomeFileS, anotherGenomeFileS, outputDirS);
    }
    
    void writeCpScore (String referenceGenomeFileS, String anotherGenomeFileS, String outputDirS) {
//        referenceGenomeFileS = "/Users/feilu/Documents/analysisL/pipelineTest/cpScore/maize_chr12.fa";
//        outputDirS = "/Users/feilu/Documents/analysisL/pipelineTest/cpScore/result";
        File outputDir = new File(outputDirS);
        outputDir.mkdir();
        int fragmentSize = 100000;
        HashByteByteMap ascIIByteMap = BaseEncoder.getAscIIBaseCodingMap();
        FastaByte f = new FastaByte(referenceGenomeFileS);
        HashMap<Integer, String> chrSeqMap = new HashMap();
        for (int i = 0; i < f.getSeqNumber(); i++) {
            chrSeqMap.put(Integer.valueOf(f.getDescription(i)), f.getSeq(i));
        }
        Set<Map.Entry<Integer, String>> chrSeqset = chrSeqMap.entrySet();
        System.out.println("Writing CpScore by chromosomes...");
        long start = System.nanoTime();
        chrSeqset.parallelStream().forEach(entry -> {
            int chr = entry.getKey();
            String seq = entry.getValue();
            byte[] bArray = seq.getBytes();
            for (int i = 0; i < bArray.length; i++) {
                bArray[i] = ascIIByteMap.get(bArray[i]);
            }
            StringBuilder sb  = new StringBuilder();
            sb.append("chr").append(PStringUtils.getNDigitNumber(3, chr)).append("_CpScore.txt.gz");
            String outfileS = new File(outputDir, sb.toString()).getAbsolutePath();
            int[][] bound = PArrayUtils.getSubsetsIndicesBySubsetSize(seq.length(), fragmentSize);
            try {
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                bw.write("CpScore");
                bw.newLine();
                for (int i = 0; i < bound.length; i++) {
                    int intervalSize = bound[i][1] - bound[i][0];
                    TIntArrayList[] kmerCountList = new TIntArrayList[intervalSize];
                    for (int j = 0; j < kmerCountList.length; j++) kmerCountList[j] = new TIntArrayList();
                    int startIndex = bound[i][0]-kmerLength+1;
                    if (startIndex < 0) startIndex = 0;
                    int endIndex = bound[i][1];
                    if (endIndex-1+kmerLength > seq.length()) endIndex = seq.length()-kmerLength+1;
                    int mark = startIndex;
                    boolean flag = false;
                    for (int j = startIndex; j < endIndex; j++) {
                        flag = false;
                        for (int k = mark; k < j+kmerLength; k++) {
                            if (bArray[k] >3) {
                                j = k;
                                flag = true;
                                break;
                            }
                        }
                        if (flag) {
                            mark = j + 1;
                            continue;
                        }
                        else {
                            mark = j + kmerLength;
                        }
                        int count = 0;
                        
                        if (intMaps != null) {
                            int barcodeIndex = ReferenceKmerLib.getBarcodeIndex(bArray, j, barcodeLength);
                            int query = BaseEncoder.getIntSeqFromSubBaseCodingArray(bArray, j, j + kmerLength);
                            count = intMaps[barcodeIndex].get(query);
                        }
                        else if (longMaps != null) {
                            int barcodeIndex = ReferenceKmerLib.getBarcodeIndex(bArray, j, barcodeLength);
                            long query = BaseEncoder.getLongFromSubBaseCodingArray(bArray, j, j + kmerLength);
                            count = longMaps[barcodeIndex].get(query);
                        }
                        int offSet = bound[i][0]-j;
                        if (offSet > 0) {
                            for (int k = j; k < j+kmerLength-offSet; k++) {
                                int a = offSet+k-bound[i][0];
                                kmerCountList[offSet+k-bound[i][0]].add(count);
                            }
                        }
                        else {
                            int end = j+kmerLength;
                            if (end > bound[i][1]) end = bound[i][1];
                            for (int k = j; k < end; k++) {
                                kmerCountList[k-bound[i][0]].add(count);
                            }
                        }
                    }
                    for (int j = 0; j < kmerCountList.length; j++) {
                        int[] kmerCount = kmerCountList[j].toArray();
                        sb = new StringBuilder();
                        if (kmerCount.length == 0) {
                            sb.append("NA");
                        }
                        else {
                            float aver = 0;
                            for (int k = 0; k < kmerCount.length; k++) {
                                aver+=(float)((double)kmerCount[k]/kmerCount.length);
                            }
                            sb.append(aver);
                        }
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                    if ((bound[i][1])%50000000 == 0) System.out.println(String.valueOf(bound[i][1])+" sites on chr "+String.valueOf(chr)+" are finished.");
                }
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        StringBuilder sb = new StringBuilder();
        sb.append("Writing copy number score from ").append(anotherGenomeFileS).append(" is finished. Time span: ").append(Benchmark.getTimeSpanSeconds(start)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
        System.out.println(sb.toString());
    }
          
    void countKmers (String inputGenomeFileS) {
//        inputGenomeFileS = "/Users/feilu/Documents/analysisL/pipelineTest/cpScore/maize_chr12.fa";
        FastaByte f = new FastaByte (inputGenomeFileS);
        
        HashByteByteMap ascIIByteMap = BaseEncoder.getAscIIBaseCodingMap();
        System.out.println("Start counting kmers in genome " + inputGenomeFileS);
        long start = System.nanoTime();
        if (intMaps != null) {
            for (int i = 0; i < f.getSeqNumber(); i++) {
                byte[] bArray = f.getSeq(i).getBytes();
                for (int j = 0; j < bArray.length; j++) {
                    bArray[j] = ascIIByteMap.get(bArray[j]);
                }
                int mark = 0;
                boolean flag = false;
                for (int j = 0; j < bArray.length-kmerLength+1; j++) {
                    flag = false;
                    for (int k = mark; k < j+kmerLength; k++) {
                        if (bArray[k] >3) {
                            j = k;
                            flag = true;
                            break;
                        }
                    }
                    if (flag) {
                        mark = j + 1;
                        continue;
                    }
                    else {
                        mark = j + kmerLength;
                    }
                    int barcodeIndex = ReferenceKmerLib.getBarcodeIndex(bArray, j, barcodeLength);
                    int kmerV = BaseEncoder.getIntSeqFromSubBaseCodingArray(bArray, j, j + kmerLength);
                    if (intMaps[barcodeIndex].containsKey(kmerV)) intMaps[barcodeIndex].addValue(kmerV, 1);
                    int rKmerV = BaseEncoder.getIntReverseComplement(kmerV, kmerLength);
                    if (intMaps[barcodeIndex].containsKey(rKmerV)) intMaps[barcodeIndex].addValue(rKmerV, 1);          
                    int pos = j+1;
                    if (pos%50000000 == 0) {
                        System.out.println(inputGenomeFileS + ". Chromosome: "+f.getDescription(i)+". Length = "+String.valueOf(bArray.length)+"bp. Position: "+String.valueOf(pos));
                    }
                }
            }
        }
        else if (longMaps != null) {
            for (int i = 0; i < f.getSeqNumber(); i++) {
                byte[] bArray = f.getSeq(i).getBytes();
                for (int j = 0; j < bArray.length; j++) {
                    bArray[j] = ascIIByteMap.get(bArray[j]);
                }
                int mark = 0;
                boolean flag = false;
                for (int j = 0; j < bArray.length-kmerLength+1; j++) {
                    flag = false;
                    for (int k = mark; k < j+kmerLength; k++) {
                        if (bArray[k] >3) {
                            j = k;
                            flag = true;
                            break;
                        }
                    }
                    if (flag) {
                        mark = j + 1;
                        continue;
                    }
                    else {
                        mark = j + kmerLength;
                    }
                    int barcodeIndex = ReferenceKmerLib.getBarcodeIndex(bArray, j, barcodeLength);
                    long kmerV = BaseEncoder.getLongFromSubBaseCodingArray(bArray, j, j + kmerLength);
                    if (longMaps[barcodeIndex].containsKey(kmerV)) {
                        longMaps[barcodeIndex].addValue(kmerV, 1);
                    }
                    long rKmerV = BaseEncoder.getLongReverseComplement(kmerV, kmerLength);
                    if (longMaps[barcodeIndex].containsKey(rKmerV)) longMaps[barcodeIndex].addValue(rKmerV, 1);
                    int pos = j+1;
                    if (pos%50000000 == 0) {
                        System.out.println(inputGenomeFileS + ". Chromosome: "+f.getDescription(i)+". Length = "+String.valueOf(bArray.length)+"bp. Position: "+String.valueOf(pos));
                    }
                }
            }
        }
        StringBuilder sb = new StringBuilder();
        sb.append("Couting kmer is finished. Time span: ").append(Benchmark.getTimeSpanSeconds(start)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
        System.out.println(sb.toString());
    }
    
    void initializeKmerCountMap (String libFileS) {
//        libFileS = "/Users/feilu/Documents/analysisL/pipelineTest/cpScore/kmerLib/kmerLib.bin"; 
        System.out.println("Reading kmer library from " + libFileS);
        long start = System.nanoTime();
        try {
            DataInputStream dis = IOUtils.getBinaryReader(libFileS);
            kmerLength = dis.readInt();
            barcodeLength = dis.readInt();
            int setSize = (int)(Math.pow(4, barcodeLength));
            if (kmerLength == 32) {
                longMaps = new HashLongIntMap[setSize];
            }
            else if (kmerLength == 16) {
                intMaps = new HashIntIntMap[setSize];
            }
            int mapSize = -1;
            for (int i = 0; i < setSize; i++) {
                mapSize = dis.readInt();
                System.out.println("Map size of index "+ String.valueOf(i)+ ": "+mapSize);
                int[] values = new int[mapSize];
                if (kmerLength == 32) {
                    long[] keys = new long[mapSize];
                    for (int j = 0; j < mapSize; j++) {
                        keys[j] = dis.readLong();
                    }
                    longMaps[i] = HashLongIntMaps.newMutableMap(keys, values);
                }
                else if (kmerLength == 16) {
                    int[] keys = new int[mapSize];
                    for (int j = 0; j < mapSize; j++) {
                        keys[j] = dis.readInt();
                    }
                    intMaps[i] = HashIntIntMaps.newMutableMap(keys, values);
                }
            }
            dis.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println(Benchmark.getTimeSpanSeconds(start)+" seconds used to initialize the KmerCount map");
    }
}
