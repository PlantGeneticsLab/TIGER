package pgl.app.cpScore;

import com.koloboke.collect.map.hash.HashByteByteMap;
import com.koloboke.collect.set.hash.HashIntSet;
import com.koloboke.collect.set.hash.HashIntSets;
import com.koloboke.collect.set.hash.HashLongSet;
import com.koloboke.collect.set.hash.HashLongSets;
import java.io.DataOutputStream;
import pgl.infra.dna.BaseEncoder;
import pgl.infra.dna.FastaByte;

import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOUtils;

public class ReferenceKmerLib {
    int kmerLength = -1;
    HashIntSet[] intSets = null;
    HashLongSet[] longSets = null;
    int barcodeLength = 3;
    int setSize = (int)(Math.pow(4, barcodeLength));
    
    public ReferenceKmerLib(int kmerLength, String inputGenomeFileS) {
        this.kmerLength = kmerLength;
        this.createKmerSet(inputGenomeFileS);
    }

    void createKmerSet (String referenceGenomeFileS) {
//        referenceGenomeFileS = "/Users/feilu/Documents/analysisL/pipelineTest/cpScore/maize_chr12.fa";
//        inputGenomeFileS = "/Users/feilu/Documents/database/maize/reference/AGPv4/maizeAGPv4.fa";
//        inputGenomeFileS = "/Users/feilu/Documents/database/maize/reference/download/Zea_mays.AGPv4.dna.chromosome.10.fa.gz";
        FastaByte f = new FastaByte(referenceGenomeFileS);
        HashByteByteMap ascIIByteMap = BaseEncoder.getAscIIBaseCodingMap();
        System.out.println("Building kmer library from reference...");
        System.out.println("KmerLength = "+String.valueOf(kmerLength)+ " bp");
        long start = System.nanoTime();
        if (kmerLength == 16) {
            this.intSets = this.getIntKmerSet(f, ascIIByteMap);
        }
        else if (kmerLength == 32) {
            this.longSets = this.getLongKmerSets(f, ascIIByteMap);
        }
        System.out.println(Benchmark.getTimeSpanSeconds(start)+" seconds used to build Kmer library");
    }
    
    void writeBinaryFile (String libFileS) {
//        libFileS = "/Users/feilu/Documents/analysisL/pipelineTest/cpScore/kmerLib/kmerLib.bin";
        try {
            DataOutputStream dos  = IOUtils.getBinaryWriter(libFileS);
            dos.writeInt(kmerLength);
            dos.writeInt(barcodeLength);
            if (kmerLength == 16) {
                for (int i = 0; i < intSets.length; i++) {
                    int[] kmerD = intSets[i].toIntArray();
                    dos.writeInt(kmerD.length);
                    for (int j = 0; j < kmerD.length; j++) {
                        dos.writeInt(kmerD[j]);
                    }
                }
            }
            else if (kmerLength == 32) {
                for (int i = 0; i < longSets.length; i++) {
                    long[] kmerD = longSets[i].toLongArray();
                    dos.writeInt(kmerD.length);
                    for (int j = 0; j < kmerD.length; j++) {
                        dos.writeLong(kmerD[j]);
                    }
                }               
            }
            dos.flush();
            dos.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("Kmer library is written to " + libFileS);
    }
    
    public static int getBarcodeIndex (byte[] seqByte, int startIndex, int barcodeLength) {
        int index = 0;
        for (int i = 0; i < barcodeLength; i++) {
            int base = 0;
            if (i == 0) {
                base = seqByte[i+startIndex];
            }
            else {
                base = (int)(Math.pow(4, i))*seqByte[i+startIndex];
            }
            index+=base;
        }
        return index;
    }
    
    private HashLongSet[] getLongKmerSets (FastaByte f, HashByteByteMap ascIIByteMap) {
        int genomeSize = (int)f.getTotalSeqLength();
        HashLongSet[] kmerSets = new HashLongSet[setSize];
        for (int i = 0; i < setSize; i++) {
            kmerSets[i] = HashLongSets.newMutableSet(genomeSize/setSize);
        }
        int barcodeIndex = 0;
        for (int k = 0; k < f.getSeqNumber(); k++) {
            String seq = f.getSeq(k);
            byte[] bArray = seq.getBytes();
            for (int i = 0; i < bArray.length; i++) {
                bArray[i] = ascIIByteMap.get(bArray[i]);
            }
            int mark = 0;
            boolean flag = false;
            for (int i = 0; i < bArray.length-kmerLength+1; i++) {
                flag = false;
                for (int j = mark; j < i+kmerLength; j++) {
                    if (bArray[j] >3) {
                        i = j;
                        flag = true;
                        break;
                    }
                }
                if (flag) {
                    mark = i + 1;
                    continue;
                }
                else {
                    mark = i + kmerLength;
                }
                barcodeIndex = this.getBarcodeIndex(bArray, i, barcodeLength);
                long kmerL = BaseEncoder.getLongFromSubBaseCodingArray(bArray, i, i + kmerLength);
                if (!kmerSets[barcodeIndex].contains(kmerL)) kmerSets[barcodeIndex].add(kmerL);
                int pos = i+1;
                if (pos%50000000 == 0) {
                    long total = 0;
                    for (int u = 0; u < kmerSets.length; u++) {
                        total += (long)kmerSets[u].size();
                    }
                    System.out.println("Chromosome: "+f.getDescription(k)+". Length = "+String.valueOf(bArray.length)+"bp. Position: "+String.valueOf(pos) + ". Kmer set size: " + String.valueOf(total));
                }
            }
        }
        return kmerSets;
    }
    
    private HashIntSet[] getIntKmerSet (FastaByte f, HashByteByteMap ascIIByteMap) {
        int genomeSize = (int)f.getTotalSeqLength();
        HashIntSet[] kmerSets = new HashIntSet[setSize];
        for (int i = 0; i < setSize; i++) {
            kmerSets[i] = HashIntSets.newMutableSet(genomeSize/setSize);
        }
        int barcodeIndex = 0;
        for (int k = 0; k < f.getSeqNumber(); k++) {
            String seq = f.getSeq(k);
            byte[] bArray = seq.getBytes();
            for (int i = 0; i < bArray.length; i++) {
                bArray[i] = ascIIByteMap.get(bArray[i]);
            }
            int mark = 0;
            boolean flag = false;
            for (int i = 0; i < bArray.length-kmerLength+1; i++) {
                flag = false;
                for (int j = mark; j < i+kmerLength; j++) {
                    if (bArray[j] >3) {
                        i = j;
                        flag = true;
                        break;
                    }
                }
                if (flag) {
                    mark = i + 1;
                    continue;
                }
                else {
                    mark = i + kmerLength;
                }
                barcodeIndex = this.getBarcodeIndex(bArray, i, barcodeLength);
                int kmerL = BaseEncoder.getIntSeqFromSubBaseCodingArray(bArray, i, i + kmerLength);
                if (!kmerSets[barcodeIndex].contains(kmerL)) kmerSets[barcodeIndex].add(kmerL);
                int pos = i+1;
                if (pos%50000000 == 0) {
                    long total = 0;
                    for (int u = 0; u < kmerSets.length; u++) {
                        total += (long)kmerSets[u].size();
                    }
                    System.out.println("Chromosome: "+f.getDescription(k)+". Length = "+String.valueOf(bArray.length)+"bp. Position: "+String.valueOf(pos) + ". Kmer set size: " + String.valueOf(total));
                }
            }
        }
        return kmerSets;
    }
}

