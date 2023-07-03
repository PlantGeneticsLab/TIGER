package pgl.app.fastCall2;

import com.koloboke.collect.map.hash.HashByteByteMap;
import com.koloboke.collect.map.hash.HashByteByteMaps;
import pgl.AppUtils;
import pgl.infra.dna.BaseEncoder;
import pgl.infra.dna.allele.AlleleEncoder;
import pgl.infra.utils.*;
import java.util.*;

public class FastCall2 {
    //Current step ID of the pipeline
    int step = Integer.MIN_VALUE;
    //genome block size for variation discovery
    static int disBinSize = 5000000;
    //genome block size for genotype scanning
    static int scanBinSize = 20000000;
    // maximum of the number of alternative alleles
    static int maxAltNum = 2;
    static int maxIndelLength = 63;
    //A, C, G, T, -, +
    static final byte[] pileupAlleleAscIIs = {65, 67, 71, 84, 45, 43};

    static final HashByteByteMap pileupAscIIToAlleleCodingMap =
            HashByteByteMaps.getDefaultFactory().withDefaultValue((byte)-1).newImmutableMap(pileupAlleleAscIIs, AlleleEncoder.alleleCodings);

    public FastCall2 (String[] args) {
        long timeStart = System.nanoTime();
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-step")) {
                this.step = Integer.parseInt(args[i+1]);
                break;
            }
        }
        if (step == 1) {
            System.out.println("Running step 1 of FastCall 2...");
            new DiscoverVariation(args);
        }
        else if (step == 2) {
            System.out.println("Running step 2 of FastCall 2...");
            new BuildVariationLibrary(args);
        }
        else if (step == 3) {
            System.out.println("Running step 3 of FastCall 2...");
            new ScanGenotype(args);
        }
        else {
            System.out.println("Input errors in setting steps of FastCall 2. Programs stops.");
            System.exit(0);
        }
        StringBuilder sb = new StringBuilder("FastCall2 is finished in ");
        sb.append((float)Benchmark.getTimeSpanHours(timeStart)).append(" hours.");
        System.out.println(sb.toString());
    }

    public FastCall2 (String parameterFileS) {
        this.runSteps(parameterFileS);
    }

    private void runSteps(String parameterFileS) {
        long timeStart = System.nanoTime();
        Dyad<List<String>, List<String>> d = AppUtils.getParameterList(parameterFileS);
        List<String> pLineList = d.getFirstElement();
        List<String> sLineList = d.getSecondElement();
        this.step = Integer.parseInt(sLineList.get(0).split("\\s+")[1]);
        if (step == 1) {
            System.out.println("Running step 1...");
            new DiscoverVariation(pLineList);
        }
        else if (step == 2) {
            System.out.println("Running step 2...");
            new BuildVariationLibrary(pLineList);
        }
        else if (step == 3) {
            System.out.println("Running step 3...");
            new ScanGenotype(pLineList);
        }
        StringBuilder sb = new StringBuilder("FastCall2 is finished in ");
        sb.append((float)Benchmark.getTimeSpanHours(timeStart)).append(" hours.");
        System.out.println(sb.toString());
    }

    //Variation discovery is performed at the individual level, so the default bin size is larger
    static Dyad<int[][], int[]> getBinsDiscovery (int regionStart, int regionEnd) {
        int actualChrLength = regionEnd - regionStart;
        //starting from actual genome position
        int[][] binBound = PArrayUtils.getSubsetsIndicesBySubsetSize (actualChrLength, disBinSize);
        int[] binStarts = new int[binBound.length];
        for (int i = 0; i < binBound.length; i++) {
            binBound[i][0] = binBound[i][0]+regionStart;
            binBound[i][1] = binBound[i][1]+regionStart;
            binStarts[i] = binBound[i][0];
        }
        return new Dyad<>(binBound, binStarts);
    }

    //Genotyping/scanning is performed at the population level (>thousands), so the default bin size is smaller
    static Dyad<int[][], int[]> getBinsScanning (int regionStart, int regionEnd) {
        int actualChrLength = regionEnd - regionStart;
        //starting from actual genome position
        int[][] binBound = PArrayUtils.getSubsetsIndicesBySubsetSize (actualChrLength, scanBinSize);
        int[] binStarts = new int[binBound.length];
        for (int i = 0; i < binBound.length; i++) {
            binBound[i][0] = binBound[i][0]+regionStart;
            binBound[i][1] = binBound[i][1]+regionStart;
            binStarts[i] = binBound[i][0];
        }
        return new Dyad<>(binBound, binStarts);
    }

    //the first 23 bits are used to record position in a bin (maximum size = 2^23 - 1 = 8,388,607)
    //the next 3 bits are used to record alleles
    //the last 6 bits are used record the length of indel (max size = 2^6 -1 = 63)
    //Above three components are called a "pack"
    static int getCodedAllelePack(int binStart, int position, byte alleleCoding, int indelLength) {
        int v = (position-binStart) << 9;
        v = v + (alleleCoding << 6);
        if (indelLength > maxIndelLength) indelLength = maxIndelLength;
        return (v + indelLength);
    }

    /**
     * Convert Indel sequence to an array of long
     * @param indelSeq
     * @param indelLength
     * @return
     */
    static long[] getIndelSeqL (String indelSeq, int indelLength) {
        byte[] seqCodings = BaseEncoder.convertToBaseCodingArray(indelSeq.getBytes());
        long[] seqL = null;
        if (seqCodings.length > BaseEncoder.longChunkSize) {
            seqL = new long[2];
            seqL[0] = BaseEncoder.getLongSeqFromSubBaseCodingArray(seqCodings,0, BaseEncoder.longChunkSize);
            seqL[1] = BaseEncoder.getLongSeqFromSubBaseCodingArray(seqCodings, BaseEncoder.longChunkSize, indelLength);
        }
        else {
            seqL = new long[1];
            seqL[0] = BaseEncoder.getLongSeqFromBaseCodingArray(seqCodings);
        }
        return seqL;
    }

    /**
     * Convert to Indel sequence
     * @param seqL
     * @param indelLength
     * @return
     */
    static String getIndelSeq (long[] seqL, int indelLength) {
        return BaseEncoder.getSequenceFromLongs(seqL).substring(0, indelLength);
    }

    static int getAllelePosition(int codedAllelePack, int binStart) {
        int v = (codedAllelePack >>> 9) + binStart;
        return v;
    }

    static short getAlleleCodingLength(byte alleleCoding, int indelLength) {
        int v = (alleleCoding << 6);
        if (indelLength > maxIndelLength) indelLength = maxIndelLength;
        return (short)(v + indelLength);
    }

    static short getAlleleCodingLength(int codedAllelePack) {
        return (short)(codedAllelePack&511);
    }

    static byte getAlleleCoding(int codedAllelePack) {
        int v = (511 & codedAllelePack) >>> 6;
        return (byte)v;
    }

    static char getAlleleBase (int codedAllelePack) {
        return AlleleEncoder.alleleCodingToBaseMap.get(getAlleleCoding(codedAllelePack));
    }

    static byte getIndelLength (int codedAllelePack) {
        return (byte)(maxIndelLength & codedAllelePack);
    }

    static byte getAlleleCodingFromAlleleCodingLength(short alleleCodingLength) {
        int v = (alleleCodingLength>>>6);
        return (byte)v;
    }

    static char getAlleleBaseFromAlleleCodingLength(short alleleCodingLength) {
        return AlleleEncoder.alleleCodingToBaseMap.get(getAlleleCodingFromAlleleCodingLength(alleleCodingLength));
    }

    static byte getIndelLengthFromAlleleCodingLength(short alleleCodingLength) {
        return (byte)(maxIndelLength & alleleCodingLength);
    }

}
