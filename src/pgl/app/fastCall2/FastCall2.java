package pgl.app.fastCall2;

import com.koloboke.collect.map.hash.HashByteByteMap;
import com.koloboke.collect.map.hash.HashByteByteMaps;
import pgl.AppUtils;
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

    static int getCodedPosAlleleIndelLength(int binStart, int position, byte alleleCoding, int indelLength) {
        int v = (position-binStart) << 9;
        v = v + (alleleCoding << 6);
        if (indelLength > 63) indelLength = 63;
        return (v + indelLength);
    }

    static int getAllelePosition(int codedPosAlleleIndelLength, int binStart) {
        int v = (codedPosAlleleIndelLength >>> 9) + binStart;
        return v;
    }

    static short getCodedAllele (byte alleleCoding, int indelLength) {
        int v = (alleleCoding << 6);
        if (indelLength > 63) indelLength = 63;
        return (short)(v + indelLength);
    }

    static short getCodedAllele (int codedPosAlleleIndelLength) {
        return (short)(codedPosAlleleIndelLength&511);
    }

    static byte getAlleleCoding(int codedPosAlleleIndelLength) {
        int v = (511 & codedPosAlleleIndelLength) >>> 6;
        return (byte)v;
    }

    static char getAlleleBase (int codedPosAlleleIndelLength) {
        return AlleleEncoder.alleleCodingToBaseMap.get(getAlleleCoding(codedPosAlleleIndelLength));
    }

    static byte getIndelLength (int codedPosAlleleIndelLength) {
        return (byte)(63 & codedPosAlleleIndelLength);
    }

    static byte getAlleleCodingFromCodedAllele(short codedAllele) {
        int v = (codedAllele>>>6);
        return (byte)v;
    }

    static char getAlleleBaseFromCodedAllele (short codedAllele) {
        return AlleleEncoder.alleleCodingToBaseMap.get(getAlleleCodingFromCodedAllele(codedAllele));
    }

    static byte getIndelLengthFromCodedAllele (short codedAllele) {
        return (byte)(63 & codedAllele);
    }

}
