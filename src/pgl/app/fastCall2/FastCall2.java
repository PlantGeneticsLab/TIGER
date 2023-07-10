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
    //genome block size for variation discovery. The max bin size should be less than 2^23. see {@link AllelePackage}
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

}
