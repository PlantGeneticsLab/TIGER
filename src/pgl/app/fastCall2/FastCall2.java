package pgl.app.fastCall2;

import com.koloboke.collect.map.hash.HashByteByteMap;
import com.koloboke.collect.map.hash.HashByteByteMaps;
import pgl.AppUtils;
import pgl.infra.dna.allele.AlleleEncoder;
import pgl.infra.utils.*;
import java.util.*;

public class FastCall2 {

    String[] toolNames = {"disc","blib", "scan"};

    String currentTool = null;

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
            if (args[i].equals("-module")) {
                this.currentTool = args[i+1];
                break;
            }
        }
        if (currentTool.equals("disc")) {
            System.out.println("Discovering genetic variation...");
            new DiscoverVariation(args);
        }
        else if (currentTool.equals("blib")) {
            System.out.println("Building the library of genetic variation from all samples...");
            new BuildVariationLibrary(args);
        }
        else if (currentTool.equals("vlib")) {
            System.out.println("View the library of genetic variation in text");
            new ViewVariationLibrary(args);
        }
        else if (currentTool.equals("clib")) {
            System.out.println("customize the library of genetic variation from position list");
            new CustomizeVariationLibrary(args);
        }
        else if (currentTool.equals("scan")) {
            System.out.println("Genotyping samples based on the variation library...");
            new ScanGenotype(args);
        }
        else {
            System.out.println("Input errors in setting steps of FastCall 2. Programs stops.");
            System.exit(0);
        }
        StringBuilder sb = new StringBuilder("FastCall 2 is finished in ");
        sb.append((float)Benchmark.getTimeSpanHours(timeStart)).append(" hours.");
        System.out.println(sb.toString());
        System.out.println();
    }

    static Dyad<int[][], int[]> getBins (int regionStart, int regionEnd, int binSize) {
        int actualChrLength = regionEnd - regionStart;
        //starting from actual genome position
        int[][] binBound = PArrayUtils.getSubsetsIndicesBySubsetSize (actualChrLength, binSize);
        int[] binStarts = new int[binBound.length];
        for (int i = 0; i < binBound.length; i++) {
            binBound[i][0] = binBound[i][0]+regionStart;
            binBound[i][1] = binBound[i][1]+regionStart;
            binStarts[i] = binBound[i][0];
        }
        return new Dyad<>(binBound, binStarts);
    }

    static void removeFirstPositionSign (StringBuilder sb) {
        char charToRemove = '^';
        for (int i = 0; i < sb.length(); i++) {
            if (sb.charAt(i) == charToRemove) {
                sb.deleteCharAt(i);
                sb.deleteCharAt(i);
                i--;
            }
        }
    }
}
