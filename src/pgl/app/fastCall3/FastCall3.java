package pgl.app.fastCall3;

import com.koloboke.collect.map.hash.HashByteByteMap;
import com.koloboke.collect.map.hash.HashByteByteMaps;
import pgl.infra.dna.allele.AlleleEncoder;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.Dyad;
import pgl.infra.utils.PArrayUtils;

/**
 * FastCall3: A high-performance pipeline for genetic variant discovery and genotyping.
 * <p>
 * This class serves as the main entry point for the FastCall3 pipeline, which provides
 * a comprehensive solution for processing next-generation sequencing data to identify
 * and analyze genetic variations. The pipeline is designed for high-throughput
 * processing of large-scale genomic datasets with a focus on accuracy and efficiency.
 *
 * <p>Main components of the pipeline include:
 * <ul>
 *   <li><b>disc</b>: Discovers genetic variations from aligned sequencing data</li>
 *   <li><b>blib</b>: Builds a comprehensive variation library from multiple samples</li>
 *   <li><b>vlib</b>: Views and inspects the variation library in text format</li>
 *   <li><b>clib</b>: Customizes the variation library based on specific positions</li>
 *   <li><b>scan</b>: Performs genotyping of samples using the variation library</li>
 * </ul>
 *
 * <p>Key features:
 * <ul>
 *   <li>Efficient processing of large genomic datasets</li>
 *   <li>Support for both SNP and indel variant calling</li>
 *   <li>Parallel processing capabilities for improved performance</li>
 *   <li>Configurable parameters for variant calling quality control</li>
 *   <li>Modular design for flexible pipeline configuration</li>
 * </ul>
 *
 * @author Fei Lu
 * @version 3.0
 * @since 1.0
 * @see DiscoverVariationF3
 * @see BuildVariationLibraryF3
 * @see ViewVariationLibraryF3
 * @see CustomizeVariationLibraryF3
 * @see ScanGenotypeF3
 */
public class FastCall3 {

    String[] toolNames = {"disc","blib", "scan"};

    String currentTool = null;

    //genome block size for variation discovery. The max bin size should be less than 2^23. see {@link AllelePackage}
    static int disBinSize = 5000000;
    //genome block size for genotype scanning
    static int scanBinSize = 5000000;
    // maximum of the number of alternative alleles
    static int maxAltNum = 2;
    static int maxIndelLength = 63;
    //A, C, G, T, -, +
    static final byte[] pileupAlleleAscIIs = {65, 67, 71, 84, 45, 43};

    static final HashByteByteMap pileupAscIIToAlleleCodingMap =
            HashByteByteMaps.getDefaultFactory().withDefaultValue((byte)-1).newImmutableMap(pileupAlleleAscIIs, AlleleEncoder.alleleCodings);

    public FastCall3(String[] args) {
        long timeStart = System.nanoTime();
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("-mod")) {
                this.currentTool = args[i+1];
                break;
            }
        }
        if (currentTool.equals("disc")) {
            System.out.println("Discovering genetic variation...");
            new DiscoverVariationF3(args);
        }
        else if (currentTool.equals("blib")) {
            System.out.println("Building the library of genetic variation from all samples...");
            new BuildVariationLibraryF3(args);
        }
        else if (currentTool.equals("vlib")) {
            System.out.println("View the library of genetic variation in text");
            new ViewVariationLibraryF3(args);
        }
        else if (currentTool.equals("clib")) {
            System.out.println("customize the library of genetic variation from position list");
            new CustomizeVariationLibraryF3(args);
        }
        else if (currentTool.equals("scan")) {
            System.out.println("Genotyping samples based on the variation library...");
            new ScanGenotypeF3(args);
        }
        else {
            System.out.println("Input errors in setting steps of FastCall 3. Programs stops.");
            System.exit(0);
        }
        StringBuilder sb = new StringBuilder("FastCall 3 is finished in ");
        sb.append((float)Benchmark.getTimeSpanHours(timeStart)).append(" hours.");
        System.out.println(sb.toString());
        System.out.println();
    }

    /**
     * Divides a given genomic region into bins of a specified size.
     * @param regionStart the starting position of the genomic region
     * @param regionEnd the ending position of the genomic region
     * @param binSize the size of each bin
     * @return a Dyad containing two arrays: the first array contains the
     *         bounds of each bin, and the second array contains the
     *         starting positions of each bin
     */
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

    /**
     * Given a StringBuilder, remove the first position sign "^" and the character right after it.
     * This is used to remove the first position sign from a pileup string.
     * @param sb the StringBuilder to operate on
     */
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
