package pgl.app.fastCall3;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import pgl.AppAbstract;
import pgl.PGLAPPEntrance;
import pgl.PGLConstraints;
import pgl.infra.dna.FastaBit;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.Dyad;
import pgl.infra.utils.IOUtils;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.*;
import java.util.concurrent.atomic.LongAdder;

/**
 * A class for building a variation library from genomic data in the FastCall3 pipeline.
 * <p>
 * This class extends {@link AppAbstract} to provide functionality for processing genomic data,
 * identifying variations, and building a comprehensive variation library. It supports multi-threaded
 * processing and can handle specific genomic regions for targeted analysis.
 *
 * <p>Key features include:
 * <ul>
 *   <li>Processing of reference genome and individual genotype data</li>
 *   <li>Support for parallel processing of multiple taxa</li>
 *   <li>Filtering of variations based on minor allele occurrence threshold</li>
 *   <li>Generation of a variation library for downstream analysis</li>
 * </ul>
 *
 * @author Fei Lu
 * @version 3.0
 * @since 1.0
 */
class BuildVariationLibraryF3 extends AppAbstract {

    //Reference genome file with an index file (.fai). The reference should be in Fasta format. Chromosomes are labled as 1-based numbers (1,2,3,4,5...).
    String referenceFileS = null;
    //Current chromosome for variation calling
    short chrom = Short.MIN_VALUE;
    //Starting position of the specified region for variation calling, inclusive
    int regionStart = Integer.MIN_VALUE;
    //Ending position the specified regionfor variation calling, exclusive
    int regionEnd = Integer.MIN_VALUE;
    //Minor allele occurance threshhold, representing the minimum number of taxa where the minor allele exist. It is 2 by default.
    int maoThresh = 2;
    //Number of threads (taxa number to be processed at the same time)
    int threadsNum = PGLConstraints.parallelLevel;
    //Individual genotype directory
    String ingDirS = null;
    //Variation library directory
    String vLibDirS = null;

    String[] taxaNames = null;

    public BuildVariationLibraryF3(String[] args) {
        this.creatAppOptions();
        this.retrieveAppParameters(args);
        this.mkLibrary();
    }


    @Override
    public void creatAppOptions() {
        options.addOption("app", true, "App name.");
        options.addOption("mod", true, "Module name of FastCall 3.");
        options.addOption("a", true, "Reference genome file with an index file (.fai). The reference should be in Fasta format. " +
            "Chromosomes are labelled as numbers (1,2,3,4,5...).");
        options.addOption("b", true, "Chromosome or region on which genotyping will be performed (e.g. chromosome 1 " +
            "is designated as 1. Region 1bp to 100000bp on chromosome 1 is 1:1,100000).");
        options.addOption("c", true, "Minor allele occurrence threshold, representing the minimum number of " +
            "taxa where the minor allele exist. It is 2 by default.");
        options.addOption("d", true, "Number of threads (taxa number to be processed at the same time). It is 32 by default.");
        options.addOption("e", true, "Individual genotype directory.");
        options.addOption("f", true, "Variation library directory.");
    }


    @Override
    public void retrieveAppParameters(String[] args) {
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine line = parser.parse(options, args);
            String inOpt = null;
            this.referenceFileS = line.getOptionValue("a");
            String[] tem = line.getOptionValue("b").split(":");
            this.chrom = Short.parseShort(tem[0]);
            long start = System.nanoTime();
            System.out.println("Reading reference genome from "+ referenceFileS);
            FastaBit genomeFa = new FastaBit(referenceFileS);
            System.out.println("Reading reference genome took " + String.format("%.2f", Benchmark.getTimeSpanSeconds(start)) + "s");
            int chromIndex = genomeFa.getIndexByDescription(String.valueOf(this.chrom));
            if (tem.length == 1) {
                this.regionStart = 1;
                this.regionEnd = genomeFa.getSeqLength(chromIndex)+1;
            }
            else if (tem.length == 2) {
                tem = tem[1].split(",");
                this.regionStart = Integer.parseInt(tem[0]);
                this.regionEnd = Integer.parseInt(tem[1])+1;
            }
            inOpt = line.getOptionValue("c");
            if (inOpt != null) {
                this.maoThresh = Integer.parseInt(inOpt);
                inOpt = null;
            }
            inOpt = line.getOptionValue("d");
            if (inOpt != null) {
                this.threadsNum = Integer.parseInt(inOpt);
                inOpt = null;
            }
            this.ingDirS = line.getOptionValue("e");
            this.vLibDirS = line.getOptionValue("f");
        }
        catch(Exception e) {
            e.printStackTrace();
            System.out.println("\nThere are input errors in the command line. Program stops.");
            this.printInstructionAndUsage();
            System.exit(0);
        }
    }

    @Override
    public void printInstructionAndUsage() {
        System.out.println(PGLAPPEntrance.getTIGERIntroduction());
        System.out.println("Below are the commands of module \"blib\" in FastCall 3.");
        this.printUsage();
    }

/* <<<<<<<<<<<<<<  ✨ Windsurf Command ⭐ >>>>>>>>>>>>>>>> */
    /**
     * Generate a genetic variation library.
     *
     * This function is used to generate a genetic variation library from
     * individual genotypes. It takes the input directory of individual genotypes,
     * the output directory for the library, and the number of threads as parameters.
     * The library is generated in a threaded manner, with each thread building a
     * subset of the library. The library is written to a binary file in the output
     * directory.
     *
     * @param ingDirS the input directory of individual genotypes
     * @param vLibDirS the output directory for the library
     * @param threadsNum the number of threads used to build the library
     */
/* <<<<<<<<<<  109de986-8334-4f96-8bc9-ed2e33d18cf4  >>>>>>>>>>> */
    private void mkLibrary () {
        List<File> ingTaxaDirList = IOUtils.getDirListInDir(this.ingDirS);
        taxaNames = new String[ingTaxaDirList.size()];
        for (int i = 0; i < ingTaxaDirList.size(); i++) {
            taxaNames[i] = ingTaxaDirList.get(i).getName();
        }
        Arrays.sort(taxaNames);
        Dyad<int[][], int[]> d = FastCall3.getBins(this.regionStart, this.regionEnd, FastCall3.disBinSize);
        int[][] binBound = d.getFirstElement();
        int[] binStarts = d.getSecondElement();
        StringBuilder sb = new StringBuilder();
        List<Future<IndividualGenotypeF3>> futureList = new ArrayList<>();
        List<IndividualGenotypeF3> ingList = new ArrayList<>();
        List<VariationLibraryF3> vlList = new ArrayList<>();
        for (int i = 0; i < binBound.length; i++) {
            futureList.clear();
            ingList.clear();
            try {
                LongAdder counter = new LongAdder();
                ExecutorService pool = Executors.newFixedThreadPool(this.threadsNum);
                sb.setLength(0);
                sb.append(chrom).append("_").append(binBound[i][0]).append("_").append(binBound[i][1]).append(".ing.gz");
                for (int j = 0; j < taxaNames.length; j++) {
                    String fileS = new File (ingTaxaDirList.get(j), sb.toString()).getAbsolutePath();
                    TaxonRead tr = new TaxonRead(fileS);
                    Future<IndividualGenotypeF3> result = pool.submit(tr);
                    futureList.add(result);
                }
                pool.shutdown();
                pool.awaitTermination(Long.MAX_VALUE, TimeUnit.SECONDS);

                for (int j = 0; j < futureList.size(); j++) {
                    IndividualGenotypeF3 ing = futureList.get(j).get();
                    if (ing == null) continue;
                    ingList.add(ing);
                }
                Collections.sort(ingList);
            }
            catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
            VariationLibraryF3 vl = new VariationLibraryF3(ingList, maoThresh, FastCall3.maxAltNum, chrom, binStarts[i]);
            vlList.add(vl);
        }
        VariationLibraryF3 chromVl = VariationLibraryF3.getInstance(vlList);
        File f = new File (vLibDirS);
        f.mkdir();
        sb.setLength(0);
        sb.append(chrom).append("_").append(this.regionStart).append("_").append(regionEnd).append(".lib.gz");
        chromVl.writeBinaryFileS(new File (f, sb.toString()).getAbsolutePath());
        System.out.println("Building genetic variation library is finished.");
    }

    /**
     * A Callable class for reading individual genotype data from a file.
     * <p>
     * This class implements {@link Callable} to enable parallel reading of individual
     * genotype data files. It reads and parses genotype information from a specified file
     * and returns an {@link IndividualGenotypeF3} object containing the processed data.
     * 
     * <p>Usage example:
     * <pre>
     * ExecutorService executor = Executors.newFixedThreadPool(threadsNum);
     * Future<IndividualGenotypeF3> future = executor.submit(new TaxonRead(genotypeFile));
     * IndividualGenotypeF3 genotype = future.get();
     * </pre>
     *
     * @see IndividualGenotypeF3
     * @see java.util.concurrent.Callable
     */
    class TaxonRead implements Callable<IndividualGenotypeF3> {
        /** The file path to read the individual genotype data from */
        String fileS = null;

        /**
         * Constructs a new TaxonRead instance for the specified genotype file.
         *
         * @param fileS the path to the file containing individual genotype data
         */
        public TaxonRead (String fileS) {
            this.fileS = fileS;
        }

        /**
         * Read individual genotype from a file.
         * @return the read individual genotype
         * @throws Exception if any error occurs
         */
        @Override
        public IndividualGenotypeF3 call() throws Exception {
            File f = new File (fileS);
            if (!f.exists()) return null;
            IndividualGenotypeF3 ing = new IndividualGenotypeF3(this.fileS);
            return ing;
        }
    }
}
