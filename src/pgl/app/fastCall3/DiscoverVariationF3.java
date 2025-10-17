package pgl.app.fastCall3;

import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import pgl.AppAbstract;
import pgl.PGLAPPEntrance;
import pgl.PGLConstraints;
import pgl.infra.dna.BaseEncoder;
import pgl.infra.dna.FastaBit;
import pgl.infra.dna.allele.AlleleEncoder;
import pgl.infra.utils.*;
import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.File;
import java.io.InputStreamReader;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.LongAdder;

/**
 * A comprehensive variant discovery tool for the FastCall3 pipeline.
 * <p>
 * This class extends {@link AppAbstract} to provide functionality for discovering
 * genetic variations from aligned sequencing data. It processes BAM files, applies
 * various quality filters, and identifies SNPs and indels with high confidence.
 *
 * <p>Key features include:
 * <ul>
 *   <li>Support for both whole-genome and targeted region analysis</li>
 *   <li>Configurable quality filters for base quality, mapping quality, and read depth</li>
 *   <li>Parallel processing of multiple samples</li>
 *   <li>Sophisticated variant calling with support for heterozygous and homozygous variants</li>
 *   <li>Handling of both SNPs and indels</li>
 *   <li>Base Alignment Quality (BAQ) computation</li>
 * </ul>
 *
 * <p>The class implements several quality control parameters:
 * <ul>
 *   <li>Mapping Quality (MQ) threshold: {\@code mappingQThresh}</li>
 *   <li>Base Quality (BQ) threshold: {\@code baseQThresh}</li>
 *   <li>Minimum Read Depth Count (MDC): {\@code mdcThresh}</li>
 *   <li>Homozygous Ratio (HoR): {\@code horThresh}</li>
 *   <li>Heterozygous Ratio (HeR): {\@code herThresh}</li>
 *   <li>Third Allele Depth Ratio (TDR): {\@code tdrThresh}</li>
 * </ul>
 *
 * @author Fei Lu
 * @version 3.0
 * @since 1.0
 * @see pgl.AppAbstract
 * @see pgl.infra.dna.FastaBit
 * @see pgl.infra.dna.allele.AlleleEncoder
 */
class DiscoverVariationF3 extends AppAbstract {
    //Reference genome file with an index file (.fai). The reference should be in Fasta format. Chromosomes are labled as numbers (1,2,3,4,5...).
    String referenceFileS = null;
    //The taxaBamMap file contains information of taxon and its corresponding bam files. The bam file should have .bai file in the same folder.
    String taxaBamMapFileS = null;
    //The path of samtools
    String samtoolsPath = null;
    //Individual genotype output directory
    String outputDirS = null;
    //The switch of base alignment quality (BAQ) computation, 0 is disabled and 1 is enabled.
    String baqMode = "-B ";
    //Minimum mapping quality (MQ) for an alignment to be used for variation calling.
    int mappingQThresh = 30;
    //Minimum base quality (BQ) for a base to be used for variation calling.
    int baseQThresh = 20;
    //Minimum read depth count (MDC) for variation calling, meaning that sites with depth lower than the minimum will not be taken into account for variation discovery.
    int mdcThresh = 2;
    //Minimum read depth ratio (MiDR) for variation calling, meaning that sites with depth lower than the MiDR of the individual sequencing coverage will not be considered for variation discovery.
    double mindrThresh = 0.2;
    //Maximum read depth ratio (MaDR) for variation calling, meaning that sites with depth higher than the MaDR of the individual sequencing coverage will not be considered for variation discovery.
    double maxdrTrresh = 3;
    //Homozygous ratio (HoR) for variation calling, meaning that the depth of alternative allele is greater than HoR are considered to homozygous.
    double horThresh = 0.8;
    //Heterozygous ratio (HeR) for variation calling, meaning that the depth of alternative allele is greater than HR and less than (1-HR) are considered to be hets.
    double herThresh = 0.35;
    //Third allele depth ratio (TDR) for variation calling. If the depth of the third allele is greater than TDR by the individual coverage, the site will be ignored. Otherwise, the third allele will be considered as sequencing error.
    double tdrTresh = 0.2;
    //Current chromosome for variation calling
    int chrom = Integer.MIN_VALUE;
    //Starting position of the specified region for variation calling, inclusive
    int regionStart = Integer.MIN_VALUE;
    //Ending position the specified region for variation calling, exclusive
    int regionEnd = Integer.MIN_VALUE;
    //Number of threads (taxa number to be processed at the same time)
    int threadsNum = PGLConstraints.parallelLevel;
    //Two many indels meaning alignment error
    int indelTypeThresh = 1;

    HashMap<String, String[]> taxaBamPathMap = null;
    HashMap<String, Double> taxaCoverageMap = null;
    String[] taxaNames = null;

    public DiscoverVariationF3(String[] args) {
        this.creatAppOptions();
        this.retrieveAppParameters(args);
        this.variationDiscovery();
    }

    @Override
    public void creatAppOptions() {
        options.addOption("app", true, "App name.");
        options.addOption("mod", true, "Module name of FastCall 3.");
        options.addOption("a", true, "Reference genome file with an index file (.fai). The reference should be in Fasta format. " +
                "Chromosomes are labeled as numbers (1,2,3,4,5...). It is recommended to use reference chromosome while perform variation discovery " +
                "for each chromosome because loading the reference genome would be much faster.");
        options.addOption("b", true, "The taxaBamMap file contains information of taxon and its corresponding bam files. " +
                "The bam file should have .bai file in the same folder.");
        options.addOption("c", true, "The switch of base alignment quality (BAQ) computation, 0 is disabled and 1 is enabled. It is 0 by default.");
        options.addOption("d", true, "Minimum mapping quality (MQ) for an alignment to be used for variation calling. It is 30 by default.");
        options.addOption("e", true, "Minimum base quality (BQ) for a base to be used for variation calling. It is 20 by default.");
        options.addOption("f", true, "Minimum read depth count (MDC) for variation calling, meaning that sites with depth lower than " +
                "the minimum will not be taken into account for variation discovery. It is 2 by default.");
        options.addOption("g", true, "Minimum read depth ratio (MiDR) for variation calling, meaning that sites with depth lower than the " +
            "MiDR of the individual sequencing coverage will not be considered for variation discovery. It is 0.2 by default.");
        options.addOption("h", true, "Maximum read depth ratio (MaDR) for variation calling, meaning that sites with depth higher than " +
            "the MaDR of the individual sequencing coverage will not be considered for variation discovery. It is 3 by default.");
        options.addOption("i", true, "Homozygous ratio (HoR) for variation calling, meaning that the depth of alternative allele is " +
            "greater than HoR are considered to homozygous. It is 0.8 by default.");
        options.addOption("j", true, "Heterozygous ratio (HeR) for variation calling, meaning that the depth of alternative allele is " +
            "greater than HeR and less than (1-HeR) are considered to be hets. It is 0.35 by default.");
        options.addOption("k", true, "Third allele depth ratio (TDR) for variation calling. If the depth of the third allele is greater " +
            "than TDR by the individual coverage, the site will be ignored. Otherwise, the third allele will be considered as sequencing error. It is 0.2 by default.");
        options.addOption("l", true, "Chromosome or region on which genotyping will be performed (e.g. chromosome 1 is designated as 1. " +
            "Region 1bp to 100000bp on chromosome 1 is 1:1,100000)");
        options.addOption("m", true, "Number of threads (taxa number to be processed at the same time). It is 32 by default.");
        options.addOption("n", true, "Individual genotype output directory.");
        options.addOption("o", true, "The path of samtools.");
    }

    @Override
    public void retrieveAppParameters(String[] args) {
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine line = parser.parse(options, args);
            String inOpt = null;
            this.referenceFileS = line.getOptionValue("a");
            this.taxaBamMapFileS = line.getOptionValue("b");
            inOpt = line.getOptionValue("c");
            if (inOpt != null) {
                if (inOpt.equals("1")) {
                    this.baqMode ="";
                }
                else if (inOpt.equals("0")) {
                    this.baqMode = "-B ";
                }
                else {
                    System.out.println("Warning: Incorrect input for -c option. The BAQ computation is disabled by default");
                }
                inOpt = null;
            }
            inOpt = line.getOptionValue("d");
            if (inOpt != null) {
                this.mappingQThresh = Integer.parseInt(inOpt);
                inOpt = null;
            }
            inOpt = line.getOptionValue("e");
            if (inOpt != null) {
                this.baseQThresh = Integer.parseInt(inOpt);
                inOpt = null;
            }
            inOpt = line.getOptionValue("f");
            if (inOpt != null) {
                this.mdcThresh = Integer.parseInt(inOpt);
                inOpt = null;
            }
            inOpt = line.getOptionValue("g");
            if (inOpt != null) {
                this.mindrThresh = Double.parseDouble(inOpt);
                inOpt = null;
            }
            inOpt = line.getOptionValue("h");
            if (inOpt != null) {
                this.maxdrTrresh = Double.parseDouble(inOpt);
                inOpt = null;
            }
            inOpt = line.getOptionValue("i");
            if (inOpt != null) {
                this.horThresh = Double.parseDouble(inOpt);
                inOpt = null;
            }
            inOpt = line.getOptionValue("j");
            if (inOpt != null) {
                this.herThresh = Double.parseDouble(inOpt);
                inOpt = null;
            }
            inOpt = line.getOptionValue("k");
            if (inOpt != null) {
                this.tdrTresh = Double.parseDouble(inOpt);
                inOpt = null;
            }
            String[] tem = line.getOptionValue("l").split(":");
            this.chrom = Integer.parseInt(tem[0]);
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
            inOpt = line.getOptionValue("m");
            if (inOpt != null) {
                this.threadsNum = Integer.parseInt(inOpt);
                inOpt = null;
            }
            this.outputDirS = line.getOptionValue("n");
            this.samtoolsPath = line.getOptionValue("o");
            this.parseTaxaBamMap(this.taxaBamMapFileS);
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
        System.out.println("Below are the commands of module \"disc\" in FastCall 3.");
        this.printUsage();
    }

    private void variationDiscovery () {
        File outDir = new File (this.outputDirS);
        outDir.mkdirs();
        File[] taxaOutDirs = new File[this.taxaNames.length];
        for (int i = 0; i < taxaNames.length; i++) {
            File f = new File (outDir, taxaNames[i]);
            f.mkdir();
            taxaOutDirs[i] = f;
        }
        Dyad<int[][], int[]> d = FastCall3.getBins(this.regionStart, this.regionEnd, FastCall3.disBinSize);
        int[][] binBound = d.getFirstElement();
        int[] binStarts = d.getSecondElement();
        try {
            LongAdder counter = new LongAdder();
            ExecutorService pool = Executors.newFixedThreadPool(this.threadsNum);
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < taxaNames.length; i++) {
                String[] bamFiles = this.taxaBamPathMap.get(taxaNames[i]);
                sb.setLength(0);
                sb.append(this.samtoolsPath).append(" mpileup --no-output-ends ").append(this.baqMode).append("-q ").append(this.mappingQThresh).append(" -Q ").append(this.baseQThresh).append(" -f ").append(this.referenceFileS);
                for (int j = 0; j < bamFiles.length; j++) {
                    sb.append(" ").append(bamFiles[j]);
                }
                sb.append(" -r ").append(this.chrom).append(":").append(this.regionStart).append("-").append(regionEnd-1);
                String command = sb.toString();
//                System.out.println(command);
                TaxonCall tc = new TaxonCall(command, binBound, binStarts, taxaNames[i], taxaOutDirs[i], counter);
                Future<TaxonCall> f = pool.submit(tc);
            }
            pool.shutdown();
            pool.awaitTermination(Long.MAX_VALUE, TimeUnit.SECONDS);
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        System.out.println("Individual genotype of "+ String.valueOf(this.taxaNames.length)+" taxa is completed.");
        System.out.println("Variantion discovery is finished.");
    }

    /**
     * A Callable class for processing individual taxon data during variant calling.
     * <p>
     * This class is responsible for processing sequence data for a single taxon,
     * including base calling, quality filtering, and variant identification.
     * It implements {@link Callable} to enable parallel processing of multiple taxa.
     *
     * <p>Key responsibilities include:
     * <ul>
     *   <li>Executing samtools commands to process BAM files</li>
     *   <li>Processing pileup data to identify variants</li>
     *   <li>Applying quality filters to called variants</li>
     *   <li>Tracking allele counts and variant statistics</li>
     *   <li>Handling both SNP and indel variants</li>
     * </ul>
     *
     * <p>The class maintains several important data structures:
     * <ul>
     *   <li>{@code alleleCount}: Counts of each allele at the current position</li>
     *   <li>{@code insertionLengthSet}, {@code deletionLengthSet}: Track observed indel lengths</li>
     *   <li>{@code indelSeqL}, {@code indelSeq}: Store sequence information for indels</li>
     *   <li>{@code binBound}, {@code binStarts}: Manage genomic binning for efficient processing</li>
     * </ul>
     *
     * @see java.util.concurrent.Callable
     * @see DiscoverVariationF3
     */
    class TaxonCall implements Callable<TaxonCall> {
        String command = null;
        int[][] binBound = null;
        int[] binStarts = null;
        String taxon = null;
        double taxonCoverage;
        File outDir = null;
        LongAdder counter = null;

        int currentPos = Integer.MIN_VALUE;
        int currentDepth = Integer.MIN_VALUE;
        StringBuilder baseSb = new StringBuilder();
        List<String> lList = new ArrayList<>();
        int[] alleleCount = new int[FastCall3.pileupAlleleAscIIs.length];
        IntOpenHashSet insertionLengthSet = new IntOpenHashSet();
        IntOpenHashSet deletionLengthSet = new IntOpenHashSet();
        boolean ifWrite = false;
        byte altAllele = Byte.MIN_VALUE;
        int indelLength = Integer.MIN_VALUE;
        long[] indelSeqL = new long[2];
        String indelSeq = null;
        int altAlleleDepth = Integer.MIN_VALUE;
        DataOutputStream dos = null;
        int currentBinIndex = Integer.MIN_VALUE;

        public TaxonCall (String command, int[][] binBound, int[] binStarts, String taxon, File outDir, LongAdder counter) {
            this.command = command;
            this.binBound = binBound;
            this.binStarts = binStarts;
            this.taxon = taxon;
            this.taxonCoverage = taxaCoverageMap.get(taxon);
            this.outDir = outDir;
            this.counter = counter;
        }

        /**
         * Reset variables for processing a new pileup line.
         * <p>
         * This method is called at the beginning of each pileup line.
         * It resets the StringBuilder {@code baseSb} to empty, resets
         * the current depth to 0, and sets the boolean flag {@code ifWrite}
         * to false. This is necessary to prevent data from previous lines
         * from interfering with the processing of the current line.
         */
        private void initialize1 () {
            this.baseSb.setLength(0);
            this.currentDepth = 0;
            ifWrite = false;
        }

/* <<<<<<<<<<<<<<  ✨ Windsurf Command ⭐ >>>>>>>>>>>>>>>> */
        /**
         * <p>
         * This method is called after the {@code initialize1} method.
         * It resets the following variables:
         * <ul>
         *   <li>{@code alleleCount}: an array for counting the occurrence of each allele.</li>
         *   <li>{@code insertionLengthSet}, {@code deletionLengthSet}: two sets for storing the lengths of observed insertions and deletions.</li>
         *   <li>{@code indelLength}: the length of the current indel.</li>
         * </ul>
         */
/* <<<<<<<<<<  e4719dd8-60ce-4ab5-9805-28bff15e9d7f  >>>>>>>>>>> */
        private void initialize2 () {
            Arrays.fill(alleleCount, 0);
            insertionLengthSet.clear();
            deletionLengthSet.clear();
            indelLength = 0;
        }

        /**
         * Processes a single line from the pileup file and returns true if the site is valid and false otherwise.
         * <p>
         * This method processes a single line from the pileup file. It first checks if the site
         * is at a valid position (i.e. the base is A, C, G, T, or N). It then checks if the
         * site is at a valid depth (i.e. the depth is greater than or equal to the specified
         * minimum depth and less than or equal to the specified maximum depth). It then checks
         * if the site is at a valid depth ratio (i.e. the ratio of the site depth to the
         * taxon coverage is greater than or equal to the specified minimum depth ratio and
         * less than or equal to the specified maximum depth ratio). It then checks if the
         * site is at a valid heterozygous depth ratio (i.e. the ratio of the second highest
         * allele to the site depth is greater than or equal to the specified minimum
         * heterozygous depth ratio and less than or equal to the specified maximum
         * heterozygous depth ratio).
         * <p>
         * If the site is valid, this method sets the alternative allele and its depth and
         * sets the {@code ifWrite} flag to true. If the site is not valid, it sets the
         * {@code ifWrite} flag to false.
         *
         * @param line a single line from the pileup file
         * @return true if the site is valid and false otherwise
         */
        public boolean processPileupLine (String line) {
            lList = PStringUtils.fastSplit(line);
            if (Arrays.binarySearch(BaseEncoder.bases, lList.get(2).charAt(0)) < 0) return false;
            this.initialize1();
            currentPos = Integer.parseInt(lList.get(1));
            for (int i = 3; i < lList.size(); i+=3) {
                currentDepth+=Integer.parseInt(lList.get(i));
                baseSb.append(lList.get(i+1));
            }
            if (currentDepth < mdcThresh) return false;
            if (currentDepth < 1) return false;
            double siteDepthRatio = (double)currentDepth/this.taxonCoverage;
            if (siteDepthRatio < mindrThresh) return false;
            if (siteDepthRatio > maxdrTrresh) return false;
//            FastCall3.removeFirstPositionSign(baseSb);
            String baseS = baseSb.toString().toUpperCase();
            byte[] baseB = baseS.getBytes();
            this.initialize2();
            int index = 0;
            int vCnt = 0;
            for (int i = 0; i < baseB.length; i++) {
                byte alleleCoding = FastCall3.pileupAscIIToAlleleCodingMap.get(baseB[i]);
                index = Arrays.binarySearch(AlleleEncoder.alleleCodings, alleleCoding);
                if (index < 0) continue;
                if (index > 3) {
                    int startIndex = i+1;
                    int endIndex = i+2;
                    for (int j = i+2; j < baseB.length; j++) {
                        if (baseB[j] > 57) {
                            endIndex = j;
                            break;
                        }
                    }
                    baseSb.setLength(0);
                    for (int j = startIndex; j < endIndex; j++) {
                        baseSb.append((char)baseB[j]);
                    }
                    int length = Integer.parseInt(baseSb.toString());
                    if (index == 5) insertionLengthSet.add(length);
                    else deletionLengthSet.add(length);
                    indelSeq = baseS.substring(endIndex, endIndex+length);
                    i+=baseSb.length();
                    i+=length+1;
                }
                alleleCount[index]++;
                vCnt++;
            }
            if (vCnt == 0) return false;
            if (insertionLengthSet.size()+deletionLengthSet.size() > indelTypeThresh) return false;
            int[] alleleCountDesendingIndex = PArrayUtils.getIndicesByDescendingValue(alleleCount);
            double alleleDepthRatio = (double)alleleCount[alleleCountDesendingIndex[0]]/currentDepth;
            if (alleleDepthRatio < herThresh) return false;
            else if (alleleDepthRatio > 1 - herThresh && alleleDepthRatio < horThresh) return false;
            if (alleleCount[alleleCountDesendingIndex[1]] != 0) {
                alleleDepthRatio = (double)alleleCount[alleleCountDesendingIndex[1]]/currentDepth;
                if (alleleDepthRatio > tdrTresh) return false;
            }
            ifWrite = true;
            this.altAllele = AlleleEncoder.alleleCodings[alleleCountDesendingIndex[0]];
            this.altAlleleDepth = alleleCount[alleleCountDesendingIndex[0]];
            if (this.altAllele == AlleleEncoder.alleleCodings[5] && insertionLengthSet.size() > 0) {
                indelLength = insertionLengthSet.toArray(new int[insertionLengthSet.size()])[0];
            }
            else if (this.altAllele == AlleleEncoder.alleleCodings[4] && deletionLengthSet.size() > 0) {
                indelLength = deletionLengthSet.toArray(new int[deletionLengthSet.size()])[0];
            }
            return ifWrite;
        }

        /**
         * Close the DataOutputStream and flush it. This is used to mark the end
         * of a bin in the output file.
         */
        public void closeDos () {
            if (currentBinIndex < 0) return;
            try {
                dos.writeInt(Integer.MAX_VALUE);
                dos.flush();
                dos.close();
            }
            catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        }

        /**
         * Sets the DataOutputStream to write to. This is called by the
         * TaxonCall.writeVariants() method before writing a variant to the
         * output file. This method is used to mark the beginning of a new bin
         * in the output file.
         * @return nothing
         */
        public void setDos () {
            int binIndex = Arrays.binarySearch(binStarts, this.currentPos);
            if (binIndex < 0) binIndex = -binIndex-2;
            if (binIndex != currentBinIndex) {
                if (currentBinIndex > -1) {
                    this.closeDos();
                }
                StringBuilder sb = new StringBuilder();
                sb.append(chrom).append("_").append(binBound[binIndex][0]).append("_").append(binBound[binIndex][1]).append(".ing.gz");
                String outfileS = new File (outDir, sb.toString()).getAbsolutePath();
//                System.out.println(outfileS);
//                dos = IOUtils.getBinaryGzipWriter(outfileS);
                dos = IOUtils.getBigBufferedBinaryGzipWriter(outfileS); // TODO: better buffered for easy write io load.
                try {
                    dos.writeUTF(this.taxon);
                    dos.writeShort((short)chrom);
                    dos.writeInt(binBound[binIndex][0]);
                    dos.writeInt(binBound[binIndex][1]);
                }
                catch (Exception e) {
                    e.printStackTrace();
                    System.exit(1);
                }
                currentBinIndex = binIndex;
            }
        }

        /**
         * Writes the variants to the output file.
         * @return nothing
         */
        public void writeVariants () {
            this.setDos();
            try {
                int allelePack = AllelePackageF3.getAllelePack(binStarts[currentBinIndex], currentPos, altAllele, indelLength);
                dos.writeInt(allelePack);
            }
            catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        }

        /**
         * Executes the command to call the genotypes for a given taxon. This
         * method reads the output from the command, processes the output line by
         * line, and writes the results to the output file. The output file is
         * written in the format specified by the {@link AllelePackageF3}
         * class.
         * @return nothing
         * @throws Exception if there is a problem executing the command
         */
        @Override
        public TaxonCall call() throws Exception {
            Runtime rt = Runtime.getRuntime();
            Process p = rt.exec(command);

            try {
                p.getOutputStream().close();
            } catch (Exception ignore) {}

            Thread errDrainer = new Thread(() -> {
                try (BufferedReader ebr = new BufferedReader(new InputStreamReader(p.getErrorStream()))) {
                    String errLine;
                    while ((errLine = ebr.readLine()) != null) {
                        if (errLine.startsWith("[mpileup]") && errLine.contains("samples in") && errLine.contains("input files")) {
                            continue;
                        }
                        System.err.println("[samtools std err]: " + errLine);
                    }
                } catch (Exception ex) {
                    ex.printStackTrace();
                }
            }, "samtools-stderr-drainer");
            errDrainer.setDaemon(true);
            errDrainer.start();

            try (BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()))) {
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    if(!this.processPileupLine(temp)) continue;
                    this.writeVariants();
                }
                this.closeDos();
            }

            boolean finished = p.waitFor(5, TimeUnit.HOURS);    // set 5 hours for one sample
            if (!finished) {
                p.destroyForcibly();
                throw new RuntimeException("Command timed out: " + command);
            }

            errDrainer.join();

            int exit = p.exitValue();
            if (exit != 0) {
                throw new RuntimeException("Command failed (exit=" + exit + "): " + command);
            }
            System.out.println("Individual genotype is completed for taxon "+ this.taxon);

            counter.increment();
            int count = counter.intValue();
            if (count%50 == 0) {
                System.out.println("Variation calling has been performed for "+ String.valueOf(count)+ " taxa.");
            }
            return null;
        }
    }

        /**
         * Parse the taxa-bam map file. The taxa-bam map file is a tab-delimited
         * text file. Each line represents a taxon, with the first column being
         * the taxon name, the second column being the coverage of the taxon,
         * and the third and beyond columns being the paths to the bam files
         * for the taxon. This method reads the file and stores the data in two
         * maps, one for the bam paths and one for the coverage. The bam paths
         * are stored in a sorted array.
         * @param taxaBamMapFileS the path to the taxa-bam map file
         */
    private void parseTaxaBamMap(String taxaBamMapFileS) {
        this.taxaBamPathMap = new HashMap<>();
        this.taxaCoverageMap = new HashMap<>();
        try {
            BufferedReader br = IOUtils.getTextReader(taxaBamMapFileS);
            String temp = br.readLine();
            ArrayList<String> taxaList = new ArrayList();
            ArrayList<String> pathList = new ArrayList();
            int nBam = 0;
            while ((temp = br.readLine()) != null) {
                String[] tem = temp.split("\t");
                taxaList.add(tem[0]);
                String[] bams = new String[tem.length-2] ;
                for (int i = 0; i < bams.length; i++) {
                    bams[i] = tem[i+2];
                    pathList.add(bams[i]);
                }
                Arrays.sort(bams);
                taxaBamPathMap.put(tem[0], bams);
                taxaCoverageMap.put(tem[0], Double.valueOf(tem[1]));
                nBam+=bams.length;
            }
            taxaNames = taxaList.toArray(new String[taxaList.size()]);
            Arrays.sort(taxaNames);
            HashSet<String> taxaSet = new HashSet<>(taxaList);
            if (taxaSet.size() != taxaNames.length) {
                System.out.println("Taxa names are not unique. Programs quits");
                System.exit(0);
            }
            System.out.println("Created TaxaBamMap from " + taxaBamMapFileS);
            System.out.println("Taxa number:\t"+String.valueOf(taxaNames.length));
            System.out.println("Bam file number in TaxaBamMap:\t"+String.valueOf(nBam));
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
}
