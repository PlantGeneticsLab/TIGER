package pgl.app.fastCall3;

import com.koloboke.collect.map.IntDoubleMap;
import com.koloboke.collect.map.hash.HashIntDoubleMaps;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.io.FileUtils;
import pgl.AppAbstract;
import pgl.PGLAPPEntrance;
import pgl.PGLConstraints;
import pgl.infra.dna.FastaBit;
import pgl.infra.dna.FastaRecordBit;
import pgl.infra.dna.allele.AlleleEncoder;
import pgl.infra.utils.*;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.LongAdder;

import static cern.jet.math.Arithmetic.factorial;

/**
 * A high-performance genotyping module for the FastCall3 pipeline that processes BAM alignments
 * to call genetic variants and generate genotype calls across multiple samples.
 *
 * <p>This class extends {@link AppAbstract} and implements parallel processing of genomic data
 * to efficiently genotype variants from next-generation sequencing data. It supports both single
 * and multi-threaded operations for optimal performance on various hardware configurations.
 *
 * <p>Key features include:
 * <ul>
 *   <li>Parallel processing of multiple samples using thread pools</li>
 *   <li>Support for both single and multi-sample variant calling</li>
 *   <li>Integration with samtools for efficient BAM file processing</li>
 *   <li>Generation of VCF format output with comprehensive variant annotations</li>
 *   <li>Quality filtering based on mapping quality and base quality scores</li>
 *   <li>Support for both SNP and indel variant calling</li>
 * </ul>
 *
 * <p>The genotyping workflow consists of:
 * <ol>
 *   <li>Processing input BAM files and reference genome</li>
 *   <li>Scanning alignments to identify candidate variants</li>
 *   <li>Computing genotype likelihoods and making genotype calls</li>
 *   <li>Generating VCF output with variant and genotype information</li>
 * </ol>
 *
 *
 * @author Fei Lu
 * @version 3.0
 * @since 1.0
 * @see AppAbstract
 * @see IndividualGenotypeF3
 * @see IndividualCountF3
 */
class ScanGenotypeF3 extends AppAbstract {
    //Reference genome file with an index file (.fai). The reference should be in Fasta format. Chromosomes are labled as 1-based numbers (1,2,3,4,5...).
    String referenceFileS = null;
    //The taxaBamMap file contains information of taxon and its corresponding bam files. The bam file should have .bai file in the same folder.
    String taxaBamMapFileS = null;
    //The genetic variation library file
    String libFileS = null;
    int chrom = Integer.MIN_VALUE;
    //Starting position of the specified region for variation calling, inclusive
    int regionStart = Integer.MIN_VALUE;
    //Ending position the specified regionfor variation calling, exclusive
    int regionEnd = Integer.MIN_VALUE;
    //The switch of base alignment quality (BAQ) computaiton, 0 is diabled and 1 is enbabled.
    String baqMode = "-B ";
    //Minimum mapping quality (MQ) for an alignment to be used for genotyping.
    int mappingQThresh = 30;
    //Minimum base quality (BQ) for a base to be used for genotyping.
    int baseQThresh = 20;
    //combined: sequencing error and alignment error
    double combinedErrorRate = 0.05;
    //Minimum read depth ratio (MiDR) relative to the taxon's coverage. Same default as disc.
    // DP=1 is not filtered by an absolute MDC cutoff: at ~3x it is mostly Poisson sampling, not deletion.
    double mindrThresh = 0.2;
    //The path of samtools
    String samtoolsPath = null;
    //VCF output directory
    String outputDirS = null;
    //scan3: previous scan output directory that contains unfiltered indiCounts
    String inputScanDirS = null;
    //scan3: map of merged taxon -> source taxa
    String mergeMapFileS = null;
    //When true, apply MiDR when writing VCF (scan, scan2, scan3, and scan4).
    boolean filterDepthAtGenotype = true;
    String scanModuleName = "scan";
    //Number of threads (taxa number to be processed at the same time)
    int threadsNum = PGLConstraints.parallelLevel;

    String[] subDirS = {"indiVCF", "indiCounts", "VCF"};

    HashMap<String, List<String>> taxaBamsMap = null;
    HashMap<String, Double> taxaCoverageMap = null;
    String[] taxaNames = null;

    IntDoubleMap factorialMap = null;
    int maxFactorial = 150;

    String vLibPosFileS = null;

    FastaBit genomeFa = null;
    int chromIndex = Integer.MIN_VALUE;
    VariationLibraryF3 vl = null;
    int vlStartIndex = Integer.MIN_VALUE;
    int vlEndIndex = Integer.MIN_VALUE;
    HashMap<Integer, String> posRefMap = new HashMap<>();
    //HashMap would be faster to locate alleles, larger memory though
    HashMap<Integer, int[]> posAllelePackMap = new HashMap<>();
    int[] positions = null;
    int vlBinStartIndex = 0;
    int vlBinEndIndex = 0;

    public ScanGenotypeF3(String[] args) {
        this.creatAppOptions();
        this.retrieveAppParameters(args);

        this.mkDir();
        this.processVariationLibrary();
        this.creatFactorialMap();

        /*
        Output by individual allele count, fast
         */
        this.scanIndiCountsByThreadPool();
        this.mkFinalVCFFromIndiCounts();
    }

    public ScanGenotypeF3(String[] args, String stepFlag) {
        if ("iac".equals(stepFlag)) {
            this.runScanIac(args);
            return;
        }
        if ("3".equals(stepFlag)) {
            this.runScan3(args);
            return;
        }
        if ("2".equals(stepFlag) || "4".equals(stepFlag)) {
            this.creatAppOptions();
            this.retrieveStep2Parameters(args);
            this.scanModuleName = "4".equals(stepFlag) ? "scan4" : "scan2";

            this.processVariationLibrary();
            this.creatFactorialMap();

            FastaRecordBit frb = genomeFa.getFastaRecordBit(chromIndex);
            positions = new int[vlEndIndex-vlStartIndex];
            for (int i = vlStartIndex; i < vlEndIndex; i++) {
                posRefMap.put(vl.positions[i], String.valueOf(frb.getBase(vl.positions[i]-1)));
                posAllelePackMap.put(vl.positions[i], vl.getAllelePacks(i));
                positions[i-vlStartIndex] = vl.positions[i];
            }

            if ("4".equals(stepFlag)) {
                this.mkFinalVCFStreamingFromIndiCounts();
            }
            else {
                this.mkFinalVCFFromIndiCounts();
            }
            return;
        }
        new ScanGenotypeF3(args);
    }

    /**
     * Same pileup counting as {@code scan}, but stops after writing {@code indiCounts/*.iac.gz}.
     */
    private void runScanIac(String[] args) {
        this.scanModuleName = "scan_iac";
        this.creatAppOptions();
        this.retrieveAppParameters(args);
        this.mkDir();
        this.processVariationLibrary();
        this.scanIndiCountsByThreadPool();
        System.err.println("Individual allele counts written to " + new File(this.outputDirS, subDirS[1]).getAbsolutePath());
    }

    private void runScan3(String[] args) {
        this.scanModuleName = "scan3";
        this.filterDepthAtGenotype = true;
        this.creatAppOptions();
        this.retrieveScan3Parameters(args);
        new File(this.outputDirS).mkdirs();
        this.processVariationLibrary();
        this.creatFactorialMap();
        FastaRecordBit frb = genomeFa.getFastaRecordBit(chromIndex);
        positions = new int[vlEndIndex-vlStartIndex];
        for (int i = vlStartIndex; i < vlEndIndex; i++) {
            posRefMap.put(vl.positions[i], String.valueOf(frb.getBase(vl.positions[i]-1)));
            posAllelePackMap.put(vl.positions[i], vl.getAllelePacks(i));
            positions[i-vlStartIndex] = vl.positions[i];
        }
        this.mkFinalVCFStreamingFromMergedIndiCounts();
    }

    public void retrieveStep2Parameters(String[] args) {
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine line = parser.parse(options, args);
            String inOpt = null;

            this.outputDirS = line.getOptionValue("k");

            String[] tem = line.getOptionValue("d").split(":");
            this.chrom = Integer.parseInt(tem[0]);
            this.referenceFileS = line.getOptionValue("a");
            long start = System.nanoTime();
            System.err.println("Reading reference genome from "+ referenceFileS);
            this.genomeFa = new FastaBit(referenceFileS);
            System.err.println("Reading reference genome took " + String.format("%.2f", Benchmark.getTimeSpanSeconds(start)) + "s");
            this.chromIndex = genomeFa.getIndexByDescription(String.valueOf(this.chrom));
            if (tem.length == 1) {
                this.regionStart = 1;
                this.regionEnd = genomeFa.getSeqLength(chromIndex)+1;
            }
            else if (tem.length == 2) {
                tem = tem[1].split(",");
                this.regionStart = Integer.parseInt(tem[0]);
                this.regionEnd = Integer.parseInt(tem[1])+1;
            }

            this.taxaBamMapFileS = line.getOptionValue("b");

            this.libFileS = line.getOptionValue("c");

            inOpt = line.getOptionValue("j");
            if (inOpt != null) {
                this.threadsNum = Integer.parseInt(inOpt);
                inOpt = null;
            }
            inOpt = line.getOptionValue("h");
            if (inOpt != null) {
                this.combinedErrorRate = Double.parseDouble(inOpt);
                inOpt = null;
            }
            this.parseDepthCoverageOptions(line);
        }
        catch(Exception e) {
            e.printStackTrace();
            System.err.println("\nThere are input errors in the command line. Program stops.");
            this.printInstructionAndUsage();
            System.exit(0);
        }
        this.parseTaxaBamMap(this.taxaBamMapFileS);
    }

    private void retrieveScan3Parameters(String[] args) {
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine line = parser.parse(options, args);
            this.referenceFileS = line.getOptionValue("a");
            this.mergeMapFileS = line.getOptionValue("p");
            this.libFileS = line.getOptionValue("c");
            this.inputScanDirS = line.getOptionValue("o");
            this.outputDirS = line.getOptionValue("k");
            if (this.mergeMapFileS == null || this.inputScanDirS == null || this.outputDirS == null) {
                throw new IllegalArgumentException("scan3 requires -o (input scan dir), -p (merge map), and -k (output dir).");
            }
            String[] tem = line.getOptionValue("d").split(":");
            this.chrom = Integer.parseInt(tem[0]);
            long start = System.nanoTime();
            System.err.println("Reading reference genome from "+ referenceFileS);
            this.genomeFa = new FastaBit(referenceFileS);
            System.err.println("Reading reference genome took " + String.format("%.2f", Benchmark.getTimeSpanSeconds(start)) + "s");
            this.chromIndex = genomeFa.getIndexByDescription(String.valueOf(this.chrom));
            if (tem.length == 1) {
                this.regionStart = 1;
                this.regionEnd = genomeFa.getSeqLength(chromIndex)+1;
            }
            else if (tem.length == 2) {
                tem = tem[1].split(",");
                this.regionStart = Integer.parseInt(tem[0]);
                this.regionEnd = Integer.parseInt(tem[1])+1;
            }
            String inOpt = line.getOptionValue("h");
            if (inOpt != null) {
                this.combinedErrorRate = Double.parseDouble(inOpt);
            }
            inOpt = line.getOptionValue("j");
            if (inOpt != null) {
                this.threadsNum = Integer.parseInt(inOpt);
            }
            this.parseDepthCoverageOptions(line);
            this.taxaBamMapFileS = line.getOptionValue("b");
            if (this.taxaBamMapFileS != null) {
                this.parseTaxaBamMap(this.taxaBamMapFileS);
            }
            this.parseMergeMap(this.mergeMapFileS);
        }
        catch(Exception e) {
            e.printStackTrace();
            System.err.println("\nThere are input errors in the command line. Program stops.");
            this.printInstructionAndUsage();
            System.exit(0);
        }
    }

    /**
     * taxaCombineMap columns: mergedTaxon, coverageToken, sourceTaxon1, sourceTaxon2, ...
     * A taxon that is not merged is listed as its own single source
     * (mergedTaxon\t.\tmergedTaxon), so the same file can drive both the merged
     * columns and the pass-through columns. Coverage is the sum of the source
     * coverages from taxaBamMap (-b). A numeric second column is still accepted
     * as an explicit coverage.
     */
    private LinkedHashMap<String, List<String>> mergedSourceMap = new LinkedHashMap<>();
    //source taxon -> the merged taxon that consumed it, used to catch double counting
    private HashMap<String, String> sourceOwnerMap = new HashMap<>();

    private void parseMergeMap(String mergeMapFileS) {
        HashMap<String, Double> sourceCoverageMap = this.taxaCoverageMap;
        this.taxaCoverageMap = new HashMap<>();
        this.mergedSourceMap.clear();
        this.sourceOwnerMap.clear();
        try {
            BufferedReader br = IOUtils.getTextReader(mergeMapFileS);
            String temp = br.readLine();
            if (temp == null) {
                throw new IllegalArgumentException("Merge map is empty: " + mergeMapFileS);
            }
            if (!this.isMergeMapHeader(temp)) {
                this.addMergeMapLine(temp, sourceCoverageMap);
            }
            while ((temp = br.readLine()) != null) {
                if (temp.trim().isEmpty() || this.isMergeMapHeader(temp)) continue;
                this.addMergeMapLine(temp, sourceCoverageMap);
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        this.taxaNames = this.mergedSourceMap.keySet().toArray(new String[0]);
        Arrays.sort(this.taxaNames);
        System.err.println("scan3 merge groups:\t" + this.taxaNames.length);
        for (String name : this.taxaNames) {
            System.err.println("  " + name + " coverage=" + this.taxaCoverageMap.get(name)
                    + " sources=" + this.mergedSourceMap.get(name));
        }
    }

    private boolean isMergeMapHeader(String line) {
        if (line.startsWith("#")) return true;
        String first = line.split("\t", 2)[0].trim().toLowerCase();
        return first.equals("taxa") || first.equals("merged") || first.equals("mergedtaxa") || first.startsWith("merged");
    }

    /**
     * A numeric or '.' second column is coverage (old merge map). Otherwise every column
     * after the merged name is a source taxon, and coverage is the sum from -b.
     */
    private boolean isCoverageToken(String token) {
        if (token == null) return true;
        String t = token.trim();
        if (t.isEmpty() || t.equals(".")) return true;
        try {
            Double.parseDouble(t);
            return true;
        }
        catch (NumberFormatException e) {
            return false;
        }
    }

    private void addMergeMapLine(String line, HashMap<String, Double> sourceCoverageMap) {
        String[] tem = line.split("\t");
        if (tem.length < 2) {
            throw new IllegalArgumentException("Merge map line needs a merged taxon and at least one source taxon:\n" + line);
        }
        String merged = tem[0].trim();
        int sourceStart = this.isCoverageToken(tem[1]) ? 2 : 1;
        List<String> sources = new ArrayList<>();
        for (int i = sourceStart; i < tem.length; i++) {
            if (tem[i].trim().isEmpty()) continue;
            String src = tem[i].trim();
            if (sources.contains(src)) {
                throw new IllegalArgumentException("Source taxon " + src + " is listed more than once for merged taxon "
                        + merged + ", which would double its allele counts:\n" + line);
            }
            sources.add(src);
        }
        if (sources.isEmpty()) {
            throw new IllegalArgumentException("Merged taxon " + merged + " needs at least one source taxon:\n" + line);
        }
        double coverage;
        if (sourceStart == 1 || tem[1].trim().isEmpty() || tem[1].trim().equals(".")) {
            if (sourceCoverageMap == null) {
                throw new IllegalArgumentException("Merged taxon " + merged + " needs -b taxaBamMap to sum source coverages.");
            }
            coverage = 0;
            for (String src : sources) {
                Double c = sourceCoverageMap.get(src);
                if (c == null) {
                    throw new IllegalArgumentException("Source taxon " + src + " not found in taxaBamMap.");
                }
                coverage += c;
            }
        }
        else {
            coverage = Double.parseDouble(tem[1].trim());
        }
        for (String src : sources) {
            String owner = this.sourceOwnerMap.put(src, merged);
            if (owner != null && !owner.equals(merged)) {
                System.err.println("Warning: source taxon " + src + " is used by both " + owner
                        + " and " + merged + "; it will appear in two output columns.");
            }
        }
        this.mergedSourceMap.put(merged, sources);
        this.taxaCoverageMap.put(merged, coverage);
    }

    /**
     * scan3: reads the unfiltered source .iac files a window of sites at a time, sums the
     * allele counts of each merge group on the fly, and writes the merged genotypes to
     * stdout. No merged .iac is written; the source indiCounts/ stays untouched.
     */
    private void mkFinalVCFStreamingFromMergedIndiCounts () {
        final int window = FastCall3.scanVcfWindowSites;
        IndividualCountF3.SiteStream[] streams = null;
        ExecutorService readPool = null;
        try {
            SimpleDateFormat sdf = new SimpleDateFormat("MM/dd/yyyy HH:mm:ss.SSS");
            Date dt = new Date();
            String S = sdf.format(dt);
            OutputStream stdoutNoClose = new FilterOutputStream(System.out) {
                @Override
                public void close() throws IOException {
                    flush();
                }
            };
            BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(stdoutNoClose, StandardCharsets.UTF_8), 65536);
            bw.write("##fileformat=VCFv4.1\n");
            bw.write("##fileDate="+S.split(" ")[0]+"\n");
            bw.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
            bw.write("##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the reference and alternate alleles in the order listed\">\n");
            bw.write("##FORMAT=<ID=GL,Number=G,Type=Integer,Description=\"Genotype likelihoods for 0/0, 0/1, 1/1, or  0/0, 0/1, 0/2, 1/1, 1/2, 2/2 if 2 alt alleles\">\n");
            bw.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n");
            bw.write("##INFO=<ID=NZ,Number=1,Type=Integer,Description=\"Number of taxa with called genotypes\">\n");
            bw.write("##INFO=<ID=AD,Number=.,Type=Integer,Description=\"Total allelelic depths in order listed starting with REF\">\n");
            bw.write("##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Numbers of allele occurence across taxa in order listed\">\n");
            bw.write("##INFO=<ID=IS,Number=.,Type=Integer,Description=\"Indel sequence of ALT alleles in order listed\">\n");
            bw.write("##INFO=<ID=GN,Number=.,Type=Integer,Description=\"Number of taxa with genotypes AA,AB,BB or AA,AB,AC,BB,BC,CC if 2 alt alleles\">\n");
            bw.write("##INFO=<ID=HT,Number=1,Type=Integer,Description=\"Number of heterozygotes\">\n");
            bw.write("##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Minor allele frequency\">\n");
            bw.write("##ALT=<ID=D,Description=\"Deletion\">\n");
            bw.write("##ALT=<ID=I,Description=\"Insertion\">\n");
            Dyad<int[][], int[]> d = FastCall3.getBins(this.regionStart, this.regionEnd, FastCall3.scanBinSize);
            int[][] binBound = d.getFirstElement();
            int[] binStarts = d.getSecondElement();
            StringBuilder sb = new StringBuilder("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT");
            for (int i = 0; i < taxaNames.length; i++) {
                sb.append("\t").append(taxaNames[i]);
            }
            bw.write(sb.toString());
            bw.newLine();

            File inCounts = new File(this.inputScanDirS, subDirS[1]);
            if (!inCounts.isDirectory()) {
                throw new IllegalStateException("scan3 input indiCounts directory does not exist: " + inCounts.getAbsolutePath());
            }
            String inCountsFolderS = inCounts.getAbsolutePath();
            // A source taxon can feed more than one merge group; read each source once.
            List<String> sourceNames = new ArrayList<>();
            HashMap<String, Integer> sourceIndexMap = new HashMap<>();
            for (int t = 0; t < taxaNames.length; t++) {
                for (String src : this.mergedSourceMap.get(taxaNames[t])) {
                    if (!sourceIndexMap.containsKey(src)) {
                        sourceIndexMap.put(src, sourceNames.size());
                        sourceNames.add(src);
                    }
                }
            }
            int nTaxa = taxaNames.length;
            int[][] groupSourceIndices = new int[nTaxa][];
            for (int t = 0; t < nTaxa; t++) {
                List<String> sources = this.mergedSourceMap.get(taxaNames[t]);
                int[] indices = new int[sources.size()];
                for (int s = 0; s < indices.length; s++) {
                    indices[s] = sourceIndexMap.get(sources.get(s));
                }
                groupSourceIndices[t] = indices;
            }
            int nSources = sourceNames.size();
            int readThreads = Math.max(1, Math.min(this.threadsNum, nSources));
            readPool = Executors.newFixedThreadPool(readThreads);
            int offset = 0;
            for (int i = 0; i < binStarts.length; i++) {
                sb.setLength(0);
                sb.append(chrom).append("_").append(binBound[i][0]).append("_").append(binBound[i][1]).append(".iac.gz");
                String iacName = sb.toString();
                streams = new IndividualCountF3.SiteStream[nSources];
                int positionNum = -1;
                int nMissing = 0;
                for (int s = 0; s < nSources; s++) {
                    String fileS = new File(new File(inCountsFolderS, sourceNames.get(s)), iacName).getAbsolutePath();
                    IndividualCountF3.SiteStream st = IndividualCountF3.SiteStream.open(fileS);
                    if (st == null) {
                        System.err.println("Warning: "+ fileS +" does not exist");
                        nMissing++;
                        continue;
                    }
                    if (positionNum < 0) {
                        positionNum = st.positionNum;
                    }
                    else if (st.positionNum != positionNum) {
                        throw new IllegalStateException("Site-count mismatch in " + fileS + ": "
                                + st.positionNum + " vs " + positionNum);
                    }
                    streams[s] = st;
                }
                if (positionNum < 0) {
                    // No source file for this bin at all: keep the positions, call every taxon missing.
                    positionNum = this.getBinSiteNum(binBound[i][0], binBound[i][1]);
                }
                if (offset + positionNum > positions.length) {
                    throw new ArrayIndexOutOfBoundsException("iac sites " + positionNum
                            + " at offset " + offset + " exceed library positions " + positions.length);
                }

                final IndividualCountF3.SiteStream[] binStreams = streams;
                int remaining = positionNum;
                int localOffset = 0;
                while (remaining > 0) {
                    final int w = Math.min(window, remaining);
                    final short[][][] sourceCounts = new short[w][nSources][];
                    List<Future<?>> reads = new ArrayList<>(nSources);
                    for (int s = 0; s < nSources; s++) {
                        final int source = s;
                        reads.add(readPool.submit(() -> {
                            IndividualCountF3.SiteStream st = binStreams[source];
                            for (int k = 0; k < w; k++) {
                                sourceCounts[k][source] = (st == null) ? null : st.readSite();
                            }
                            return null;
                        }));
                    }
                    for (int r = 0; r < reads.size(); r++) {
                        reads.get(r).get();
                    }

                    final int binOffset = offset + localOffset;
                    final String[] vcfRecords = new String[w];
                    List<Integer> indexList = PArrayUtils.getIndexList(w);
                    indexList.parallelStream().forEach(index -> {
                        int posIdx = binOffset + index;
                        StringBuilder vsb = new StringBuilder();
                        int currentPosition = positions[posIdx];
                        vsb.append(chrom).append("\t").append(currentPosition).append("\t").append(chrom).append("-").append(currentPosition)
                                .append("\t").append(posRefMap.get(currentPosition)).append("\t");
                        int[] altAlleles = posAllelePackMap.get(currentPosition);
                        for (int j = 0; j < altAlleles.length; j++) {
                            vsb.append(AllelePackageF3.getAlleleBase(altAlleles[j])).append(",");
                        }
                        vsb.deleteCharAt(vsb.length()-1).append("\t.\t.\t");
                        List<short[]> siteCountsList = new ArrayList<>(nTaxa);
                        for (int t = 0; t < nTaxa; t++) {
                            siteCountsList.add(sumSourceAlleleCounts(sourceCounts[index], groupSourceIndices[t]));
                        }
                        vsb.append(this.getInfoAndGenotypes(siteCountsList, altAlleles));
                        vcfRecords[index] = vsb.toString();
                    });
                    for (int j = 0; j < w; j++) {
                        bw.write(vcfRecords[j]);
                        bw.newLine();
                    }
                    remaining -= w;
                    localOffset += w;
                }
                closeSiteStreams(streams);
                streams = null;
                offset += positionNum;
                if (nMissing > 0) {
                    System.err.println("Missing source iac files in this bin:\t" + nMissing);
                }
                System.err.println(iacName.split("\\.")[0]+ " genotyping is finished.");
            }
            readPool.shutdown();
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            closeSiteStreams(streams);
            if (readPool != null) {
                readPool.shutdownNow();
            }
            System.exit(1);
        }
        System.err.println("Final VCF written to stdout");
        System.err.println("Genotyping is finished.");
    }

    /**
     * Sums the allele counts of one merge group at a single site. A missing source
     * contributes zero; the site stays missing only when every source is missing.
     */
    private static short[] sumSourceAlleleCounts(short[][] sourceCounts, int[] groupSourceIndices) {
        short[] sum = null;
        for (int s = 0; s < groupSourceIndices.length; s++) {
            short[] counts = sourceCounts[groupSourceIndices[s]];
            if (counts == null) continue;
            if (sum == null) {
                sum = (groupSourceIndices.length == 1) ? counts : counts.clone();
                continue;
            }
            if (sum.length != counts.length) {
                throw new IllegalArgumentException("Allele-number mismatch when merging a site: "
                        + sum.length + " vs " + counts.length);
            }
            for (int j = 0; j < sum.length; j++) {
                int v = sum[j] + counts[j];
                sum[j] = (v > Short.MAX_VALUE) ? Short.MAX_VALUE : (short) v;
            }
        }
        return sum;
    }

    /**
     * Number of library sites in a bin, as recorded in the .iac header by scan.
     */
    private int getBinSiteNum(int binStart, int binEnd) {
        int startIdx = vl.getStartIndex(binStart);
        int endIdx = vl.getEndIndex(binEnd - 1);
        return (startIdx == Integer.MIN_VALUE || endIdx == Integer.MIN_VALUE) ? 0 : Math.max(0, endIdx - startIdx);
    }

    @Override
    public void creatAppOptions() {
        options.addOption("app", true, "App name.");
        options.addOption("mod", true, "Module name of FastCall 3.");
        options.addOption("a", true, "Reference genome file with an index file (.fai). The reference should be in Fasta format. " +
            "Chromosomes are labelled as numbers (1,2,3,4,5...). It is recommended to use reference chromosome while perform genotyping for " +
            "each chromosome because loading reference genome would be much faster.");
        options.addOption("b", true, "The taxaBamMap file contains information of taxon and its corresponding bam files. " +
            "The bam file should have .bai file in the same folder. For scan3, this optional file supplies source-taxon " +
            "coverage when the merge map uses '.' as the coverage column.");
        options.addOption("c", true, "The genetic variation library file.");
        options.addOption("d", true, "Chromosome or region on which genotyping will be performed (e.g. chromosome 1 is designated as 1. " +
            "Region 1bp to 100000bp on chromosome 1 is 1:1,100000)");
        options.addOption("e", true, "The switch of base alignment quality (BAQ) computation, 0 is disabled and 1 is enabled. It is 0 by default.");
        options.addOption("f", true, "Minimum mapping quality (MQ) for an alignment to be used for genotyping. It is 30 by default.");
        options.addOption("g", true, "Minimum base quality (BQ) for a base to be used for genotyping. It is 20 by default.");
        options.addOption("h", true, "Combined error rate of sequencing and misalignment. Heterozygous read mapping are more " +
            "likely to be genotyped as homozygote when the combined error rate is high.");
        options.addOption("i", true, "The path of samtools.");
        options.addOption("j", true, "Number of threads. It is 32 by default.");
        options.addOption("k", true, "Output directory for scan files (pos.txt; indiCounts only for scan and scan_iac). " +
            "scan3 reads its input from -o and writes no iac. The VCF is written to stdout.");
        options.addOption("m", true, "Minimum read depth ratio (MiDR) for genotyping in scan/scan2/scan3/scan4. Sites with depth lower than " +
            "the MiDR of the individual sequencing coverage are treated as missing. It is 0.2 by default. " +
            "There is no maximum depth filter (MaDR) and no absolute minimum depth (MDC) at scan; abnormally high depth is kept.");
        options.addOption("o", true, "scan3: input directory of a previous scan run (must contain indiCounts/).");
        options.addOption("p", true, "scan3: taxaCombineMap. Each line: mergedTaxon\\tsourceTaxon1[\\tsourceTaxon2...]. " +
            "At least one source taxon, so a taxon that is not merged can be listed as its own single source " +
            "(mergedTaxon\\t.\\tmergedTaxon) and passes through unchanged. Coverage is the sum of those sources from -b. " +
            "An optional numeric coverage may still be given as the second column.");
    }

    @Override
    public void retrieveAppParameters(String[] args) {
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine line = parser.parse(options, args);
            String inOpt = null;
            this.referenceFileS = line.getOptionValue("a");
            this.taxaBamMapFileS = line.getOptionValue("b");
            this.libFileS = line.getOptionValue("c");
            String[] tem = line.getOptionValue("d").split(":");
            this.chrom = Integer.parseInt(tem[0]);
            long start = System.nanoTime();
            System.err.println("Reading reference genome from "+ referenceFileS);
            this.genomeFa = new FastaBit(referenceFileS);
            System.err.println("Reading reference genome took " + String.format("%.2f", Benchmark.getTimeSpanSeconds(start)) + "s");
            this.chromIndex = genomeFa.getIndexByDescription(String.valueOf(this.chrom));
            if (tem.length == 1) {
                this.regionStart = 1;
                this.regionEnd = genomeFa.getSeqLength(chromIndex)+1;
            }
            else if (tem.length == 2) {
                tem = tem[1].split(",");
                this.regionStart = Integer.parseInt(tem[0]);
                this.regionEnd = Integer.parseInt(tem[1])+1;
            }
            inOpt = line.getOptionValue("e");
            if (inOpt != null) {
                if (inOpt.equals("1")) {
                    this.baqMode ="";
                }
                else if (inOpt.equals("0")) {
                    this.baqMode = "-B ";
                }
                else {
                    System.err.println("Warning: Incorrect input for -c option. The BAQ computation is disabled by default");
                }
                inOpt = null;
            }
            inOpt = line.getOptionValue("f");
            if (inOpt != null) {
                this.mappingQThresh = Integer.parseInt(inOpt);
                inOpt = null;
            }
            inOpt = line.getOptionValue("g");
            if (inOpt != null) {
                this.baseQThresh = Integer.parseInt(inOpt);
                inOpt = null;
            }
            inOpt = line.getOptionValue("h");
            if (inOpt != null) {
                this.combinedErrorRate = Double.parseDouble(inOpt);
                inOpt = null;
            }
            this.samtoolsPath = line.getOptionValue("i");
            inOpt = line.getOptionValue("j");
            if (inOpt != null) {
                this.threadsNum = Integer.parseInt(inOpt);
                inOpt = null;
            }
            this.outputDirS = line.getOptionValue("k");
            this.parseDepthCoverageOptions(line);
        }
        catch(Exception e) {
            e.printStackTrace();
            System.err.println("\nThere are input errors in the command line. Program stops.");
            this.printInstructionAndUsage();
            System.exit(0);
        }
        this.parseTaxaBamMap(this.taxaBamMapFileS);
    }

    /**
     * Parse optional MiDR. Default 0.2, same as disc. Scan does not apply MaDR or MDC.
     */
    private void parseDepthCoverageOptions(CommandLine line) {
        String inOpt = line.getOptionValue("m");
        if (inOpt != null) {
            this.mindrThresh = Double.parseDouble(inOpt);
        }
    }

    /**
     * Return true if site depth passes MiDR relative to the taxon's mean coverage.
     * Depth=0 is already missing upstream. High depth is kept: disc already filtered
     * discovery, and reads matching the library are treated as valid even if unusually deep.
     */
    private boolean passDepthCoverageFilter(int siteDepth, double taxonCoverage) {
        if (taxonCoverage <= 0) return true;
        double siteDepthRatio = (double) siteDepth / taxonCoverage;
        if (siteDepthRatio < this.mindrThresh) return false;
        return true;
    }

    @Override
    public void printInstructionAndUsage() {
        System.err.println(PGLAPPEntrance.getTIGERIntroduction());
        System.err.println("Below are the commands of module \"" + this.scanModuleName + "\" in FastCall 3.");
        this.printUsage();
    }

/* <<<<<<<<<<<<<<  ✨ Windsurf Command ⭐ >>>>>>>>>>>>>>>> */
    /**
     * Generates a final VCF file from individual count files.
     *
     * @author Fei Lu
     */
/* <<<<<<<<<<  3c9f6e1b-fde9-4fc1-855d-77758a6cbaf1  >>>>>>>>>>> */
    public void mkFinalVCFFromIndiCounts () {
        try {
            SimpleDateFormat sdf = new SimpleDateFormat("MM/dd/yyyy HH:mm:ss.SSS");
            Date dt = new Date();
            String S = sdf.format(dt);
            OutputStream stdoutNoClose = new FilterOutputStream(System.out) {
                @Override
                public void close() throws IOException {
                    flush();
                }
            };
            BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(stdoutNoClose, StandardCharsets.UTF_8), 65536);
            bw.write("##fileformat=VCFv4.1\n");
            bw.write("##fileDate="+S.split(" ")[0]+"\n");
            bw.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
            bw.write("##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the reference and alternate alleles in the order listed\">\n");
            bw.write("##FORMAT=<ID=GL,Number=G,Type=Integer,Description=\"Genotype likelihoods for 0/0, 0/1, 1/1, or  0/0, 0/1, 0/2, 1/1, 1/2, 2/2 if 2 alt alleles\">\n");
            bw.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n");
            bw.write("##INFO=<ID=NZ,Number=1,Type=Integer,Description=\"Number of taxa with called genotypes\">\n");
            bw.write("##INFO=<ID=AD,Number=.,Type=Integer,Description=\"Total allelelic depths in order listed starting with REF\">\n");
            bw.write("##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Numbers of allele occurence across taxa in order listed\">\n");
            bw.write("##INFO=<ID=IS,Number=.,Type=Integer,Description=\"Indel sequence of ALT alleles in order listed\">\n");
            bw.write("##INFO=<ID=GN,Number=.,Type=Integer,Description=\"Number of taxa with genotypes AA,AB,BB or AA,AB,AC,BB,BC,CC if 2 alt alleles\">\n");
            bw.write("##INFO=<ID=HT,Number=1,Type=Integer,Description=\"Number of heterozygotes\">\n");
            bw.write("##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Minor allele frequency\">\n");
            bw.write("##ALT=<ID=D,Description=\"Deletion\">\n");
            bw.write("##ALT=<ID=I,Description=\"Insertion\">\n");
            Dyad<int[][], int[]> d = FastCall3.getBins(this.regionStart, this.regionEnd, FastCall3.scanBinSize);
            int[][] binBound = d.getFirstElement();
            int[] binStarts = d.getSecondElement();
                StringBuilder sb = new StringBuilder("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT");
            for (int i = 0; i < taxaNames.length; i++) {
                sb.append("\t").append(taxaNames[i]);
            }
            bw.write(sb.toString());
            bw.newLine();
            List<Future<IndividualCountF3>> futureList = new ArrayList<>();
            List<IndividualCountF3> incList = new ArrayList<>();
//            vlBinStartIndex = 0;
//            vlBinEndIndex = 0;
            int offset = 0;
            for (int i = 0; i < binStarts.length; i++) {
                futureList.clear();
                incList.clear();
                String indiCountFolderS = new File(outputDirS, subDirS[1]).getAbsolutePath();
                try {
                    LongAdder counter = new LongAdder();
                    ExecutorService pool = Executors.newFixedThreadPool(this.threadsNum);
                    sb.setLength(0);
                    sb.append(chrom).append("_").append(binBound[i][0]).append("_").append(binBound[i][1]).append(".iac.gz");
                    for (int j = 0; j < taxaNames.length; j++) {
                        String indiTaxonDirS = new File (indiCountFolderS, taxaNames[j]).getAbsolutePath();
                        String fileS = new File (indiTaxonDirS, sb.toString()).getAbsolutePath();
                        TaxonCountRead tr = new TaxonCountRead(fileS);
                        Future<IndividualCountF3> result = pool.submit(tr);
                        futureList.add(result);
                    }
                    pool.shutdown();
                    pool.awaitTermination(Long.MAX_VALUE, TimeUnit.SECONDS);
                    for (int j = 0; j < futureList.size(); j++) {
                        IndividualCountF3 inc = futureList.get(j).get();
                        if (inc == null) continue;
                        incList.add(inc);
                    }
                    Collections.sort(incList);
                }
                catch (Exception e) {
                    e.printStackTrace();
                    System.exit(1);
                }

                HashMap<String, IndividualCountF3> incByTaxon = new HashMap<>();
                for (int j = 0; j < incList.size(); j++) {
                    incByTaxon.put(incList.get(j).taxonName, incList.get(j));
                }
                final HashMap<String, IndividualCountF3> incMap = incByTaxon;

                final int binSiteCount = incList.get(0).alleleNum.length;
                final int binOffset = offset;
                List<Integer> indexList = PArrayUtils.getIndexList(binSiteCount);
//                vlBinEndIndex = vlBinStartIndex + incList.get(0).alleleNum.length;
//                List<Integer> indexList = new ArrayList<>();
//                for (int j = 0; j < incList.get(0).alleleNum.length; j++) {
//                    indexList.add(j);
//                }

                String[] vcfRecords = new String[indexList.size()];
                indexList.parallelStream().forEach(index -> {
                    int posIdx = binOffset + index;
                    if (posIdx < 0 || posIdx >= positions.length) {
                        throw new ArrayIndexOutOfBoundsException("positions index " + posIdx + " out of " + positions.length);
                    }
                    StringBuilder vsb = new StringBuilder();
                    int currentPosition = positions[posIdx];
//                    StringBuilder vsb = new StringBuilder();
//                    int currentPosition = positions[index+vlBinStartIndex];

                    vsb.append(chrom).append("\t").append(currentPosition).append("\t").append(chrom).append("-").append(currentPosition)
                            .append("\t").append(posRefMap.get(currentPosition)).append("\t");
                    int[] altAlleles = posAllelePackMap.get(currentPosition);
                    for (int j = 0; j < altAlleles.length; j++) {
                        vsb.append(AllelePackageF3.getAlleleBase(altAlleles[j])).append(",");
                    }
                    vsb.deleteCharAt(vsb.length()-1).append("\t.\t.\t");
                    List<short[]> siteCountsList = new ArrayList<>();
                    for (int j = 0; j < taxaNames.length; j++) {
                        IndividualCountF3 inc = incMap.get(taxaNames[j]);
                        siteCountsList.add(inc == null ? null : inc.alleleCounts[index]);
                    }
                    vsb.append(this.getInfoAndGenotypes(siteCountsList, altAlleles));
                    vcfRecords[index] = vsb.toString();
                });
                for (int j = 0; j < vcfRecords.length; j++) {
                    bw.write(vcfRecords[j]);
                    bw.newLine();
                }

                offset += binSiteCount;
//                vlBinStartIndex = vlBinEndIndex;

                System.err.println(sb.toString().split("\\.")[0]+ " genotyping is finished.");
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
//        this.deleteTemperateFile();
        System.err.println("Final VCF written to stdout");
        System.err.println("Genotyping is finished.");
    }

    /**
     * Same genotypes as {@link #mkFinalVCFFromIndiCounts()}, but reads {@code .iac.gz}
     * a window of sites at a time, writes those VCF lines, then discards the AD.
     */
    public void mkFinalVCFStreamingFromIndiCounts () {
        final int window = FastCall3.scanVcfWindowSites;
        IndividualCountF3.SiteStream[] streams = null;
        try {
            SimpleDateFormat sdf = new SimpleDateFormat("MM/dd/yyyy HH:mm:ss.SSS");
            Date dt = new Date();
            String S = sdf.format(dt);
            OutputStream stdoutNoClose = new FilterOutputStream(System.out) {
                @Override
                public void close() throws IOException {
                    flush();
                }
            };
            BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(stdoutNoClose, StandardCharsets.UTF_8), 65536);
            bw.write("##fileformat=VCFv4.1\n");
            bw.write("##fileDate="+S.split(" ")[0]+"\n");
            bw.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
            bw.write("##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the reference and alternate alleles in the order listed\">\n");
            bw.write("##FORMAT=<ID=GL,Number=G,Type=Integer,Description=\"Genotype likelihoods for 0/0, 0/1, 1/1, or  0/0, 0/1, 0/2, 1/1, 1/2, 2/2 if 2 alt alleles\">\n");
            bw.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n");
            bw.write("##INFO=<ID=NZ,Number=1,Type=Integer,Description=\"Number of taxa with called genotypes\">\n");
            bw.write("##INFO=<ID=AD,Number=.,Type=Integer,Description=\"Total allelelic depths in order listed starting with REF\">\n");
            bw.write("##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Numbers of allele occurence across taxa in order listed\">\n");
            bw.write("##INFO=<ID=IS,Number=.,Type=Integer,Description=\"Indel sequence of ALT alleles in order listed\">\n");
            bw.write("##INFO=<ID=GN,Number=.,Type=Integer,Description=\"Number of taxa with genotypes AA,AB,BB or AA,AB,AC,BB,BC,CC if 2 alt alleles\">\n");
            bw.write("##INFO=<ID=HT,Number=1,Type=Integer,Description=\"Number of heterozygotes\">\n");
            bw.write("##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Minor allele frequency\">\n");
            bw.write("##ALT=<ID=D,Description=\"Deletion\">\n");
            bw.write("##ALT=<ID=I,Description=\"Insertion\">\n");
            Dyad<int[][], int[]> d = FastCall3.getBins(this.regionStart, this.regionEnd, FastCall3.scanBinSize);
            int[][] binBound = d.getFirstElement();
            int[] binStarts = d.getSecondElement();
            StringBuilder sb = new StringBuilder("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT");
            for (int i = 0; i < taxaNames.length; i++) {
                sb.append("\t").append(taxaNames[i]);
            }
            bw.write(sb.toString());
            bw.newLine();

            String indiCountFolderS = new File(outputDirS, subDirS[1]).getAbsolutePath();
            int nTaxa = taxaNames.length;
            int readThreads = Math.max(1, Math.min(this.threadsNum, nTaxa));
            int offset = 0;
            for (int i = 0; i < binStarts.length; i++) {
                sb.setLength(0);
                sb.append(chrom).append("_").append(binBound[i][0]).append("_").append(binBound[i][1]).append(".iac.gz");
                String iacName = sb.toString();
                streams = new IndividualCountF3.SiteStream[nTaxa];
                int positionNum = -1;
                int nMissing = 0;
                for (int t = 0; t < nTaxa; t++) {
                    String fileS = new File(new File(indiCountFolderS, taxaNames[t]), iacName).getAbsolutePath();
                    IndividualCountF3.SiteStream st = IndividualCountF3.SiteStream.open(fileS);
                    if (st == null) {
                        System.err.println("Warning: "+ fileS +" does not exist");
                        nMissing++;
                        continue;
                    }
                    if (positionNum < 0) {
                        positionNum = st.positionNum;
                    }
                    else if (st.positionNum != positionNum) {
                        throw new IllegalStateException("Site-count mismatch in " + fileS + ": "
                                + st.positionNum + " vs " + positionNum);
                    }
                    streams[t] = st;
                }
                if (positionNum < 0) {
                    throw new IllegalStateException("No individual count files found for " + iacName);
                }
                if (offset + positionNum > positions.length) {
                    throw new ArrayIndexOutOfBoundsException("iac sites " + positionNum
                            + " at offset " + offset + " exceed library positions " + positions.length);
                }

                ExecutorService readPool = Executors.newFixedThreadPool(readThreads);
                int remaining = positionNum;
                int localOffset = 0;
                while (remaining > 0) {
                    final int w = Math.min(window, remaining);
                    final short[][][] windowCounts = new short[w][nTaxa][];
                    final IndividualCountF3.SiteStream[] binStreams = streams;
                    List<Future<?>> reads = new ArrayList<>(nTaxa);
                    for (int t = 0; t < nTaxa; t++) {
                        final int taxon = t;
                        reads.add(readPool.submit(() -> {
                            IndividualCountF3.SiteStream st = binStreams[taxon];
                            for (int s = 0; s < w; s++) {
                                windowCounts[s][taxon] = (st == null) ? null : st.readSite();
                            }
                            return null;
                        }));
                    }
                    for (int r = 0; r < reads.size(); r++) {
                        reads.get(r).get();
                    }

                    final int binOffset = offset + localOffset;
                    final String[] vcfRecords = new String[w];
                    List<Integer> indexList = PArrayUtils.getIndexList(w);
                    indexList.parallelStream().forEach(index -> {
                        int posIdx = binOffset + index;
                        StringBuilder vsb = new StringBuilder();
                        int currentPosition = positions[posIdx];
                        vsb.append(chrom).append("\t").append(currentPosition).append("\t").append(chrom).append("-").append(currentPosition)
                                .append("\t").append(posRefMap.get(currentPosition)).append("\t");
                        int[] altAlleles = posAllelePackMap.get(currentPosition);
                        for (int j = 0; j < altAlleles.length; j++) {
                            vsb.append(AllelePackageF3.getAlleleBase(altAlleles[j])).append(",");
                        }
                        vsb.deleteCharAt(vsb.length()-1).append("\t.\t.\t");
                        vsb.append(this.getInfoAndGenotypes(Arrays.asList(windowCounts[index]), altAlleles));
                        vcfRecords[index] = vsb.toString();
                    });
                    for (int j = 0; j < w; j++) {
                        bw.write(vcfRecords[j]);
                        bw.newLine();
                    }
                    remaining -= w;
                    localOffset += w;
                }
                readPool.shutdown();
                readPool.awaitTermination(Long.MAX_VALUE, TimeUnit.SECONDS);
                closeSiteStreams(streams);
                streams = null;
                offset += positionNum;
                if (nMissing > 0) {
                    System.err.println("Missing individual count files in this bin:\t" + nMissing);
                }
                System.err.println(iacName.split("\\.")[0]+ " genotyping is finished.");
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            closeSiteStreams(streams);
            System.exit(1);
        }
        System.err.println("Final VCF written to stdout");
        System.err.println("Genotyping is finished.");
    }

    private static void closeSiteStreams(IndividualCountF3.SiteStream[] streams) {
        if (streams == null) {
            return;
        }
        for (int i = 0; i < streams.length; i++) {
            if (streams[i] == null) {
                continue;
            }
            try {
                streams[i].close();
            }
            catch (Exception ignore) {
            }
            streams[i] = null;
        }
    }

    /**
     * A Callable implementation for reading and processing individual count data for a single taxon.
     * This class is used by the thread pool in {@link #scanIndiCountsByThreadPool()} to parallelize
     * the reading of individual count files across multiple threads.
     *
     * <p>Each instance processes a single taxon's count file, which contains allele count information
     * for a specific genomic region. The file is expected to be in a binary format that can be
     * read by the {@link IndividualCountF3} class.
     *
     * <p>If the specified file does not exist, this class will log a warning and return null.
     * The calling code is responsible for handling null returns appropriately.
     *
     * @see IndividualCountF3
     * @see #scanIndiCountsByThreadPool()
     * @since 3.0
     */
    class TaxonCountRead implements Callable<IndividualCountF3> {
        /** The path to the individual count file to be processed */
        private final String fileS;

        /**
         * Creates a new TaxonCountRead for the specified file.
         *
         * @param fileS the path to the individual count file to be processed
         */
        public TaxonCountRead(String fileS) {
            this.fileS = fileS;
        }

        @Override
        public IndividualCountF3 call() throws Exception {
            File f = new File (fileS);
            if (!f.exists()) {
                System.err.println("Warning: "+ f.getAbsolutePath()+" does not exist");
                return null;
            }
            IndividualCountF3 inc = new IndividualCountF3(this.fileS);
            return inc;
        }
    }

        /**
         * Scan the individual count files for all taxa in the specified region in parallel across
         * multiple threads.
         *
         * <p>This method uses the {@link java.util.concurrent.ExecutorService} and
         * {@link java.util.concurrent.Future} classes to parallelize the reading of individual
         * count files across multiple threads. The total number of threads is controlled by the
         * {@link #threadsNum} field.
         *
         * <p>The method takes no arguments and returns no value. Instead, it populates the
         * {@link #posRefMap}, {@link #posAllelePackMap}, and {@link #positions} fields with the
         * reference base at each position, the allele pack at each position, and the positions
         * themselves, respectively.
         */
    public void scanIndiCountsByThreadPool () {
        FastaRecordBit frb = genomeFa.getFastaRecordBit(chromIndex);
        posRefMap = new HashMap<>();
        posAllelePackMap = new HashMap<>(vlEndIndex-vlStartIndex);
        positions = new int[vlEndIndex-vlStartIndex];
        for (int i = vlStartIndex; i < vlEndIndex; i++) {
            posRefMap.put(vl.positions[i], String.valueOf(frb.getBase(vl.positions[i]-1)));
            posAllelePackMap.put(vl.positions[i], vl.getAllelePacks(i));
            positions[i-vlStartIndex] = vl.positions[i];
        }
        Set<String> taxaSet = taxaBamsMap.keySet();
        ArrayList<String> taxaList = new ArrayList(taxaSet);
        Collections.sort(taxaList);
        Dyad<int[][], int[]> d = FastCall3.getBins(this.regionStart, this.regionEnd, FastCall3.scanBinSize);
        int[][] binBound = d.getFirstElement();
        int[] binStarts = d.getSecondElement();
        LongAdder counter = new LongAdder();
        ExecutorService pool = Executors.newFixedThreadPool(this.threadsNum);
        List<Future<IndiCount>> resultList = new ArrayList<>();
        for (int i = 0; i < taxaList.size(); i++) {
            List<String> bamPaths = taxaBamsMap.get(taxaList.get(i));
            StringBuilder sb = new StringBuilder(samtoolsPath);
            sb.append(" mpileup --no-output-ends ").append(this.baqMode).append("-q ").append(this.mappingQThresh).append(" -Q ").append(this.baseQThresh).append(" -f ").append(this.referenceFileS);
            for (int j = 0; j < bamPaths.size(); j++) {
                sb.append(" ").append(bamPaths.get(j));
            }
            sb.append(" -l ").append(vLibPosFileS).append(" -r ");
            sb.append(chrom).append(":").append(this.regionStart).append("-").append(this.regionEnd);
            String command = sb.toString();
//            System.err.println(command);
            IndiCount idv = new IndiCount(command, taxaList.get(i), binBound, binStarts, bamPaths, counter);
            Future<IndiCount> result = pool.submit(idv);
            resultList.add(result);
        }
        try {
            pool.shutdown();
            pool.awaitTermination(Long.MAX_VALUE, TimeUnit.SECONDS);
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    /**
     * A Callable implementation for processing individual sample data in parallel.
     * This class is responsible for reading BAM files, counting alleles, and writing
     * the results to individual count files in a binary format.
     *
     * <p>Each instance processes a single taxon's data across one or more genomic bins,
     * with each bin being processed as a separate task in the thread pool. The class
     * implements {@link Callable} to support parallel execution.</p>
     *
     * <p>Key responsibilities include:
     * <ul>
     *   <li>Executing samtools commands to process BAM files</li>
     *   <li>Counting alleles at each variant position</li>
     *   <li>Writing results to binary count files</li>
     *   <li>Handling missing or empty data</li>
     *   <li>Managing file I/O for output files</li>
     * </ul>
     *
     * <p>The output files are written in a binary format that can be efficiently
     * read by {@link IndividualCountF3} for downstream analysis.</p>
     *
     * @see Callable
     * @see IndividualCountF3
     * @see #scanIndiCountsByThreadPool()
     * @since 3.0
     */
    class IndiCount implements Callable<IndiCount> {
        /** The command to execute for processing BAM files */
        String command = null;
        
        /** The name of the taxon being processed */
        String taxonName = null;
        
        /** The output directory for this taxon's count files */
        String indiTaxonDirS = null;
        
        /** Array of bin boundaries [start, end] for genomic regions */
        int[][] binBound = null;
        
        /** Array of bin start positions for quick lookup */
        int[] binStarts = null;
        
        /** List of BAM file paths for this taxon */
        List<String> bamPaths = null;
        
        /** Counter for tracking progress */
        LongAdder counter = null;
        
        /** Output stream for writing binary data */
        DataOutputStream dos = null;
        
        /** Index of the current bin being processed */
        int currentBinIndex = Integer.MIN_VALUE;

        public IndiCount (String command, String taxonName, int[][] binBound, int[] binStarts, List<String> bamPaths, LongAdder counter) {
            this.command = command;
            this.taxonName = taxonName;
            this.binBound = binBound;
            this.binStarts = binStarts;
            this.bamPaths = bamPaths;
            this.counter = counter;
            String indiCountFolderS = new File(outputDirS, subDirS[1]).getAbsolutePath();
            indiTaxonDirS = new File(indiCountFolderS, taxonName).getAbsolutePath();
            new File (indiTaxonDirS).mkdir();
        }

        /**
         * Closes the DataOutputStream and flushes it. This is used to mark
         * the end of a bin in the output file.
         */
        public void closeDos () {
            try {
                dos.flush();
                dos.close();
            }
            catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        }

        /**
         * Opens a new DataOutputStream for writing binary data if the current
         * bin index is different from the previous one.
         *
         * @param queryPos The position of the query site.
         */
        public void setDos (int queryPos) {
            int binIndex = Arrays.binarySearch(binStarts, queryPos);
            if (binIndex < 0) binIndex = -binIndex-2;
            if (binIndex != currentBinIndex) {
                if (currentBinIndex > -1) this.closeDos();
                StringBuilder sb = new StringBuilder();
                sb.append(chrom).append("_").append(binBound[binIndex][0]).append("_").append(binBound[binIndex][1]).append(".iac.gz");
                String outfileS = new File (indiTaxonDirS, sb.toString()).getAbsolutePath();
                dos = IOUtils.getBinaryGzipWriter(outfileS);
                try {
                    dos.writeUTF(this.taxonName);
                    dos.writeShort((short)chrom);
                    dos.writeInt(binBound[binIndex][0]);
                    dos.writeInt(binBound[binIndex][1]);

                    int startIdx = vl.getStartIndex(binBound[binIndex][0]);
                    int endIdx = vl.getEndIndex(binBound[binIndex][1] - 1);
                    int sites = (startIdx == Integer.MIN_VALUE || endIdx == Integer.MIN_VALUE) ? 0 : Math.max(0, endIdx - startIdx);
                    dos.writeInt(sites);
//                    vlBinStartIndex = vl.getStartIndex(binBound[binIndex][0]);
//                    vlBinEndIndex = vl.getEndIndex(binBound[binIndex][1]-1);
//                    dos.writeInt(vlBinEndIndex-vlBinStartIndex);
                }
                catch (Exception e) {
                    e.printStackTrace();
                    System.exit(1);
                }
                currentBinIndex = binIndex;
            }
        }

        /**
         * Writes the allele counts at the current site to the output file.
         *
         * @param alleleCounts The allele counts, where the length of the array is
         * the number of alleles at the site, and each element is the number of
         * times the allele is observed.
         */
        public void writeAlleleCounts (int[] alleleCounts) {
            try {
                dos.writeByte((byte)alleleCounts.length);
                for (int i = 0; i < alleleCounts.length; i++) {
                    if (alleleCounts[i] < Short.MAX_VALUE) {
                        dos.writeShort((short)alleleCounts[i]);
                    }
                    else {
                        dos.writeShort(Short.MAX_VALUE);
                    }
                }
            }
            catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        }

        /**
         * Writes a missing indicator (-1) to the output file.
         *
         * This is used to indicate that the current site is missing data.
         *
         * @see #writeAlleleCounts(int[])
         */
        public void writeMissing () {
            try {
                dos.writeByte((byte)-1);
            }
            catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        }

        /**
         * Writes empty files for each bin in the genomic region.
         *
         * If a file for a bin already exists, it is skipped.
         *
         * For each file, the following information is written:
         * <ol>
         * <li>The taxon name</li>
         * <li>The chromosome number</li>
         * <li>The start and end positions of the bin</li>
         * <li>The number of sites in the bin</li>
         * <li>For each site in the bin, a missing indicator (-1)</li>
         * </ol>
         *
         * @see #writeMissing()
         */
        public void writeEmptyFiles() {
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < binStarts.length; i++) {
                sb.setLength(0);
                sb.append(chrom).append("_").append(binBound[i][0]).append("_").append(binBound[i][1]).append(".iac.gz");
                String outfileS = new File (indiTaxonDirS, sb.toString()).getAbsolutePath();
                File outf = new File (outfileS);
                if (outf.exists()) continue;
                try {
                    dos = IOUtils.getBinaryGzipWriter(outfileS);
                    dos.writeUTF(this.taxonName);
                    dos.writeShort((short)chrom);
                    dos.writeInt(binBound[i][0]);
                    dos.writeInt(binBound[i][1]);

                    int startIdx = vl.getStartIndex(binBound[i][0]);
                    int endIdx = vl.getEndIndex(binBound[i][1] - 1);
                    int len = (startIdx == Integer.MIN_VALUE || endIdx == Integer.MIN_VALUE) ? 0 : Math.max(0, endIdx - startIdx);
//                    vlBinStartIndex = vl.getStartIndex(binBound[i][0]);
//                    vlBinEndIndex = vl.getEndIndex(binBound[i][1]-1);
//                    int len = vlBinEndIndex-vlBinStartIndex;
//                    if (vlBinStartIndex == Integer.MIN_VALUE || vlBinEndIndex == Integer.MIN_VALUE) {
//                        len = 0;
//                    }

                    dos.writeInt(len);
                    for (int j = 0; j < len; j++) {
                        this.writeMissing();
                    }
                    this.closeDos();
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }

        /**
         * Calls mpileup program to calculate the allele counts for each site in this individual.
         *
         * <p>For each site, it writes a missing indicator (-1) if the site is not in the pileup file,
         * or the allele counts if the site is in the pileup file.
         *
         * @throws Exception
         * @return this object
         */
        @Override
        public IndiCount call() throws Exception {
            Runtime rt = Runtime.getRuntime();
            Process p = rt.exec(command);

            try {
                p.getOutputStream().close();
            } catch (Exception e) {}

            Thread errDrainer = new Thread(() -> {
                try (BufferedReader ebr = new BufferedReader(new InputStreamReader(p.getErrorStream()))) {
                    String errLine;
                    while ((errLine = ebr.readLine()) != null) {
                        if (errLine.startsWith("[mpileup]") && errLine.contains("samples in") && errLine.contains("input files")) {
                            continue;
                        }
                        System.err.println("[samtools std err]: " + errLine);
                    }
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            }, "samtools-stderr-drainer");
            errDrainer.setDaemon(true);
            errDrainer.start();

            try {
                String temp = null;
                BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
                DataOutputStream dis = null;
                String current = br.readLine();
                List<String> currentList = null;
                int currentPosition = -1;
                if (current != null) {
                    currentList = PStringUtils.fastSplit(current);
                    currentPosition = Integer.parseInt(currentList.get(1));
                }
                StringBuilder baseS = new StringBuilder();
                StringBuilder indelSb = new StringBuilder();
                for (int i = 0; i < positions.length; i++) {
                    this.setDos(positions[i]);
                    //at the end of pileup file
                    if (current == null) {
                        this.writeMissing();
                    }
                    else {
                        if (positions[i] == currentPosition) {
                            String ref = posRefMap.get(currentPosition);
                            int[] altAlleles = posAllelePackMap.get(currentPosition);
                            baseS.setLength(0);
                            int siteDepth = 0;
                            for (int j = 0; j < bamPaths.size(); j++) {
                                siteDepth+=Integer.parseInt(currentList.get(3+j*3));
                                baseS.append(currentList.get(4+j*3));
                            }
                            if (siteDepth == 0) {
                                this.writeMissing();
                            }
                            else {
//                                FastCall3.removeFirstPositionSign(baseS);
                                int[] alleleCounts = getAlleleCounts (altAlleles, baseS.toString().toUpperCase(), siteDepth, indelSb);
                                this.writeAlleleCounts(alleleCounts);
                            }
                            current = br.readLine();
                            if (current != null) {
                                currentList = PStringUtils.fastSplit(current);
                                currentPosition = Integer.parseInt(currentList.get(1));
                            }
                        }
                        else if (positions[i] < currentPosition) {
                            this.writeMissing();
                        }
                        else {
                            System.err.println("Current position is greater than pileup position. It should not happen. Program quits");
                            System.exit(1);
                        }
                    }
                }
                this.closeDos();
                br.close();

                boolean finished = p.waitFor(5, TimeUnit.HOURS);    // set 5 hours for one sample
                if (!finished) {
                    p.destroyForcibly();
                    throw new RuntimeException("Command timed out: " + command);
                }

                this.writeEmptyFiles();
//                StringBuilder sb = new StringBuilder();
//                sb.append("mpileup command: ").append(command).append("\n");
//                sb.append("mpileup error profile: ");
//                while ((current = bre.readLine()) != null) {
//                    sb.append(current).append("\n");
//                }
//                sb.append("Individual allele counting is completed for taxon ").append(this.taxonName).append("\n");
//                System.err.println(sb.toString());
                errDrainer.join();

                int exit = p.exitValue();
                if (exit != 0) {
                    throw new RuntimeException("Command failed (exit=" + exit + "): " + command);
                }
                System.err.println("Individual allele counting is completed for taxon " + this.taxonName);
            }
            catch (Exception ee) {
                ee.printStackTrace();
            }
            counter.increment();
            int cnt = counter.intValue();
            if (cnt%50 == 0) System.err.println("Finished individual genotype allele counting in " + String.valueOf(cnt) + " taxa. Total: " + String.valueOf(taxaBamsMap.size()));
            return this;
        }
    }

    /**
     * Delete the temporary files.
     *
     * @see #outputDirS
     * @see #subDirS
     * @see #vLibPosFileS
     */
    private void deleteTemperateFile () {
        File f1 = new File(outputDirS, subDirS[0]);
        File f2 = new File(outputDirS, subDirS[1]);
        try {
            FileUtils.cleanDirectory(new File(outputDirS, subDirS[0]));
            FileUtils.cleanDirectory(new File(outputDirS, subDirS[1]));
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        f1.delete();
        f2.delete();
        new File (vLibPosFileS).delete();
    }

    /**
     *
     * @param siteCountList a list of short arrays, each is a SNP site's count of alleles.
     * @param altAlleles an array of alternate alleles.
     * @return a string of info and genotypes.
     */
    private String getInfoAndGenotypes (List<short[]> siteCountList, int[] altAlleles) {
        StringBuilder genoSB = new StringBuilder("GT:AD:GL");
        StringBuilder infoSB = new StringBuilder();
        int dp = 0;
        int nz = 0;
        int alleleNumber = altAlleles.length+1;
        int[] adCnt = new int[alleleNumber];
        int[] acCnt = new int[alleleNumber];
        int[][] gnCnt = new int[alleleNumber][alleleNumber];
        int ht = 0;
        for (int i = 0; i < siteCountList.size(); i++) {
            if (siteCountList.get(i) == null) {
                nz++;
                genoSB.append("\t./.");
                continue;
            }
            if (this.filterDepthAtGenotype) {
                int siteDepth = 0;
                for (int j = 0; j < siteCountList.get(i).length; j++) {
                    siteDepth += siteCountList.get(i)[j];
                }
                double coverage = Double.NaN;
                if (this.taxaNames != null && i < this.taxaNames.length && this.taxaCoverageMap != null) {
                    Double cov = this.taxaCoverageMap.get(this.taxaNames[i]);
                    if (cov != null) coverage = cov;
                }
                if (!this.passDepthCoverageFilter(siteDepth, coverage)) {
                    nz++;
                    genoSB.append("\t./.");
                    continue;
                }
            }
            for (int j = 0; j < alleleNumber; j++) {
                int currentCount = siteCountList.get(i)[j];
                dp+=currentCount;
                adCnt[j] += currentCount;
            }
            String genoS = this.getGenotypeByShort(siteCountList.get(i));
            genoSB.append("\t").append(genoS);
            int g1 = genoS.charAt(0)-'0';
            int g2 = genoS.charAt(2)-'0';
            acCnt[g1]++;
            acCnt[g2]++;
            gnCnt[g1][g2]++;
            if (g1 != g2) ht++;
        }
        nz = siteCountList.size() - nz;
        int sum = 0;
        for (int i = 0; i < acCnt.length; i++) {
            sum+=acCnt[i];
        }
        float maf = (float)((double)acCnt[0]/sum);
        if (maf>0.5) maf = (float)(1-maf);
        infoSB.append("DP=").append(dp).append(";NZ=").append(nz).append(";AD=");
        for (int i = 0; i < adCnt.length; i++) {
            infoSB.append(adCnt[i]).append(",");
        }
        infoSB.deleteCharAt(infoSB.length()-1);
        infoSB.append(";AC=");
        for (int i = 0; i < acCnt.length; i++) {
            infoSB.append(acCnt[i]).append(",");
        }
        infoSB.deleteCharAt(infoSB.length()-1);
        infoSB.append(";GN=");
        for (int i = 0; i < gnCnt.length; i++) {
            for (int j = i; j < gnCnt.length; j++) {
                infoSB.append(gnCnt[i][j]).append(",");
            }
        }
        infoSB.deleteCharAt(infoSB.length()-1);
        infoSB.append(";HT=").append(ht).append(";MAF=").append(maf);
        infoSB.append("\t").append(genoSB);
        return infoSB.toString();
    }

    /**
     * Count the occurrence of each allele in the given pileup string.
     * @param altAlleles the alternative alleles to be counted
     * @param baseS the pileup string
     * @param siteDepth the total depth at the site
     * @param indelSb a StringBuilder to store the indel string
     * @return an array of length (altAlleles.length+1) where the first element is the count of the reference allele
     * and the remaining elements are the counts of the alternative alleles in the same order as in altAlleles
     */
    private int[] getAlleleCounts (int[] altAlleles, String baseS, int siteDepth, StringBuilder indelSb) {
        byte[] baseB = baseS.getBytes();
        int[] altAlleleCounts = new int[altAlleles.length];
        int index = Integer.MIN_VALUE;
        int vCnt = 0;
        for (int i = 0; i < baseB.length; i++) {
            byte queryAlleleCoding = FastCall3.pileupAscIIToAlleleCodingMap.get(baseB[i]);
            int queryIndelLength = 0;
            index = Arrays.binarySearch(AlleleEncoder.alleleCodings, queryAlleleCoding);
            if (index > 3) {
                int startIndex = i+1;
                int endIndex = i+2;
                for (int j = i+2; j < baseB.length; j++) {
                    if (baseB[j] > 57) {
                        endIndex = j;
                        break;
                    }
                }
                indelSb.setLength(0);
                for (int j = startIndex; j < endIndex; j++) {
                    indelSb.append((char)baseB[j]);
                }
                queryIndelLength = Integer.parseInt(indelSb.toString());
                // Skip length digits and indel sequence. The loop's i++ then
                // moves to the next read. An extra +1 ate that next read, and
                // when BAM pileup strings were concatenated it ate the next BAM's first base.
                i+=indelSb.length();
                i+=queryIndelLength;
            }
            for (int j = 0; j < altAlleles.length; j++) {
                if (AllelePackageF3.getAlleleCoding(altAlleles[j]) == queryAlleleCoding && AllelePackageF3.getIndelLength(altAlleles[j]) == queryIndelLength) {
                    altAlleleCounts[j]++;
                    vCnt++;
                }
            }
        }
        int[] alleleCounts = new int[altAlleles.length+1];
        alleleCounts[0] = siteDepth - vCnt;
        for (int i = 0; i < altAlleles.length; i++) {
            alleleCounts[i+1] = altAlleleCounts[i];
        }
        return alleleCounts;
    }

        /**
         * Calculate the genotype of a given site from its allele counts.
         * @param cnt the allele counts of a site.
         * @return the genotype of the site in the format of "a1/a2:cnt1,cnt2,...,cntN:likelihood1,likelihood2,...,likelihoodN",
         *         where a1 and a2 are the two alleles of the site, cnt1,cnt2,...,cntN are the counts of alleles, and likelihood1,likelihood2,...,likelihoodN are
         *         the likelihoods of the corresponding genotypes.
         */
    private String getGenotypeByShort (short[] cnt) {
        //in case some allele depth is greater than maxFactorial, to keep the original allele counts
        short[] oriCnt = null;
        int n = cnt.length*(cnt.length+1)/2;
        int[] likelihood = new int[n];
        int sum = 0;
        for (int i = 0; i < cnt.length; i++) sum+=cnt[i];
        if (sum == 0) return "./.";
        else if (sum > this.maxFactorial) {
            oriCnt = new short[cnt.length];
            System.arraycopy(cnt, 0, oriCnt, 0, cnt.length);
            double portion = (double)this.maxFactorial/sum;
            for (int i = 0; i < cnt.length; i++) {
                cnt[i] = (short)(cnt[i]*portion);
            }
            sum = this.maxFactorial;
        }
        double coe = this.factorialMap.get(sum);
        for (int i = 0; i < cnt.length; i++) coe = coe/this.factorialMap.get(cnt[i]);
        double max = Double.MAX_VALUE;
        int a1 = 0;
        int a2 = 0;
        for (int i = 0; i < cnt.length; i++) {
            for (int j = i; j < cnt.length; j++) {
                int index = (j*(j+1)/2)+i;
                double value = Double.MAX_VALUE;
                if (i == j) {
                    value = -Math.log10(coe*Math.pow((1-0.75*this.combinedErrorRate), cnt[i])*Math.pow(this.combinedErrorRate /4, (sum-cnt[i])));
                }
                else {
                    value = -Math.log10(coe*Math.pow((0.5-this.combinedErrorRate /4), cnt[i]+cnt[j])*Math.pow(this.combinedErrorRate /4, (sum-cnt[i]-cnt[j])));
                }
                if (value < max) {
                    max = value;
                    a1 = i;
                    a2 = j;
                }
                likelihood[index] = (int)Math.round(value);
            }
        }
        StringBuilder sb = new StringBuilder();
        sb.append(a1).append("/").append(a2).append(":");
        if (sum > this.maxFactorial) {
            for (int i = 0; i < oriCnt.length; i++) sb.append(oriCnt[i]).append(",");
        }
        else {
            for (int i = 0; i < cnt.length; i++) sb.append(cnt[i]).append(",");
        }
        sb.deleteCharAt(sb.length()-1); sb.append(":");
        for (int i = 0; i < likelihood.length; i++) sb.append(likelihood[i]).append(",");
        sb.deleteCharAt(sb.length()-1);
        return sb.toString();
    }

    /**
     * This method will process the variation library according to the given chromosome and region.
     * It will output a file which contains the positions of the variation library.
     * The file name is "{chrom}_{regionStart}_{regionEnd}.pos.txt", where {chrom} is the chromosome number,
     * {regionStart} and {regionEnd} are the start and end positions of the region respectively.
     * The file contains one position per line, with the format of "chrom\tPosition".
     * If the chromosome number of the library does not match the given one, or the region is incorrectly set,
     * the program will quit.
     */
    private void processVariationLibrary () {
        StringBuilder sb = new StringBuilder();
        sb.append(this.chrom).append("_").append(this.regionStart).append("_").append(regionEnd).append(".pos.txt");
        this.vLibPosFileS = new File (this.outputDirS, sb.toString()).getAbsolutePath();
        this.vl = new VariationLibraryF3(this.libFileS);
        if (this.chrom != vl.getChrom()) {
            System.err.println("The chromosome number of library and the specified one do not match. Program quits.");
            System.exit(0);
        }
        try {
            vlStartIndex = vl.getStartIndex(this.regionStart);
            vlEndIndex = vl.getEndIndex(this.regionEnd);
            if (vlStartIndex == Integer.MIN_VALUE || vlEndIndex == Integer.MIN_VALUE) {
                System.err.println("The chromosome region was incorrectly set. Program quits.");
                System.exit(0);
            }
            BufferedWriter bw = IOUtils.getTextWriter(this.vLibPosFileS);
            for (int i = vlStartIndex; i < vlEndIndex; i++) {
                sb.setLength(0);
                sb.append(this.chrom).append("\t").append(vl.getPosition(i));
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * Initialize a map to store the factorials from 0 to maxFactorial.
     * The key is the index of the factorial, and the value is the corresponding factorial.
     * For example, if maxFactorial = 5, then the map will contain the following entries:
     * <p>
     * 0 -> 1
     * 1 -> 1
     * 2 -> 2
     * 3 -> 6
     * 4 -> 24
     * 5 -> 120
     * <p>
     * This map is used to calculate the binomial coefficients.
     */
    private void creatFactorialMap () {
        this.factorialMap = HashIntDoubleMaps.getDefaultFactory().newMutableMap();
        for (int i = 0; i < this.maxFactorial+1; i++) {
            this.factorialMap.put(i, factorial(i));
        }
    }

    /**
     * This method will create a directory with the given name,
     * and create several subdirectories with the names specified in
     * the subDirS array.
     * <p>
     * The subdirectories are:
     * <p>
     * mpileup
     * indiVCF
     * VCF
     * <p>
     * The output directory name is specified by the outputDirS variable.
     */
    public void mkDir () {
        File f = new File (this.outputDirS);
        f.mkdir();
        for (int i = 0; i < subDirS.length; i++) {
            f = new File(outputDirS, subDirS[i]);
            f.mkdir();
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
        this.taxaBamsMap = new HashMap<>();
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
                List<String> bamList = Arrays.asList(bams);
                taxaBamsMap.put(tem[0], bamList);
                taxaCoverageMap.put(tem[0], Double.valueOf(tem[1]));
                nBam+=bams.length;
            }
            taxaNames = taxaList.toArray(new String[taxaList.size()]);
            Arrays.sort(taxaNames);
            HashSet<String> taxaSet = new HashSet<>(taxaList);
            if (taxaSet.size() != taxaNames.length) {
                System.err.println("Taxa names are not unique. Programs quits");
                System.exit(0);
            }
            System.err.println("Created TaxaBamMap from " + taxaBamMapFileS);
            System.err.println("Taxa number:\t"+String.valueOf(taxaNames.length));
            System.err.println("Bam file number in TaxaBamMap:\t"+String.valueOf(nBam));
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
