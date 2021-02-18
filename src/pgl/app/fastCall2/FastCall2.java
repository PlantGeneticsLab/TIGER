package pgl.app.fastCall2;

import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import pgl.AppUtils;
import pgl.PGLConstraints;
import pgl.infra.dna.FastaBit;
import pgl.infra.utils.*;

import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.File;
import java.io.InputStreamReader;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.LongAdder;

public class FastCall2 {
    //Reference genome file with an index file (.fai). The reference should be in Fasta format. Chromosomes are labled as 1-based numbers (1,2,3,4,5...).
    String referenceFileS = null;
    //The taxaRefBam file containing information of taxon and its corresponding bam files. The bam file should have .bai file in the same folder
    String taxaRefBamFileS = null;
    //The path of samtools
    String samtoolsPath = null;
    //VCF output directory
    String outputDirS = null;
    //Minimum mapping quality (MQ) for an alignment to be used for variation calling.
    int mappingQThresh = 30;
    //Minimum base quality (BQ) for a base to be used for variation calling.
    int baseQThresh = 20;
    //Minimum read depth count (MDC) for variation calling, meaning that sites with depth lower than the minimum will not be taken into account for variation discovery.
    int mdcThresh = 2;
    //Minimum read depth ratio (MiDR) for variation calling, meaning that sites with depth lower than the MiDR by the individual coverage will not be considered for variation discovery.
    double mindrThresh = 0.2;
    //Maximum read depth ratio (MaDR) for variation calling, meaning that sites with depth higher than the MaDR by the individual coverage will not be considered for variation discovery.
    double maxdrTrresh = 3;
    //Homozygous ratio (HoR) for variation calling, meaning that the depth of alternative allele is greater than HoR are considered to homozygous.
    double horThresh = 0.8;
    //Heterozygous ratio (HeR) for variation calling, meaning that the depth of alternative allele is greater than HR and less than (1-HR) are considered to be hets.
    double herThresh = 0.4;
    //Third allele depth ratio (TDR) for variation calling. If the depth of the third allele is greater than TDR by the individual coverage, the site will be ignored. Otherwise, the third allele will be considered as sequencing error.
    double tdrTresh = 0.2;
    //Current chromosome for variation calling
    int chrom = Integer.MIN_VALUE;
    //Starting position of the specified region for variation calling, inclusive
    int regionStart = Integer.MIN_VALUE;
    //Ending position the specified regionfor variation calling, exclusive
    int regionEnd = Integer.MIN_VALUE;
    //Number of threads (taxa number to be processed at the same time)
    int threadsNum = PGLConstraints.parallelLevel;
    //Current step ID of the pipeline
    int step = Integer.MIN_VALUE;
    //Two many indels meaning alignment error
    int indelTypeThresh = 1;

    HashMap<String, String[]> taxaBamPathMap = null;
    HashMap<String, Double> taxaCoverageMap = null;

    String[] taxaNames = null;
    int chromLength = Integer.MIN_VALUE;
    //genome block for individual ing
    int binSize = 5000000;
    //+, -, A, C, G, T,
    public final byte[] pileupAlleleAscIIs = {43, 45, 65, 67, 71, 84};

    public FastCall2 (String parameterFileS) {
        this.parseParameters(parameterFileS);
        if (step == 1) {
            this.variationDiscovery();
        }
        else if (step == 2) {

        }
    }

    private void variationDiscovery () {
        File outDir = new File (this.outputDirS);
        outDir.mkdir();
        File[] taxaOutDirs = new File[this.taxaNames.length];
        for (int i = 0; i < taxaNames.length; i++) {
            File f = new File (outDir, taxaNames[i]);
            f.mkdir();
            taxaOutDirs[i] = f;
        }
        int actualChrLength = this.regionEnd - this.regionStart;
        //starting from actual genome position
        int[][] binBound = PArrayUtils.getSubsetsIndicesBySubsetSize (actualChrLength, this.binSize);
        int[] binStarts = new int[binBound.length];
        for (int i = 0; i < binBound.length; i++) {
            binBound[i][0] = binBound[i][0]+regionStart;
            binBound[i][1] = binBound[i][1]+regionStart;
            binStarts[i] = binBound[i][0];
        }
        try {
            LongAdder counter = new LongAdder();
            ExecutorService pool = Executors.newFixedThreadPool(1);
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < taxaNames.length; i++) {
                String[] bamFiles = this.taxaBamPathMap.get(taxaNames[i]);
                sb.setLength(0);
                sb.append(this.samtoolsPath).append(" mpileup -A -B -q ").append(this.mappingQThresh).append(" -Q ").append(this.baseQThresh).append(" -f ").append(this.referenceFileS);
                for (int j = 0; j < bamFiles.length; j++) {
                    sb.append(" ").append(bamFiles[j]);
                }
                sb.append(" -r ").append(this.chrom).append(":").append(this.regionStart).append("-").append(regionEnd-1);
                String command = sb.toString();
                System.out.println(command);
                TaxonCall tc = new TaxonCall(command, binBound, binStarts, taxaNames[i], taxaOutDirs[i], counter);
                Future<TaxonCall> f = pool.submit(tc);
            }
            pool.shutdown();
            pool.awaitTermination(Long.MAX_VALUE, TimeUnit.SECONDS);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    class TaxonCall implements Callable <TaxonCall> {

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
        int[] alleleCount = new int[pileupAlleleAscIIs.length];
        IntOpenHashSet insertionLengthSet = new IntOpenHashSet();
        IntOpenHashSet deletionLengthSet = new IntOpenHashSet();


        public TaxonCall (String command, int[][] binBound, int[] binStarts, String taxon, File outDir, LongAdder counter) {
            this.command = command;
            this.binBound = binBound;
            this.binStarts = binStarts;
            this.taxon = taxon;
            this.taxonCoverage = taxaCoverageMap.get(taxon);
            this.outDir = outDir;
            this.counter = counter;
        }

        public DataOutputStream getSNPDos (int binIndex) {
            StringBuilder sb = new StringBuilder();
            sb.append(chrom).append("_").append(binBound[binIndex][0]).append("_").append(binBound[binIndex][1]).append(".snp.ing.gz");
            String outfileS = new File (outDir, sb.toString()).getAbsolutePath();
            System.out.println(outfileS);
            DataOutputStream dos = IOUtils.getBinaryGzipWriter(outfileS);
            return dos;
        }

        public DataOutputStream getINDDos (int binIndex) {
            StringBuilder sb = new StringBuilder();
            sb.append(chrom).append("_").append(binBound[binIndex][0]).append("_").append(binBound[binIndex][1]).append(".ind.ing.gz");
            String outfileS = new File (outDir, sb.toString()).getAbsolutePath();
            DataOutputStream dos = IOUtils.getBinaryGzipWriter(outfileS);
            return dos;
        }

        private void initialize1 () {
            this.baseSb.setLength(0);
            this.currentDepth = 0;
        }

        private void initialize2 () {
            Arrays.fill(alleleCount, 0);
            insertionLengthSet.clear();
            deletionLengthSet.clear();
        }

        public void processPileupLine (String line) {
            lList = PStringUtils.fastSplit(line);
            this.initialize1();
            currentPos = Integer.parseInt(lList.get(1));
            for (int i = 3; i < lList.size(); i+=3) {
                currentDepth+=Integer.parseInt(lList.get(i));
                baseSb.append(lList.get(i+1));
            }
            if (currentDepth < mdcThresh) return;
            double siteDepthRatio = (double)currentDepth/this.taxonCoverage;
            if (siteDepthRatio < mindrThresh) return;
            if (siteDepthRatio > maxdrTrresh) return;
            String baseS = baseSb.toString().toUpperCase();
            byte[] baseB = baseS.getBytes();
            this.initialize2();
            int index = 0;
            int vCnt = 0;
            for (int i = 0; i < baseB.length; i++) {
                index = Arrays.binarySearch(pileupAlleleAscIIs, baseB[i]);
                if (index < 0) continue;
                else if (index < 2) {
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
                    if (index == 0) insertionLengthSet.add(length);
                    else deletionLengthSet.add(length);
                    i+=baseSb.length();
                    i+=length;
                }
                alleleCount[index]++;
                vCnt++;
            }
            if (vCnt == 0) return;
            if (insertionLengthSet.size()+deletionLengthSet.size() > indelTypeThresh) return;
            int[] alleleCountDesendingIndex = PArrayUtils.getIndicesByDescendingValue(alleleCount);
            double alleleDepthRatio = (double)alleleCount[alleleCountDesendingIndex[0]]/currentDepth;
            if (alleleDepthRatio < herThresh) return;
            else if (alleleDepthRatio > 1 - herThresh && alleleDepthRatio < horThresh) return;
            if (alleleCount[alleleCountDesendingIndex[1]] != 0) {
                alleleDepthRatio = (double)alleleCount[alleleCountDesendingIndex[1]]/currentDepth;
                if (alleleDepthRatio > tdrTresh) return;
            }

        }

        @Override
        public TaxonCall call() throws Exception {
            try {
                Runtime rt = Runtime.getRuntime();
                Process p = rt.exec(command);
                String temp = null;

                BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
                DataOutputStream snpDos = null;
                DataOutputStream indDos = null;

                while ((temp = br.readLine()) != null) {
                    this.processPileupLine(temp);

                }
                br.close();

                BufferedReader bre = new BufferedReader(new InputStreamReader(p.getErrorStream()));
                while ((temp = bre.readLine()) != null) {
                    if (temp.startsWith("[m")) continue;
                    System.out.println(command);
                    System.out.println(temp);
                }
                bre.close();
                p.waitFor();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            counter.increment();
            int count = counter.intValue();
            if (count%50 == 0) {
                System.out.println("Variation calling has been performed for "+ String.valueOf(count)+ " taxa.");
            }
            return null;
        }
    }



    private void parseParameters (String parameterFileS) {
        Dyad<List<String>, List<String>> d = AppUtils.getParameterList(parameterFileS);
        List<String> pLineList = d.getFirstElement();
        List<String> sLineList = d.getSecondElement();
        this.step = Integer.parseInt(sLineList.get(0).split("\\s+")[1]);
        if (step == 1) {
            this.parseParametersStep1(pLineList);
        }
        else if (step == 2) {

        }
    }

    private void parseParametersStep1 (List<String> pLineList) {
        this.referenceFileS = pLineList.get(0);
        taxaRefBamFileS = pLineList.get(1);
        this.mappingQThresh = Integer.parseInt(pLineList.get(2));
        this.baseQThresh = Integer.parseInt(pLineList.get(3));
        this.mdcThresh = Integer.parseInt(pLineList.get(4));
        this.mindrThresh = Double.parseDouble(pLineList.get(5));
        this.maxdrTrresh = Double.parseDouble(pLineList.get(6));
        this.horThresh = Double.parseDouble(pLineList.get(7));
        this.herThresh = Double.parseDouble(pLineList.get(8));
        this.tdrTresh = Double.parseDouble(pLineList.get(9));
        String[] tem = pLineList.get(10).split(":");
        this.chrom = Integer.parseInt(tem[0]);
        long start = System.nanoTime();
        System.out.println("Reading reference genome from "+ referenceFileS);
        FastaBit genomeFa = new FastaBit(referenceFileS);
        System.out.println("Reading reference genome took " + String.format("%.2f", Benchmark.getTimeSpanSeconds(start)) + "s");
        int chromIndex = genomeFa.getIndexByName(String.valueOf(this.chrom));
        if (tem.length == 1) {
            this.regionStart = 1;
            this.regionEnd = genomeFa.getSeqLength(chromIndex)+1;
        }
        else if (tem.length == 2) {
            tem = tem[1].split(",");
            this.regionStart = Integer.parseInt(tem[0]);
            this.regionEnd = Integer.parseInt(tem[1])+1;
        }
        this.threadsNum = Integer.parseInt(pLineList.get(11));
        this.outputDirS = pLineList.get(12);
        this.samtoolsPath = pLineList.get(13);


        this.parseTaxaBamMap(this.taxaRefBamFileS);
    }

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
            System.out.println("Created TaxaBamMap from" + taxaBamMapFileS);
            System.out.println("Taxa number:\t"+String.valueOf(taxaNames.length));
            System.out.println("Bam file number in TaxaBamMap:\t"+String.valueOf(nBam));
            System.out.println();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
