package pgl.app.fastCall;

import gnu.trove.list.array.TByteArrayList;
import gnu.trove.set.hash.TByteHashSet;
import gnu.trove.set.hash.TIntHashSet;
import org.apache.commons.cli.*;
import org.apache.commons.math3.stat.inference.ChiSquareTest;
import pgl.AppAbstract;
import pgl.PGLAPPEntrance;
import pgl.infra.dna.FastaBit;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PArrayUtils;
import pgl.infra.utils.PStringUtils;
import java.io.*;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.Map.Entry;
import java.util.concurrent.*;
import java.util.concurrent.atomic.LongAdder;

import static cern.jet.math.Arithmetic.factorial;

public class FastCall extends AppAbstract {
    String referenceFileS =null;
    String taxaBamMapFileS = null;
    String vcfDirS = null;
    String samtoolsPath = null;
    int currentChr = Integer.MIN_VALUE;
    int regionStart = Integer.MIN_VALUE;
    int regionEnd = Integer.MIN_VALUE;
    int[] chroms = null;
    int[] chromLength = null;
    String[] taxaNames = null;
    double[] taxaCoverage = null;
    String[] bamPaths = null;
    int threadsNum = 32;
    int[][] bounds = null;
    HashMap<String, String> bamPathTaxamap = null;
    int[] bamTaxaIndices = null;
    HashMap<String, String[]> taxaBamPathMap = null;
    HashMap<String, Double> taxaCoverageMap = null;
    HashMap<String, String> bamPathPileupPathMap = null;
    ConcurrentHashMap<Integer, Double> factorialMap = new ConcurrentHashMap();
    int maxFactorial = 150;
    FastaBit genomeFa = null;

    String[] baseS = {"A", "C", "G", "T"};
    //A, C, G, T
    byte[] bases = {65, 67, 71, 84};
    //D, I
    byte[] possibleIndel = {68, 73};
    //A, C, D, G, I, T
    byte[] possibleAllele = {65, 67, 68, 71, 73, 84};


    int binSize = 100000;
    double individualDepthRatioThresh = 0.4;
    double individualThirdAlleleRatioThresh = 0.2;
    //Mendelian segregation test
    double segregationPValueThresh = 1;
    //Both sequencing error and alignment error
    double combinedErrorRate = 0.05;
    //Minimum number of taxa where the minor allele shows up
    int minorOccurrenceThresh = 2;

    public FastCall(String parameterFileS) {
        this.callSNP(parameterFileS);
    }

    public FastCall (String[] args) {
        this.creatAppOptions();
        this.retrieveAppParameters(args);
        this.callSNP();
    }

    @Override
    public void creatAppOptions() {
        options.addOption("app", true, "App name.");
        options.addOption("a", true, "Reference genome file with an index file (.fai). The reference should be in Fasta format. " +
            "Chromosomes are labelled as numbers (1,2,3,4,5...).");
        options.addOption("b", true, "Taxa bam information file, including the info about what bams are included for each taxon");
        options.addOption("c", true, "Combined error rate from sequencing and alignment. When the error rate is low, " +
            "heterogeneous sites are more likely to be called as heterozygous, and vice versa. It is 0.05 by default.");
        options.addOption("d", true, "Individual depth ratio. This is the depth of the  lower-depth allele vs. total depth. " +
            "When the threshold is low, heterogeneous sites are more likely to be called as heterozygous, and vice versa. It is 0.4 by default.");
        options.addOption("e", true, "P-value threshold of Mendelian segregation test. Lower threshold (e.g. 1, the test will be removed) " +
            "is recommended if rare alleles are the major interest. It is 1 by default. ");
        options.addOption("f", true, "Minor allele occurrence threshold, representing the minimum number of taxa where the minor allele exist. " +
            "It is 2 by default.");
        options.addOption("g", true, "Chromosome or region on which genotyping will be performed (e.g. chromosome 1 is designated as 1. " +
                "Region 1bp to 100000bp on chromosome 1 is 1:1,100000)");
        options.addOption("h", true, "VCF output directory");
        options.addOption("i", true, "Number of threads for pileup");
        options.addOption("j", true, "The path of samtools");
    }

    @Override
    public void retrieveAppParameters(String[] args) {
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine line = parser.parse(options, args);
            this.referenceFileS = line.getOptionValue("a");
            this.taxaBamMapFileS = line.getOptionValue("b");
            this.combinedErrorRate = Double.parseDouble(line.getOptionValue("c"));
            this.individualDepthRatioThresh = Double.parseDouble(line.getOptionValue("d"));
            this.segregationPValueThresh = Double.parseDouble(line.getOptionValue("e"));
            this.minorOccurrenceThresh = Integer.parseInt(line.getOptionValue("f"));
            if (line.getOptionValue("g").contains(":")) {
                String[] temp = line.getOptionValue("g").split(":");
                currentChr = Integer.parseInt(temp[0]);
                temp = temp[1].split(",");
                regionStart = Integer.parseInt(temp[0]);
                regionEnd = Integer.parseInt(temp[1]);
            }
            else {
                currentChr = Integer.valueOf(line.getOptionValue("g"));
            }
            this.vcfDirS = line.getOptionValue("h");
            this.threadsNum = Integer.valueOf(line.getOptionValue("i"));
            this.samtoolsPath = line.getOptionValue("j");

        }
        catch(Exception e) {
            e.printStackTrace();
            System.out.println("\nThere are input errors in the command line. Program stops.");
            this.printInstructionAndUsage();
            System.exit(0);
        }
    }
    public void printInstructionAndUsage () {
        System.out.println(PGLAPPEntrance.getTIGERIntroduction());
        System.out.println("Below are the commands of FastCall.");
        this.printUsage();
    }
    private void callSNP() {
        long start = System.nanoTime();
        System.out.println("Reading reference genome from "+ referenceFileS);
        genomeFa = new FastaBit(referenceFileS);
        System.out.println("Reading reference genome took " + String.format("%.2f", Benchmark.getTimeSpanSeconds(start)) + "s");
        genomeFa.sortByDescriptionValue();
        chroms = new int[genomeFa.getSeqNumber()];
        chromLength = new int[genomeFa.getSeqNumber()];
        for (int i = 0; i < genomeFa.getSeqNumber(); i++) {
            chroms[i] = Integer.valueOf(genomeFa.getDescription(i));
            chromLength[i] = genomeFa.getSeqLength(i);
        }
        String pileupDirS = new File(new File(vcfDirS).getParent(), "pileup").getAbsolutePath();
        new File(pileupDirS).mkdir();
        new File(vcfDirS).mkdir();
        this.getTaxaBamMap(taxaBamMapFileS);
        this.creatPileupMap(pileupDirS);
        this.creatFactorialMap();
        this.callSNPByRegion(currentChr, regionStart, regionEnd, referenceFileS, vcfDirS);
        File[] fs = new File(pileupDirS).listFiles();
        for (int i = 0; i < fs.length; i++) fs[i].delete();
        System.out.println("Variant calling completed");
    }

    public void callSNP (String parameterFileS) {
        ArrayList<String> pLineList = new ArrayList();
        try {
            BufferedReader br = IOUtils.getTextReader(parameterFileS);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("@")) continue;
                if (temp.startsWith("#")) continue;
                if (temp.isEmpty()) continue;
                pLineList.add(temp);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        String referenceFileS = pLineList.get(0);
        String taxaBamMapFileS = pLineList.get(1);
        this.combinedErrorRate = Double.parseDouble(pLineList.get(2));
        this.segregationPValueThresh = Double.parseDouble(pLineList.get(3));
        this.individualDepthRatioThresh = Double.parseDouble(pLineList.get(4));
        this.minorOccurrenceThresh = Integer.parseInt(pLineList.get(5));
        int currentChr = Integer.MIN_VALUE;
        int regionStart = Integer.MIN_VALUE;
        int regionEnd = Integer.MIN_VALUE;
        if (pLineList.get(6).contains(":")) {
            String[] temp = pLineList.get(6).split(":");
            currentChr = Integer.valueOf(temp[0]);
            temp = temp[1].split(",");
            regionStart = Integer.valueOf(temp[0]);
            regionEnd = Integer.valueOf(temp[1]);
        }
        else {
            currentChr = Integer.valueOf(pLineList.get(6));
        }
        String vcfDirS = pLineList.get(7);
        this.threadsNum = Integer.valueOf(pLineList.get(8));
        this.samtoolsPath = pLineList.get(9);
        long start = System.nanoTime();
        System.out.println("Reading reference genome from "+ referenceFileS);
        genomeFa = new FastaBit(referenceFileS);
        System.out.println("Reading reference genome took " + String.format("%.2f", Benchmark.getTimeSpanSeconds(start)) + "s");
        genomeFa.sortByDescriptionValue();
        chroms = new int[genomeFa.getSeqNumber()];
        chromLength = new int[genomeFa.getSeqNumber()];
        for (int i = 0; i < genomeFa.getSeqNumber(); i++) {
            chroms[i] = Integer.valueOf(genomeFa.getDescription(i));
            chromLength[i] = genomeFa.getSeqLength(i);
        }
        String pileupDirS = new File(new File(vcfDirS).getParent(), "pileup").getAbsolutePath();
        new File(pileupDirS).mkdir();
        new File(vcfDirS).mkdir();
        this.getTaxaBamMap(taxaBamMapFileS);
        this.creatPileupMap(pileupDirS);
        this.creatFactorialMap();
        this.callSNPByRegion(currentChr, regionStart, regionEnd, referenceFileS, vcfDirS);
        File[] fs = new File(pileupDirS).listFiles();
        for (int i = 0; i < fs.length; i++) fs[i].delete();
        System.out.println("Variant calling completed");
    }

    private void callSNPByRegion(int currentChr, int rStart, int rEnd, String referenceFileS, String vcfDirS) {
        int chrIndex = Arrays.binarySearch(chroms, currentChr);
        String chrSeq = genomeFa.getSeq(chrIndex).toUpperCase();
        int regionStart = Integer.MIN_VALUE;
        int regionEnd = Integer.MIN_VALUE;
        if (rStart == Integer.MIN_VALUE) {
            regionStart = 1;
            regionEnd = chrSeq.length();
        }
        else {
            if (rEnd > chrSeq.length() || rStart < 1) {
                System.out.println("Chromosome/region setting issues. Program quits");
                System.exit(0);
            }
            else if (rEnd <= rStart) {
                System.out.println("Chromosome/region setting issues. Program quits");
                System.exit(0);
            }
            else {
                regionStart = rStart;
                regionEnd = rEnd;
            }
        }
        this.performPileup(currentChr, regionStart, regionEnd, referenceFileS);
        String outfileS = "chr"+PStringUtils.getNDigitNumber(3, currentChr)+".vcf";
        outfileS = new File (vcfDirS, outfileS).getAbsolutePath();
        //note both start and end are inclusive
        int[][] binBound = this.creatBins(currentChr, binSize, regionStart, regionEnd);
        try {
            BufferedReader[] pileupReaders = this.getPileupReaders();
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write(this.getAnnotation(referenceFileS));
            bw.write(this.getVCFHeader());
            bw.newLine();
            List<Integer> bamIndexList = new ArrayList();
            for (int i = 0; i < this.bamPaths.length; i++) {
                bamIndexList.add(i);
            }
            List<Integer> positionIndexList = new ArrayList();
            for (int i = 0; i < binSize; i++) {
                positionIndexList.add(i);
            }
            StringBuilder[][] baseSb = this.getPopulatedBaseBuilder(1, this.binSize);
            int[][] depth = getPopulatedDepthArray(1, this.binSize);
            List<List<String>>[] pileupResult = null;
            for (int i = 0; i < binBound.length; i++) {
                long startTimePoint = System.nanoTime();
                int binStart = binBound[i][0];
                int binEnd = binBound[i][1];
                List<String>[] remainder = new List[bamPaths.length];
                for (int j = 0; j < remainder.length; j++) remainder[j] = new ArrayList<String>();
                pileupResult = this.getPileupResult(currentChr, binStart, binEnd, pileupReaders, remainder, bamIndexList);
                if (i < binBound.length - 1) {
                    this.resetBaseBuilder(baseSb);
                    this.resetDepthArray(depth);
                }
                else {
                    baseSb = this.getPopulatedBaseBuilder(binStart, binEnd);
                    depth = this.getPopulatedDepthArray(binStart, binEnd);
                }
                this.fillDepthAndBase(pileupResult, baseSb, depth, binStart);
                String[][] base = this.getBaseMatrix(baseSb);
                ArrayList<Integer> positionList = this.getPositionList(binStart, binEnd);
                String[] subVCFs = new String[positionList.size()];
                this.calculateVCF(subVCFs, positionList, currentChr, binStart, chrSeq, depth, base);
                for (int j = 0; j < positionList.size(); j++) {
                    String vcfStr = subVCFs[j];
                    if (vcfStr == null) continue;
                    bw.write(vcfStr);
                    bw.newLine();
                }
                StringBuilder sb = new StringBuilder();
                sb.append("Bin from ").append(binStart).append(" to ").append(binEnd).append(" is finished. Took ").append(Benchmark.getTimeSpanSeconds(startTimePoint)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
                System.out.println(sb.toString());
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("Chromosome " + String.valueOf(currentChr) + " is finished. File written to " + outfileS + "\n");
    }

    private void calculateVCF(String[] subVCFs, List<Integer> positionList, int currentChr, int startPos, String chrSeq, int[][] depth, String[][] base) {
        positionList.parallelStream().forEach(position -> {
            int index = position-startPos;
            byte refBase = (byte)(chrSeq.charAt(position-1));
            int baseIndex = Arrays.binarySearch(bases, refBase);
            if (baseIndex < 0) {

            }
            else {
                String vcfStr = this.getVCFStringV2(base[index], depth[index], currentChr, position, refBase); //limit on sequencing depth
                if (vcfStr != null) {
                    subVCFs[index] = vcfStr;
                }
            }
        });
    }

    /**
     * For mixed depth sequencing samples
     * @param base
     * @param depth
     * @param currentChr
     * @param position
     * @param refBase
     * @return
     */
    private String getVCFStringV2 (String[] base, int[] depth, int currentChr, int position, byte refBase) {
        int totalDepth = 0;
        for (int i = 0; i < depth.length; i++) {
            totalDepth+=depth[i];
        }
        if (totalDepth == 0) return null;
        TByteArrayList bList = new TByteArrayList();
        byte[] ba = null;
        boolean ifRecordedDeletion = false;
        TIntHashSet insertionLengthSet = new TIntHashSet();
        TIntHashSet deletionLengthSet = new TIntHashSet();
        int[][] pAlleleCount = new int[base.length][this.possibleAllele.length];
        int[] refDepth = new int[base.length];
        for (int i = 0; i < base.length; i++) {
            base[i] = base[i].toUpperCase();
            bList.clear(base[i].length());
            ba = base[i].getBytes();
            for (int j = 0; j < ba.length; j++) {
                if (ba[j] == '.') {
                    bList.add(refBase);
                }
                else if (ba[j] == ',') {
                    bList.add(refBase);
                }
                else if (ba[j] == 'A') {
                    bList.add((byte)65);
                }
                else if (ba[j] == 'C') {
                    bList.add((byte)67);
                }
                else if (ba[j] == 'G') {
                    bList.add((byte)71);
                }
                else if (ba[j] == 'T') {
                    bList.add((byte)84);
                }
                else if (ba[j] == '+') {
                    int endIndex = j+2;
                    for (int k = j+1; k < ba.length; k++) {
                        if (ba[k] > 57) {
                            endIndex = k;
                            break;
                        }
                    }
                    StringBuilder sb = new StringBuilder();
                    for (int k = j+1; k < endIndex; k++) {
                        sb.append((char)ba[k]);
                    }
                    int length = Integer.valueOf(sb.toString());
                    insertionLengthSet.add(length);
                    j+=sb.length();
                    j+=length;
                    bList.replace(bList.size()-1, (byte)73);
                }
                else if (ba[j] == '-') {
                    int endIndex = j+2;
                    for (int k = j+1; k < ba.length; k++) {
                        if (ba[k] > 57) {
                            endIndex = k;
                            break;
                        }
                    }
                    StringBuilder sb = new StringBuilder();
                    for (int k = j+1; k < endIndex; k++) {
                        sb.append((char)ba[k]);
                    }
                    int length = Integer.valueOf(sb.toString());
                    deletionLengthSet.add(length);
                    j+=sb.length();
                    j+=length;
                    bList.replace(bList.size()-1, (byte)68);
                }
                else if (ba[j] == '^') {
                    j++;
                }
                else if (ba[j] == '*') {
                    bList.add(refBase);
                    ifRecordedDeletion = true;
                }
                //N, n, $, >, <, |
                else {
                    //do nothing
                }
            }
            byte[] taxonBase = bList.toArray();
            for (int j = 0; j < taxonBase.length; j++) {
                int index = Arrays.binarySearch(this.possibleAllele, taxonBase[j]);
                pAlleleCount[i][index]++;
            }
            int altSum = 0;
            for (int j = 0; j < pAlleleCount[i].length; j++) {
                if (this.possibleAllele[j] == refBase) continue;
//                if (this.possibleAllele[j] == 68) continue;
//                if (this.possibleAllele[j] == 73) continue;
                altSum+=pAlleleCount[i][j];
            }
            refDepth[i] = depth[i] - altSum;
        }

        int[] indelTypeCount = new int[2];
        indelTypeCount[0] = deletionLengthSet.size();
        indelTypeCount[1] = insertionLengthSet.size();

        //****************************Filter1 IndelFilter*****************************************
        //Too many indels usually means mis-alignment
        if (indelTypeCount[0]+indelTypeCount[1] > 1) return null;
        //=======================================Filter1=====================================================

        //****************************Filter2 Depth_ratio_test*****************************************
        //In any individual, alt allele show up < 2 times, ignore
        //In any individual, depth > 0.2 * coverage and depth < 2*coverage. Pick up
        //In any individual, alt allele show up > individualDepthRatioThresh. Pick up
        //When depth is low, tend to have assembly errors, LTR

        TByteHashSet altAlleleSet = new TByteHashSet();
        for (int i = 0; i < pAlleleCount[0].length; i++) {
            if (possibleAllele[i] == refBase) continue;
            for (int j = 0; j < pAlleleCount.length; j++) {
                if (depth[j] < 2) {}
                else {
                    double r = depth[j]/this.taxaCoverage[j];
                    if (r < 0.2) continue;
                    if (r > 3) continue;
                    if ((double)pAlleleCount[j][i]/depth[j] > this.individualDepthRatioThresh) {
                        altAlleleSet.add(this.possibleAllele[i]);
                    }
                }
            }
        }
        byte[] altAllele = altAlleleSet.toArray();
        if (altAllele.length == 0) return null;
        Arrays.sort(altAllele);
        //=======================================Filter2=====================================================


        int[] altAllele2PAlleleIndex = new int[altAllele.length];
        for (int i = 0; i < this.possibleAllele.length; i++) {
            int index = Arrays.binarySearch(altAllele, this.possibleAllele[i]);
            if (index < 0) continue;
            altAllele2PAlleleIndex[index] = i;
        }
        int[] altAlleleTotalDepth = new int[altAllele.length];
        for (int i = 0; i < altAllele.length; i++) {
            for (int j = 0; j < base.length; j++) {
                altAlleleTotalDepth[i]+=pAlleleCount[j][altAllele2PAlleleIndex[i]];
            }
        }
        int[] altAlleleDepthDesendingIndex = PArrayUtils.getIndicesByDescendingValue(altAlleleTotalDepth);
        int refTotalDepth = 0;
        for (int i = 0; i < refDepth.length; i++) refTotalDepth+=refDepth[i];

        //****************************Filter3 third_allele_test************************************************
        //individual should not have the third allele
        if (altAllele.length > 1) {
            for (int i = 0; i < base.length; i++) {
                int[] tempCnt = new int[altAllele.length];
                for (int j = 0; j < altAllele.length; j++) {
                    tempCnt[j] = pAlleleCount[i][altAllele2PAlleleIndex[j]];
                }
                int sum = refDepth[i];
                double[] v = new double[altAllele.length+1];
                for (int j = 0; j < tempCnt.length; j++) v[j] = (double)tempCnt[j]/sum;
                v[v.length-1] = (double)refDepth[i]/sum;
                Arrays.sort(v);
                if (v[v.length-3] > individualThirdAlleleRatioThresh) return null;
            }
        }
        //===========================Filter3=========================================================

        //****************************Filter4 Segregation_test*****************************************
        long[] observed = null;
        double[] expected = null;
        double[] segregationP = null;
        if (this.segregationPValueThresh < 1) {
            observed = new long[base.length];
            expected = new double[base.length];
            segregationP = new double[altAllele.length];
            ChiSquareTest ct = new ChiSquareTest();
            int cnt = 0;
            for (int i = 0; i < altAllele.length; i++) {
                double r = (double) altAlleleTotalDepth[i] / (refTotalDepth + altAlleleTotalDepth[i]);
                for (int j = 0; j < base.length; j++) {
                    observed[j] = pAlleleCount[j][altAllele2PAlleleIndex[i]];
                    expected[j] = r;
                }
                segregationP[i] = ct.chiSquareTest(expected, observed);
                if (segregationP[i] > this.segregationPValueThresh) cnt++;
            }
            if (cnt == altAllele.length) return null;
        }
        //===========================Filter4=========================================================

        int nonMissingCnt = 0;
        int[] refAndAllelePresence = new int[altAllele.length+1];
        for (int i = 0; i < base.length; i++) {
            if (refDepth[i] != 0) {
                nonMissingCnt++;
                refAndAllelePresence[0]++;
            }
            else {
                for (int j = 0; j < altAllele.length; j++) {
                    if (pAlleleCount[i][altAllele2PAlleleIndex[j]] != 0) {
                        nonMissingCnt++;
                        break;
                    }
                }
            }
            for (int j = 0; j < altAllele.length; j++) {
                if (pAlleleCount[i][altAllele2PAlleleIndex[j]] != 0) {
                    refAndAllelePresence[j+1]++;
                }
            }
        }
        StringBuilder sb = new StringBuilder();
        sb.append(currentChr).append("\t").append(position).append("\t").append(currentChr).append("-").append(position).append("\t").append((char)refBase).append("\t");
        for (int i = 0; i < altAllele.length; i++) sb.append(String.valueOf((char)altAllele[altAlleleDepthDesendingIndex[i]])).append(",");
        sb.deleteCharAt(sb.length()-1);
        sb.append("\t.\t.\t").append("DP=").append(totalDepth).append(";AD=").append(refTotalDepth);
        for (int i = 0; i < altAllele.length; i++) sb.append(",").append(altAlleleTotalDepth[altAlleleDepthDesendingIndex[i]]);
        sb.append(";NZ=").append(nonMissingCnt).append(";AP=");
        for (int i = 0; i < refAndAllelePresence.length; i++) sb.append(refAndAllelePresence[i]).append(",");

        //****************************Filter5 Minor allele occurence *****************************************
        Arrays.sort(refAndAllelePresence);
        if (refAndAllelePresence.length < 3) {
            if (refAndAllelePresence[0] < this.minorOccurrenceThresh) return null;
        }
        else {
            if (refAndAllelePresence[refAndAllelePresence.length-2] < this.minorOccurrenceThresh) return null;
        }
        //===========================Filter5=========================================================


        sb.deleteCharAt(sb.length()-1);
        sb.append(";PV=");
        if (this.segregationPValueThresh < 1) {
            for (int i = 0; i < altAllele.length; i++) sb.append(String.format("%.2e", segregationP[i])).append(",");
        }
        else {
            sb.append("NA,");
        }
        sb.deleteCharAt(sb.length()-1);
        sb.append(";DI=");
        for (int i = 0; i < indelTypeCount.length; i++) sb.append(indelTypeCount[i]).append(",");
        sb.deleteCharAt(sb.length()-1);
        sb.append("\t").append("GT:AD:PL");

        for (int i = 0; i < base.length; i++) {
            int[] dep = new int[altAllele.length+1];
            dep[0] = refDepth[i];
            for (int j = 0; j < altAllele.length; j++) {
                dep[j+1] = pAlleleCount[i][altAllele2PAlleleIndex[altAlleleDepthDesendingIndex[j]]];
            }
            sb.append("\t").append(this.getGenotype(dep));
        }
        return sb.toString();
    }

    private ArrayList<Integer> getPositionList (int startPos, int endPos) {
        ArrayList<Integer> positionList = new ArrayList();
        for (int i = startPos; i <= endPos; i++) {
            positionList.add(i);
        }
        return positionList;
    }

    private String[][] getBaseMatrix (StringBuilder[][] baseSb) {
        String[][] base = new String[baseSb.length][baseSb[0].length];
        for (int i = 0; i < base.length; i++) {
            for (int j = 0; j < base[0].length; j++) base[i][j] = baseSb[i][j].toString();
        }
        return base;
    }

    private void fillDepthAndBase(List<List<String>>[] pileupResult, StringBuilder[][] baseSb, int[][] depth, int startPos) {
        Set<Entry<String, String[]>> entries = this.taxaBamPathMap.entrySet();
        List<Entry<String,String[]>> entryList = new ArrayList(entries);
        entryList.parallelStream().forEach(e -> {
            String taxa = e.getKey();
            int taxaIndex = Arrays.binarySearch(this.taxaNames, taxa);
            String[] bams = e.getValue();
            int[] indices = new int[bams.length];
            for (int i = 0; i < indices.length; i++) {
                indices[i] = Arrays.binarySearch(this.bamPaths, bams[i]);
            }

            int count = 0;
            String b = null;
            try {
                for (int i = 0; i < bams.length; i++) {
                    List<List<String>> lines = pileupResult[indices[i]];
                    count = lines.size();
                    b = bams[i];
                    for (int j = 0; j < lines.size(); j++) {
                        List<String> split = lines.get(j);
                        //if (split.get(2).startsWith("N") || split.get(2).startsWith("n")) continue;
                        String refB = split.get(2);
                        if (Arrays.binarySearch(baseS, refB) < 0) continue;
                        int siteIndex = Integer.valueOf(split.get(1)) - startPos;
                        depth[siteIndex][taxaIndex]+=Integer.valueOf(split.get(3));
                        baseSb[siteIndex][taxaIndex].append(split.get(4));
                    }
                }
            }
            catch (Exception ee) {
                System.out.println(b);
                System.out.println(count);
                ee.printStackTrace();
                System.exit(1);
            }
        });
    }

    private int[][] getPopulatedDepthArray (int startPos, int endPos) {
        int[][] depth = new int[endPos-startPos+1][this.taxaNames.length];
        return depth;
    }

    private void resetDepthArray (int[][] depth) {
        Arrays.stream(depth).forEach(x -> Arrays.fill(x, 0));
    }

    private StringBuilder[][] getPopulatedBaseBuilder(int startPos, int endPos) {
        StringBuilder[][] sbs = new StringBuilder[endPos-startPos+1][this.taxaNames.length];
        for (int i = 0; i < sbs.length; i++) {
            for (int j = 0; j < sbs[0].length; j++) sbs[i][j] = new StringBuilder();
        }
        return sbs;
    }

    private void resetBaseBuilder (StringBuilder[][] sbs) {
        for (int i = 0; i < sbs.length; i++) {
            for (int j = 0; j < sbs[0].length; j++) sbs[i][j].setLength(0);
        }
    }

    private List<List<String>>[] getPileupResult (int currentChr, int binStart, int binEnd, BufferedReader[] pileupReaders, List<String>[] remainders, List<Integer> bamIndexList) {
        ArrayList<String> empty = new ArrayList();
        List<List<String>>[] pileupResult = new List[bamPaths.length];
        bamIndexList.parallelStream().forEach(bamIndex -> {
            ArrayList<List<String>> lineList = new ArrayList();
            BufferedReader br = pileupReaders[bamIndex];
            List<String> remainder = remainders[bamIndex];
            boolean flag = false;
            if (remainder.size() == 0) {
                String temp = null;
                try {
                    temp = br.readLine();
                }
                catch (Exception e) {}
                if (temp != null) {
                    List<String> split = PStringUtils.fastSplit(temp, "\t");
                    int currentPos = Integer.valueOf(split.get(1));
                    if (currentPos > binEnd) {
                        remainders[bamIndex] = split;
                    }
                    else {
                        lineList.add(split);
                        flag = true;
                    }
                }
            }
            else {
                int currentPos = Integer.valueOf(remainder.get(1));
                if (currentPos <= binEnd) {
                    lineList.add(remainder);
                    flag = true;
                    remainders[bamIndex] = empty;
                }
            }
            if (flag == true) {
                try {
                    String temp;
                    while ((temp = br.readLine()) != null) {
                        List<String> split = PStringUtils.fastSplit(temp, "\t");
                        int currentPos = Integer.valueOf(split.get(1));
                        if (currentPos < binEnd) {
                            lineList.add(split);
                        }
                        else if (currentPos == binEnd){
                            lineList.add(split);
                            remainders[bamIndex] = empty;
                            break;
                        }
                        else {
                            remainders[bamIndex] = split;
                            break;
                        }
                    }
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
            }
            pileupResult[bamIndex] = lineList;
        });
        return pileupResult;
    }

    private BufferedReader[] getPileupReaders () {
        BufferedReader[] brs = new BufferedReader[this.bamPaths.length];
        try {
            for (int i = 0; i < this.bamPaths.length; i++) {
                String pileupFileS = this.bamPathPileupPathMap.get(bamPaths[i]);
                File pileupF = new File(pileupFileS);
                if (!pileupF.exists()) {
                    brs[i] = new BufferedReader(new StringReader(""), 65536);
                    continue;
                }
                brs[i] = new BufferedReader (new FileReader(pileupF), 65536);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return brs;
    }

    private HashMap<String, BufferedReader> getBamPathPileupReaderMap () {
        HashMap<String, BufferedReader> bamPathPileupReaderMap = new HashMap();
        try {
            for (int i = 0; i < this.bamPaths.length; i++) {
                String pileupFileS = this.bamPathPileupPathMap.get(bamPaths[i]);
                File pileupF = new File(pileupFileS);
                if (!pileupF.exists()) {
                    BufferedReader br = new BufferedReader(new StringReader(""), 1024);
                    bamPathPileupReaderMap.put(bamPaths[i], br);
                    continue;
                }
                BufferedReader br = new BufferedReader (new FileReader(pileupF), 1024);
                bamPathPileupReaderMap.put(bamPaths[i], br);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        return bamPathPileupReaderMap;
    }

    class RunPileup implements Callable <RunPileup> {
        String command = null;
        LongAdder counter;

        public RunPileup (String cmd, LongAdder counter) {
            this.command = cmd;
            this.counter = counter;
        }

        @Override
        public RunPileup call() throws Exception {
            try {
                Runtime rt = Runtime.getRuntime();
                Process p = rt.exec(command);
                                BufferedReader br = new BufferedReader(new InputStreamReader(p.getErrorStream()));
                                String temp = null;
                                while ((temp = br.readLine()) != null) {
                                    System.out.println(command);
                                    System.out.println(temp);
                                }
                p.waitFor();
                counter.increment();
                int count = counter.intValue();
                if (count%50 == 0) {
                    System.out.println(String.valueOf(count)+ " files pileuped.");
                }
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            return null;
        }
    }


    private void performPileup (int currentChr, int startPos, int endPos, String referenceFileS) {
        System.out.println("Pileup is being performed on chromosome "+String.valueOf(currentChr)+" from "+String.valueOf(startPos)+" to "+String.valueOf(endPos));
        long timeStart = System.nanoTime();
        List<String> bamList = Arrays.asList(bamPaths);
        try {
            LongAdder counter = new LongAdder();
            ExecutorService pool = Executors.newFixedThreadPool(this.threadsNum);
            for (int i = 0; i < bamList.size(); i++) {
                String bamFileS = bamList.get(i);
                String pileupFileS = this.bamPathPileupPathMap.get(bamFileS);
                StringBuilder sb = new StringBuilder(this.samtoolsPath);
                sb.append(" mpileup -A -B -q 30 -Q 20 -f ").append(referenceFileS).append(" ").append(bamFileS).append(" -r ");
                sb.append(currentChr).append(":").append(startPos).append("-").append(endPos).append(" -o ").append(pileupFileS);
                String command = sb.toString();
                RunPileup rp = new RunPileup(command, counter);
                Future<RunPileup> f = pool.submit(rp);
            }
            pool.shutdown();
            pool.awaitTermination(Long.MAX_VALUE, TimeUnit.SECONDS);
        }
        catch (Exception e) {
            e.printStackTrace();
        }

        System.out.println("Pileup is finished. Time took " + String.format("%.2f", Benchmark.getTimeSpanMinutes(timeStart)) + " mins");
    }

    private int[][] creatBins (int currentChr, int binSize, int regionStart, int regionEnd) {
        int[][] binBound = PArrayUtils.getSubsetsIndicesBySubsetSize(regionEnd-regionStart+1, binSize);
        for (int i = 0; i < binBound.length; i++) {
            binBound[i][0]+=regionStart;
            binBound[i][1]+=regionStart;
            binBound[i][1]--;
        }
        System.out.println("SNP calling will performed on chromosome "+String.valueOf(currentChr)+" from "+String.valueOf(regionStart)+" to "+String.valueOf(regionEnd));
        System.out.println("Chromosome " + String.valueOf(currentChr) +"  is devided into " + String.valueOf(binBound.length) + " bins. Bin size: " + String.valueOf(this.binSize) + " bp");
        return binBound;
    }

    private String getVCFHeader () {
        StringBuilder sb = new StringBuilder();
        sb.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
        for (int i = 0; i < taxaNames.length; i++) {
            sb.append("\t").append(taxaNames[i]);
        }
        return sb.toString();
    }


    private void creatFactorialMap () {
        for (int i = 0; i < this.maxFactorial+1; i++) {
            this.factorialMap.put(i, factorial(i));
        }
    }

    public String getGenotype (int[] cnt) {
        int n = cnt.length*(cnt.length+1)/2;
        int[] likelihood = new int[n];
        int sum = 0;
        for (int i = 0; i < cnt.length; i++) sum+=cnt[i];
        if (sum == 0) return "./.";
        else if (sum > this.maxFactorial) {
            double portion = (double)this.maxFactorial/sum;
            for (int i = 0; i < cnt.length; i++) {
                cnt[i] = (int)(cnt[i]*portion);
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
        for (int i = 0; i < cnt.length; i++) sb.append(cnt[i]).append(",");
        sb.deleteCharAt(sb.length()-1); sb.append(":");
        for (int i = 0; i < likelihood.length; i++) sb.append(likelihood[i]).append(",");
        sb.deleteCharAt(sb.length()-1);
        return sb.toString();
    }

    private void creatPileupMap (String pileupDirS) {
        bamPathPileupPathMap = new HashMap();
        Set<Entry<String, String[]>> entries = taxaBamPathMap.entrySet();
        for (Entry<String, String[]> e : entries) {
            String[] bams = e.getValue();
            for (String bam : bams) {
                String pileupFileS = new File (pileupDirS, new File(bam).getName().replaceFirst(".bam", ".pileup.txt")).getAbsolutePath();
                bamPathPileupPathMap.put(bam, pileupFileS);
            }
        }
    }

    private void getTaxaBamMap (String taxaBamMapFileS) {
        this.bamPathTaxamap = new HashMap();
        this.taxaBamPathMap = new HashMap();
        this.taxaCoverageMap = new HashMap();
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
                    this.bamPathTaxamap.put(bams[i], tem[0]);
                }
                Arrays.sort(bams);
                taxaBamPathMap.put(tem[0], bams);
                taxaCoverageMap.put(tem[0], Double.valueOf(tem[1]));
                nBam+=bams.length;
            }
            taxaNames = taxaList.toArray(new String[taxaList.size()]);
            Arrays.sort(taxaNames);
            HashSet<String> taxaSet = new HashSet(taxaList);
            if (taxaSet.size() != taxaNames.length) {
                System.out.println("Taxa names are not unique. Programs quits");
                System.exit(0);
            }
            this.taxaCoverage = new double[taxaNames.length];
            for (int i = 0; i < taxaCoverage.length; i++) {
                taxaCoverage[i] = this.taxaCoverageMap.get(taxaNames[i]);
            }
            this.bamPaths = pathList.toArray(new String[pathList.size()]);
            Arrays.sort(bamPaths);
            this.bamTaxaIndices = new int[bamPaths.length];
            for (int i = 0; i < bamPaths.length; i++) {
                int index = Arrays.binarySearch(taxaNames, this.bamPathTaxamap.get(bamPaths[i]));
                bamTaxaIndices[i] = index;
            }
            bounds = PArrayUtils.getSubsetsIndicesBySubsetSize(taxaNames.length, this.threadsNum);
            System.out.println("Created TaxaBamMap from" + taxaBamMapFileS);
            System.out.println("Taxa number:\t"+String.valueOf(taxaNames.length));
            System.out.println("Bam file number in TaxaBamMap:\t"+String.valueOf(nBam));
            System.out.println();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    private String getAnnotation (String referenceFileS) {
        StringBuilder sb = new StringBuilder();
        sb.append("##fileformat=VCFv4.1\n");
        SimpleDateFormat sdf = new SimpleDateFormat("MM/dd/yyyy HH:mm:ss.SSS");
        Date dt = new Date();
        String S = sdf.format(dt);
        sb.append("##fileDate=").append(S.split(" ")[0]).append("\n");
        sb.append("##reference=").append(referenceFileS).append("\n");
        sb.append("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"").append("Total depth").append("\">\n");
        sb.append("##INFO=<ID=AD,Number=2+,Type=Integer,Description=\"").append("Total allelelic depths in order listed starting with REF").append("\">\n");
        sb.append("##INFO=<ID=NZ,Number=1,Type=Integer,Description=\"").append("Number of individuals with alleles present").append("\">\n");
        sb.append("##INFO=<ID=AP,Number=2+,Type=Integer,Description=\"").append("Number of individuals in which an allele is present").append("\">\n");
        sb.append("##INFO=<ID=PV,Number=1+,Type=Float,Description=\"").append("Segreagation test P-Value of alternative alleles").append("\">\n");
        sb.append("##INFO=<ID=DI,Number=2,Type=Integer,Description=\"").append("Number of deletion and insertion type").append("\">\n");
        sb.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"").append("Genotype").append("\">\n");
        sb.append("##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"").append("Allelic depths for the reference and alternate alleles in the order listed").append("\">\n");
        sb.append("##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"").append("Genotype likelihoods for 0/0, 0/1, 1/1, or 0/0, 0/1, 0/2, 1/1, 1/2, 2/2 if 2 alt alleles").append("\">\n");
        return sb.toString();
    }

}
