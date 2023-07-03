package pgl.app.fastCall2;

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
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.Dyad;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import java.io.*;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.LongAdder;

import static cern.jet.math.Arithmetic.factorial;

class ScanGenotype extends AppAbstract {
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
    //combined: sequencing error and alignment error
    double combinedErrorRate = 0.05;
    //The path of samtools
    String samtoolsPath = null;
    //VCF output directory
    String outputDirS = null;
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
    VariationLibrary vl = null;
    int vlStartIndex = Integer.MIN_VALUE;
    int vlEndIndex = Integer.MIN_VALUE;
    HashMap<Integer, String> posRefMap = new HashMap<>();
    HashMap<Integer, short[]> posCodedAlleleMap = new HashMap<>();
    int[] positions = null;
    int vlBinStartIndex = 0;
    int vlBinEndIndex = 0;

    public ScanGenotype (String[] args) {
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

    @Override
    public void creatAppOptions() {
        options.addOption("app", true, "App name.");
        options.addOption("step", true, "Step of FastCall 2 (e.g. 1).");
        options.addOption("a", true, "Reference genome file with an index file (.fai). The reference should be in Fasta format. " +
            "Chromosomes are labelled as numbers (1,2,3,4,5...). It is recommended to use reference chromosome while perform genotyping for " +
            "each chromosome because loading reference genome would be much faster.");
        options.addOption("b", true, "The taxaBamMap file contains information of taxon and its corresponding bam files. " +
            "The bam file should have .bai file in the same folder.");
        options.addOption("c", true, "The genetic variation library file, which is from step 2.");
        options.addOption("d", true, "Chromosome or region on which genotyping will be performed (e.g. chromosome 1 is designated as 1. " +
            "Region 1bp to 100000bp on chromosome 1 is 1:1,100000)");
        options.addOption("e", true, "Combined error rate of sequencing and misalignment. Heterozygous read mapping are more " +
            "likely to be genotyped as homozygote when the combined error rate is high.");
        options.addOption("f", true, "The path of samtools.");
        options.addOption("g", true, "Number of threads. It is 32 by default.");
        options.addOption("h", true, "The directory of VCF output.");
    }

    @Override
    public void retrieveAppParameters(String[] args) {
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine line = parser.parse(options, args);
            this.referenceFileS = line.getOptionValue("a");
            this.taxaBamMapFileS = line.getOptionValue("b");
            this.libFileS = line.getOptionValue("c");
            String[] tem = line.getOptionValue("d").split(":");
            this.chrom = Integer.parseInt(tem[0]);
            long start = System.nanoTime();
            System.out.println("Reading reference genome from "+ referenceFileS);
            this.genomeFa = new FastaBit(referenceFileS);
            System.out.println("Reading reference genome took " + String.format("%.2f", Benchmark.getTimeSpanSeconds(start)) + "s");
            this.chromIndex = genomeFa.getIndexByName(String.valueOf(this.chrom));
            if (tem.length == 1) {
                this.regionStart = 1;
                this.regionEnd = genomeFa.getSeqLength(chromIndex)+1;
            }
            else if (tem.length == 2) {
                tem = tem[1].split(",");
                this.regionStart = Integer.parseInt(tem[0]);
                this.regionEnd = Integer.parseInt(tem[1])+1;
            }
            this.combinedErrorRate = Double.parseDouble(line.getOptionValue("e"));
            this.samtoolsPath = line.getOptionValue("f");
            this.threadsNum = Integer.parseInt(line.getOptionValue("g"));
            this.outputDirS = line.getOptionValue("h");
        }
        catch(Exception e) {
            e.printStackTrace();
            System.out.println("\nThere are input errors in the command line. Program stops.");
            this.printInstructionAndUsage();
            System.exit(0);
        }
        this.parseTaxaBamMap(this.taxaBamMapFileS);
    }

    @Override
    public void printInstructionAndUsage() {
        System.out.println(PGLAPPEntrance.getTIGERIntroduction());
        System.out.println("Below are the commands of Step 3 of FastCall 2.");
        this.printUsage();
    }

    public ScanGenotype (List<String> pLineList) {
        this.parseParameters(pLineList);
        this.mkDir();
        this.processVariationLibrary();
        this.creatFactorialMap();

        /*
        Output by individual allele count, fast
         */
        this.scanIndiCountsByThreadPool();
        this.mkFinalVCFFromIndiCounts();
    }

    public void mkFinalVCFFromIndiCounts () {
        String outfileS = new File(outputDirS, subDirS[2]).getAbsolutePath();
        outfileS = new File(outfileS, "chr"+PStringUtils.getNDigitNumber(3, chrom)+".vcf").getAbsolutePath();
        try {
            SimpleDateFormat sdf = new SimpleDateFormat("MM/dd/yyyy HH:mm:ss.SSS");
            Date dt = new Date();
            String S = sdf.format(dt);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("##fileformat=VCFv4.1\n");
            bw.write("##fileDate="+S.split(" ")[0]+"\n");
            bw.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
            bw.write("##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the reference and alternate alleles in the order listed\">\n");
            bw.write("##FORMAT=<ID=GL,Number=G,Type=Integer,Description=\"Genotype likelihoods for 0/0, 0/1, 1/1, or  0/0, 0/1, 0/2, 1/1, 1/2, 2/2 if 2 alt alleles\">\n");
            bw.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n");
            bw.write("##INFO=<ID=NZ,Number=1,Type=Integer,Description=\"Number of taxa with called genotypes\">\n");
            bw.write("##INFO=<ID=AD,Number=.,Type=Integer,Description=\"Total allelelic depths in order listed starting with REF\">\n");
            bw.write("##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Numbers of allele occurence across taxa in order listed\">\n");
            bw.write("##INFO=<ID=IL,Number=.,Type=Integer,Description=\"Indel length of ALT alleles in order listed\">\n");
            bw.write("##INFO=<ID=GN,Number=.,Type=Integer,Description=\"Number of taxa with genotypes AA,AB,BB or AA,AB,AC,BB,BC,CC if 2 alt alleles\">\n");
            bw.write("##INFO=<ID=HT,Number=1,Type=Integer,Description=\"Number of heterozygotes\">\n");
            bw.write("##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Minor allele frequency\">\n");
            bw.write("##ALT=<ID=D,Description=\"Deletion\">\n");
            bw.write("##ALT=<ID=I,Description=\"Insertion\">\n");
            Dyad<int[][], int[]> d = FastCall2.getBinsScanning(this.regionStart, this.regionEnd);
            int[][] binBound = d.getFirstElement();
            int[] binStarts = d.getSecondElement();
                StringBuilder sb = new StringBuilder("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT");
            for (int i = 0; i < taxaNames.length; i++) {
                sb.append("\t").append(taxaNames[i]);
            }
            bw.write(sb.toString());
            bw.newLine();
            List<Future<IndividualCount>> futureList = new ArrayList<>();
            List<IndividualCount> incList = new ArrayList<>();
            vlBinStartIndex = 0;
            vlBinEndIndex = 0;
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
                        Future<IndividualCount> result = pool.submit(tr);
                        futureList.add(result);
                    }
                    pool.shutdown();
                    pool.awaitTermination(Long.MAX_VALUE, TimeUnit.SECONDS);
                    for (int j = 0; j < futureList.size(); j++) {
                        IndividualCount inc = futureList.get(j).get();
                        if (inc == null) continue;
                        incList.add(inc);
                    }
                    Collections.sort(incList);
                }
                catch (Exception e) {
                    e.printStackTrace();
                    System.exit(1);
                }
                vlBinEndIndex = vlBinStartIndex + incList.get(0).alleleNum.length;
                List<Integer> indexList = new ArrayList<>();
                for (int j = 0; j < incList.get(0).alleleNum.length; j++) {
                    indexList.add(j);
                }
                String[] vcfRecords = new String[indexList.size()];
                indexList.parallelStream().forEach(index -> {
                    StringBuilder vsb = new StringBuilder();
                    int currentPosition = positions[index+vlBinStartIndex];
                    vsb.append(chrom).append("\t").append(currentPosition).append("\t").append(chrom).append("-").append(currentPosition)
                            .append("\t").append(posRefMap.get(currentPosition)).append("\t");
                    short[] codedAlts = posCodedAlleleMap.get(currentPosition);
                    for (int j = 0; j < codedAlts.length; j++) {
                        vsb.append(FastCall2.getAlleleBaseFromAlleleCodingLength(codedAlts[j])).append(",");
                    }
                    vsb.deleteCharAt(vsb.length()-1).append("\t.\t.\t");
                    List<short[]> siteCountsList = new ArrayList<>();
                    for (int j = 0; j < incList.size(); j++) {
                        siteCountsList.add(incList.get(j).alleleCounts[index]);
                    }
                    vsb.append(this.getInfoAndGenotypes(siteCountsList, codedAlts));
                    vcfRecords[index] = vsb.toString();
                });
                for (int j = 0; j < vcfRecords.length; j++) {
                    bw.write(vcfRecords[j]);
                    bw.newLine();
                }
                vlBinStartIndex = vlBinEndIndex;
                System.out.println(sb.toString().split("\\.")[0]+ " genotyping is finished.");
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        this.deleteTemperateFile();
        System.out.println("Final VCF is completed at " + outfileS);
        System.out.println("Step 3 is finished.");
    }

    class TaxonCountRead implements Callable<IndividualCount> {
        String fileS;
        public TaxonCountRead (String fileS) {
            this.fileS = fileS;
        }

        @Override
        public IndividualCount call() throws Exception {
            File f = new File (fileS);
            if (!f.exists()) {
                System.out.println("Warning: "+ f.getAbsolutePath()+" does not exist");
                return null;
            }
            IndividualCount inc = new IndividualCount(this.fileS);
            return inc;
        }
    }

    public void scanIndiCountsByThreadPool () {
        FastaRecordBit frb = genomeFa.getFastaRecordBit(chromIndex);
        posRefMap = new HashMap<>();
        posCodedAlleleMap = new HashMap<>();
        List<String> altList = new ArrayList<>();
        positions = new int[vlEndIndex-vlStartIndex];
        for (int i = vlStartIndex; i < vlEndIndex; i++) {
            posRefMap.put(vl.positions[i], String.valueOf(frb.getBase(vl.positions[i]-1)));
            posCodedAlleleMap.put(vl.positions[i], vl.alleleCodingLengths[i]);
            positions[i-vlStartIndex] = vl.positions[i];
        }
        Set<String> taxaSet = taxaBamsMap.keySet();
        ArrayList<String> taxaList = new ArrayList(taxaSet);
        Collections.sort(taxaList);
        Dyad<int[][], int[]> d = FastCall2.getBinsScanning(this.regionStart, this.regionEnd);
        int[][] binBound = d.getFirstElement();
        int[] binStarts = d.getSecondElement();
        LongAdder counter = new LongAdder();
        ExecutorService pool = Executors.newFixedThreadPool(this.threadsNum);
        List<Future<IndiCount>> resultList = new ArrayList<>();
        for (int i = 0; i < taxaList.size(); i++) {
            List<String> bamPaths = taxaBamsMap.get(taxaList.get(i));
            StringBuilder sb = new StringBuilder(samtoolsPath);
            sb.append(" mpileup -A -B -q 20 -Q 20 -f ").append(this.referenceFileS);
            for (int j = 0; j < bamPaths.size(); j++) {
                sb.append(" ").append(bamPaths.get(j));
            }
            sb.append(" -l ").append(vLibPosFileS).append(" -r ");
            sb.append(chrom).append(":").append(this.regionStart).append("-").append(this.regionEnd);
            String command = sb.toString();
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

    class IndiCount implements Callable<IndiCount> {
        String command = null;
        String taxonName = null;
        String indiTaxonDirS = null;
        int[][] binBound = null;
        int[] binStarts = null;
        List<String> bamPaths = null;
        LongAdder counter = null;
        DataOutputStream dos = null;
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
                    vlBinStartIndex = vl.getStartIndex(binBound[binIndex][0]);
                    vlBinEndIndex = vl.getEndIndex(binBound[binIndex][1]-1);
                    dos.writeInt(vlBinEndIndex-vlBinStartIndex);
                }
                catch (Exception e) {
                    e.printStackTrace();
                    System.exit(1);
                }
                currentBinIndex = binIndex;
            }
        }

        public void writeAlleleCounts (int[] alleleCounts) {
            try {
                dos.writeByte((byte)alleleCounts.length);
                for (int i = 0; i < alleleCounts.length; i++) {
                    dos.writeShort((short)alleleCounts[i]);
                }
            }
            catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        }

        public void writeMissing () {
            try {
                dos.writeByte((byte)-1);
            }
            catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        }

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
                    vlBinStartIndex = vl.getStartIndex(binBound[i][0]);
                    vlBinEndIndex = vl.getEndIndex(binBound[i][1]-1);
                    int len = vlBinEndIndex-vlBinStartIndex;
                    if (vlBinStartIndex == Integer.MIN_VALUE || vlBinEndIndex == Integer.MIN_VALUE) {
                        len = 0;
                    }
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

        @Override
        public IndiCount call() throws Exception {
            try {
                Runtime rt = Runtime.getRuntime();
                Process p = rt.exec(command);
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
                StringBuilder baseSb = new StringBuilder();
                for (int i = 0; i < positions.length; i++) {
                    this.setDos(positions[i]);
                    if (current == null) {
                        this.writeMissing();
                    }
                    else {
                        if (positions[i] == currentPosition) {
                            String ref = posRefMap.get(currentPosition);
                            short[] codedAltAlleles = posCodedAlleleMap.get(currentPosition);
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
                                int[] alleleCounts = getAlleleCounts (codedAltAlleles, baseS.toString().toUpperCase(), siteDepth, baseSb);
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
                            System.out.println("Current position is greater than pileup position. It should not happen. Program quits");
                            System.exit(1);
                        }
                    }
                }
                this.closeDos();
                br.close();
                p.waitFor();
                this.writeEmptyFiles();
                System.out.println("Individual allele counting is completed for taxon "+ this.taxonName);
            }
            catch (Exception ee) {
                ee.printStackTrace();
            }
            counter.increment();
            int cnt = counter.intValue();
            if (cnt%50 == 0) System.out.println("Finished individual genotype allele counting in " + String.valueOf(cnt) + " taxa. Total: " + String.valueOf(taxaBamsMap.size()));
            return this;
        }
    }

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

    private String getInfoAndGenotypes (List<short[]> siteCountList, short[] codedAlts) {
        StringBuilder genoSB = new StringBuilder("GT:AD:GL");
        StringBuilder infoSB = new StringBuilder();
        int dp = 0;
        int nz = 0;
        int alleleNumber = codedAlts.length+1;
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
        infoSB.append(";IL=");
        for (int i = 0; i < codedAlts.length; i++) {
            infoSB.append(FastCall2.getIndelLengthFromAlleleCodingLength(codedAlts[i])).append(",");
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

    private String getInfo (String[] genoArray, byte[] codedAlts) {
        int dp = 0;
        int nz = 0;
        int nAlt = codedAlts.length;
        int[] adCnt = new int[1+nAlt];
        int[] acCnt = new int[1+nAlt];
        int[][] gnCnt = new int[1+nAlt][1+nAlt];
        int ht = 0;
        List<String> tempList = null;
        List<String> temList = null;
        for (int i = 0; i < genoArray.length; i++) {
            if (genoArray[i].startsWith(".")) {
                nz++;
                continue;
            }
            tempList = PStringUtils.fastSplit(genoArray[i], ":");
            temList = PStringUtils.fastSplit(tempList.get(1), ",");
            for (int j = 0; j < temList.size(); j++) {
                int c = Integer.parseInt(temList.get(j));
                dp+=c;
                adCnt[j] += c;
            }
            temList = PStringUtils.fastSplit(tempList.get(0), "/");
            for (int j = 0; j < temList.size(); j++) {
                int c = Integer.parseInt(temList.get(j));
                acCnt[c]++;
            }
            int index1 = Integer.parseInt(temList.get(0));
            int index2 = Integer.parseInt(temList.get(1));
            gnCnt[index1][index2]++;
            if (index1 != index2) ht++;
        }
        nz = genoArray.length - nz;
        int sum = 0;
        for (int i = 0; i < acCnt.length; i++) {
            sum+=acCnt[i];
        }
        float maf = (float)((double)acCnt[0]/sum);
        if (maf>0.5) maf = (float)(1-maf);
        StringBuilder sb = new StringBuilder();
        sb.append("DP=").append(dp).append(";NZ=").append(nz).append(";AD=");
        for (int i = 0; i < adCnt.length; i++) {
            sb.append(adCnt[i]).append(",");
        }
        sb.deleteCharAt(sb.length()-1);
        sb.append(";AC=");
        for (int i = 0; i < acCnt.length; i++) {
            sb.append(acCnt[i]).append(",");
        }
        sb.deleteCharAt(sb.length()-1);
        sb.append(";IL=");
        for (int i = 0; i < codedAlts.length; i++) {
            sb.append(FastCall2.getIndelLengthFromAlleleCodingLength(codedAlts[i])).append(",");
        }
        sb.deleteCharAt(sb.length()-1);
        sb.append(";GN=");
        for (int i = 0; i < gnCnt.length; i++) {
            for (int j = i; j < gnCnt.length; j++) {
                sb.append(gnCnt[i][j]).append(",");
            }
        }
        sb.deleteCharAt(sb.length()-1);
        sb.append(";HT=").append(ht).append(";MAF=").append(maf);
        return sb.toString();
    }

    private int[] getAlleleCounts (short[] codedAltAlleles, String baseS, int siteDepth, StringBuilder baseSb) {
        byte[] baseB = baseS.getBytes();
        int[] altAlleleCounts = new int[codedAltAlleles.length];
        int index = Integer.MIN_VALUE;
        int vCnt = 0;
        for (int i = 0; i < baseB.length; i++) {
            byte alleleByte = FastCall2.pileupAscIIToAlleleCodingMap.get(baseB[i]);
            index = Arrays.binarySearch(AlleleEncoder.alleleCodings, alleleByte);
            short queryAlleleShort = -1;
            if (index < 0) continue;
            //weird sign of "^" before a char
            if (i > 0 && baseB[i-1] == 94) continue;
            if (index < 4) {
                queryAlleleShort = FastCall2.getAlleleCodingLength(alleleByte, 0);
            }
            else {
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
                queryAlleleShort = FastCall2.getAlleleCodingLength(alleleByte, length);
                i+=baseSb.length();
                i+=length;
            }

            for (int j = 0; j < codedAltAlleles.length; j++) {
                if (codedAltAlleles[j] == queryAlleleShort) {
                    altAlleleCounts[j]++;
                    vCnt++;
                }
            }
        }
        int[] alleleCounts = new int[codedAltAlleles.length+1];
        alleleCounts[0] = siteDepth - vCnt;
        for (int i = 0; i < codedAltAlleles.length; i++) {
            alleleCounts[i+1] = altAlleleCounts[i];
        }
        return alleleCounts;
    }

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

    private void processVariationLibrary () {
        StringBuilder sb = new StringBuilder();
        sb.append(this.chrom).append("_").append(this.regionStart).append("_").append(regionEnd).append(".pos.txt");
        this.vLibPosFileS = new File (this.outputDirS, sb.toString()).getAbsolutePath();
        this.vl = new VariationLibrary(this.libFileS);
        if (this.chrom != vl.getChrom()) {
            System.out.println("The chromosome number of library and the specified one do not match. Program quits.");
            System.exit(0);
        }
        try {
            vlStartIndex = vl.getStartIndex(this.regionStart);
            vlEndIndex = vl.getEndIndex(this.regionEnd);
            if (vlStartIndex == Integer.MIN_VALUE || vlEndIndex == Integer.MIN_VALUE) {
                System.out.println("The chromosome region was incorrectly set. Program quits.");
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

    private void creatFactorialMap () {
        this.factorialMap = HashIntDoubleMaps.getDefaultFactory().newMutableMap();
        for (int i = 0; i < this.maxFactorial+1; i++) {
            this.factorialMap.put(i, factorial(i));
        }
    }

    public void mkDir () {
        File f = new File (this.outputDirS);
        f.mkdir();
        for (int i = 0; i < subDirS.length; i++) {
            f = new File(outputDirS, subDirS[i]);
            f.mkdir();
        }
    }

    private void parseParameters (List<String> pLineList) {
        this.referenceFileS = pLineList.get(0);
        taxaBamMapFileS = pLineList.get(1);
        this.libFileS = pLineList.get(2);
        String[] tem = pLineList.get(3).split(":");
        this.chrom = Integer.parseInt(tem[0]);
        long start = System.nanoTime();
        System.out.println("Reading reference genome from "+ referenceFileS);
        genomeFa = new FastaBit(referenceFileS);
        System.out.println("Reading reference genome took " + String.format("%.2f", Benchmark.getTimeSpanSeconds(start)) + "s");
        chromIndex = genomeFa.getIndexByName(String.valueOf(this.chrom));
        if (tem.length == 1) {
            this.regionStart = 1;
            this.regionEnd = genomeFa.getSeqLength(chromIndex)+1;
        }
        else if (tem.length == 2) {
            tem = tem[1].split(",");
            this.regionStart = Integer.parseInt(tem[0]);
            this.regionEnd = Integer.parseInt(tem[1])+1;
        }
        this.combinedErrorRate = Double.parseDouble(pLineList.get(4));
        this.samtoolsPath = pLineList.get(5);
        this.threadsNum = Integer.parseInt(pLineList.get(6));
        this.outputDirS = pLineList.get(7);
        this.parseTaxaBamMap(this.taxaBamMapFileS);
    }

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
