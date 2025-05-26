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

class BuildVariationLibrary extends AppAbstract {

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

    public BuildVariationLibrary(String[] args) {
        this.creatAppOptions();
        this.retrieveAppParameters(args);
        this.mkLibrary();
    }

    @Override
    public void creatAppOptions() {
        options.addOption("app", true, "App name.");
        options.addOption("mod", true, "Module name of FastCall 2.");
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
        System.out.println("Below are the commands of module \"blib\" in FastCall 2.");
        this.printUsage();
    }

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
        List<Future<IndividualGenotype>> futureList = new ArrayList<>();
        List<IndividualGenotype> ingList = new ArrayList<>();
        List<VariationLibrary> vlList = new ArrayList<>();
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
                    Future<IndividualGenotype> result = pool.submit(tr);
                    futureList.add(result);
                }
                pool.shutdown();
                pool.awaitTermination(Long.MAX_VALUE, TimeUnit.SECONDS);

                for (int j = 0; j < futureList.size(); j++) {
                    IndividualGenotype ing = futureList.get(j).get();
                    if (ing == null) continue;
                    ingList.add(ing);
                }
                Collections.sort(ingList);
            }
            catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
            VariationLibrary vl = new VariationLibrary (ingList, maoThresh, FastCall3.maxAltNum, chrom, binStarts[i]);
            vlList.add(vl);
        }
        VariationLibrary chromVl = VariationLibrary.getInstance(vlList);
        File f = new File (vLibDirS);
        f.mkdir();
        sb.setLength(0);
        sb.append(chrom).append("_").append(this.regionStart).append("_").append(regionEnd).append(".lib.gz");
        chromVl.writeBinaryFileS(new File (f, sb.toString()).getAbsolutePath());
        System.out.println("Building genetic variation library is finished.");
    }

    class TaxonRead implements Callable<IndividualGenotype> {
        String fileS;
        public TaxonRead (String fileS) {
            this.fileS = fileS;
        }

        @Override
        public IndividualGenotype call() throws Exception {
            File f = new File (fileS);
            if (!f.exists()) return null;
            IndividualGenotype ing = new IndividualGenotype(this.fileS);
            return ing;
        }
    }
}
