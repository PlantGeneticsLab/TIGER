package pgl.app.fastCall2;

import com.mysql.cj.x.protobuf.MysqlxDatatypes;
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

public class BuildVariationLibrary {

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

    int maxAltNum = 2;

    String[] taxaNames = null;

    public BuildVariationLibrary(List<String> pLineList) {
        this.parseParameters(pLineList);
        this.mkLibrary();
    }

    private void mkLibrary () {
        List<File> ingTaxaDirList = IOUtils.getDirListInDir(this.ingDirS);
        taxaNames = new String[ingTaxaDirList.size()];
        for (int i = 0; i < ingTaxaDirList.size(); i++) {
            taxaNames[i] = ingTaxaDirList.get(i).getName();
        }
        Arrays.sort(taxaNames);
        Dyad<int[][], int[]> d = FastCall2.getBins(this.regionStart, this.regionEnd);
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
            }
            VariationLibrary vl = new VariationLibrary (ingList, maoThresh, maxAltNum, chrom, binStarts[i]);
            vlList.add(vl);
        }
        VariationLibrary chromVl = VariationLibrary.getInstance(vlList);
        File f = new File (vLibDirS);
        f.mkdir();
        sb.setLength(0);
        sb.append(chrom).append("_").append(this.regionStart).append("_").append(regionEnd).append(".lib.gz");
        chromVl.writeBinaryFileS(new File (f, sb.toString()).getAbsolutePath());
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

    private void parseParameters (List<String> pLineList) {
        this.referenceFileS = pLineList.get(0);
        String[] tem = pLineList.get(1).split(":");
        this.chrom = Short.parseShort(tem[0]);
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
        this.maoThresh = Integer.parseInt(pLineList.get(2));
        this.threadsNum = Integer.parseInt(pLineList.get(3));
        this.ingDirS = pLineList.get(4);
        this.vLibDirS = pLineList.get(5);
    }
}
