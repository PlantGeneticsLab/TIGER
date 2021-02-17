package pgl.app.fastCall2;

import pgl.AppUtils;
import pgl.PGLConstraints;
import pgl.infra.utils.Dyad;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PArrayUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.File;
import java.util.*;

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
    //Heterozygous ratio (HR) for variation calling, meaning that the depth of alternative allele is greater than HR and less than (1-HR) are considered to be hets.
    double hrThresh = 0.4;
    //Third allele depth ratio (TDR) for variation calling. If the depth of the third allele is greater than TDR by the individual coverage, the site will be ignored. Otherwise, the third allele will be considered as sequencing error.
    double tdrTresh = 0.2;
    //Current chromosome for variation calling
    int chr = Integer.MIN_VALUE;
    //Starting position of the specified region for variation calling, inclusive
    int regionStart = Integer.MIN_VALUE;
    //Ending position the specified regionfor variation calling, exclusive
    int regionEnd = Integer.MIN_VALUE;
    //Number of threads (taxa number to be processed at the same time)
    int threadsNum = PGLConstraints.parallelLevel;
    //Current step ID of the pipeline
    int currentStep = Integer.MIN_VALUE;

    HashMap<String, String[]> taxaBamPathMap = null;
    HashMap<String, Double> taxaCoverageMap = null;

    String[] taxaNames = null;

    public FastCall2 (String parameterFileS) {
        this.parseParameters(parameterFileS);
        if (currentStep == 1) {
            this.variationDiscovery();
        }
        else if (currentStep == 2) {

        }
    }

    private void variationDiscovery () {

    }

    private void parseParameters (String parameterFileS) {
        Dyad<List<String>, List<String>> d = AppUtils.getParameterList(parameterFileS);
        List<String> pLineList = d.getFirstElement();
        List<String> sLineList = d.getSecondElement();
        this.currentStep = Integer.parseInt(sLineList.get(0).split("\\s+")[1]);
        if (currentStep == 1) {
            this.parseParametersStep1(pLineList);
        }
        else if (currentStep == 2) {

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
        this.hrThresh = Double.parseDouble(pLineList.get(7));
        this.tdrTresh = Double.parseDouble(pLineList.get(8));
        String[] tem = pLineList.get(9).split(":");
        this.chr = Integer.parseInt(tem[0]);
        if (tem.length == 2) {
            tem = tem[1].split(",");
            this.regionStart = Integer.parseInt(tem[0]);
            this.regionEnd = Integer.parseInt(tem[1]);
        }
        this.threadsNum = Integer.parseInt(pLineList.get(10));
        this.outputDirS = pLineList.get(11);
        this.samtoolsPath = pLineList.get(11);
        this.getTaxaBamMap(this.taxaRefBamFileS);
    }

    private void getTaxaBamMap (String taxaBamMapFileS) {
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
