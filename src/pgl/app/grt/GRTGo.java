/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.app.grt;

import pgl.infra.pos.ChrPos;
import java.io.File;
import java.util.List;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.CLIInterface;

/**
 *
 * @author feilu
 */
public class GRTGo implements CLIInterface {
    Options options = new Options();
    HelpFormatter optionFormat = new HelpFormatter();
    String introduction = null;
    String mode = null;
    String workingDirS = null;
    String barcodeFileS = null;
    String libraryFqMapFileS = null;
    String enzymeCutterF = null;
    String enzymeCutterR = null;
    String bwaPath = null;
    String referenceFileS = null;
    int numThreads = 32;
    int minReadCount = 3;
    int minReadCountAlt = 10;
    int minMappingQ = 30;
    int maxMappingLength = 1000;
    int maxDivergence = 7;
    int tagIdentifyThreshold = 1;
    
    LibraryInfo li = null;
    String[] subDirS = {"tagsBySample","tagsLibrary","alignment", "rawGenotype", "filteredGenotype", "queryGenotype"};
    
    public GRTGo (String[] args) {
        this.createOptions();
        introduction = this.createIntroduction();
        this.retrieveParameters (args);
    }
    
    @Override
    public void retrieveParameters (String[] args) {
        long start = System.nanoTime();
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine line = parser.parse(options, args);
            mode = line.getOptionValue("m");
            workingDirS = line.getOptionValue("w");
            barcodeFileS = line.getOptionValue("b");
            libraryFqMapFileS = line.getOptionValue("f");
            enzymeCutterF = line.getOptionValue("ef");
            enzymeCutterR = line.getOptionValue("er");
            String temp = line.getOptionValue("t");
            int numCores = Runtime.getRuntime().availableProcessors();
            if (temp != null) {
                int inputThreads = Integer.parseInt(temp);
                numThreads = inputThreads;
                if (inputThreads < 0) numThreads = numCores;
            }
            if (numThreads > numCores) numThreads = numCores;
            temp = line.getOptionValue("mc");
            if (temp != null) {
                this.minReadCount = Integer.parseInt(temp);
            }
            temp = line.getOptionValue("mca");
            if (temp != null) {
                this.minReadCountAlt = Integer.parseInt(temp);
            }
            this.referenceFileS = line.getOptionValue("g");
            this.bwaPath = line.getOptionValue("bwa");
            temp = line.getOptionValue("mq");
            if (temp != null) {
                this.minMappingQ = Integer.parseInt(temp);
            }
            temp = line.getOptionValue("ml");
            if (temp != null) {
                this.maxMappingLength = Integer.parseInt(temp);
            }
            temp = line.getOptionValue("md");
            if (temp != null) {
                this.maxDivergence = Integer.parseInt(temp);
            }
            temp = line.getOptionValue("it");
            if (temp != null) {
                this.tagIdentifyThreshold = Integer.parseInt(temp);
            }
        }
        catch(Exception e) {
            e.printStackTrace();
            System.exit(1); 
        }
        this.creatDirectories();
        if (mode == null) {
            this.printIntroductionAndUsage();
            System.exit(1); 
            return;
        }
        else if (mode.equals("pf")) {
            if (workingDirS == null || barcodeFileS == null || libraryFqMapFileS == null || enzymeCutterF == null || enzymeCutterR == null) {
                this.printIntroductionAndUsage();
                System.exit(1); 
                return;
            }
            System.out.println("Start parsing fastq files based on barcodes...");
            li = new LibraryInfo(barcodeFileS, libraryFqMapFileS, enzymeCutterF, enzymeCutterR);
            String tagBySampleDirS = new File (this.workingDirS, this.subDirS[0]).getAbsolutePath();
            TagParser tp = new TagParser(li);
            tp.setThreads(numThreads);
            tp.parseFastq(tagBySampleDirS);
            tp.compressTagsBySample(tagBySampleDirS);
            System.out.println("Parsing fastq is complemeted in " + String.format("%.4f", Benchmark.getTimeSpanHours(start)) + " hours");
        }
        else if (mode.equals("mt")) {
            if (workingDirS == null) {
                this.printIntroductionAndUsage();
                System.exit(1); 
                return;
            }
            System.out.println("Start merging tags in each individual TagAnnotations file...");
            String tagBySampleDirS = new File (this.workingDirS, this.subDirS[0]).getAbsolutePath();
            String tagLibraryDirS = new File (this.workingDirS, this.subDirS[1]).getAbsolutePath();
            String mergedTagCountFileS = new File(tagLibraryDirS, "tag.tas").getAbsolutePath();
            TagMerger tm = new TagMerger(tagBySampleDirS, mergedTagCountFileS, minReadCount);
            System.out.println("Merging tags is complemeted in " + String.format("%.4f", Benchmark.getTimeSpanHours(start)) + " hours");
        }
        else if (mode.equals("at")) {
            if (workingDirS == null || this.referenceFileS == null || bwaPath == null) {
                this.printIntroductionAndUsage();
                System.exit(1); 
                return;
            }
            System.out.println("Start aligning tags using bwa...");
            String tagLibraryDirS = new File (this.workingDirS, this.subDirS[1]).getAbsolutePath();
            String mergedTagAnnotationFileS = new File(tagLibraryDirS, "tag.tas").getAbsolutePath();
            String alignmentDirS = new File (this.workingDirS, this.subDirS[2]).getAbsolutePath();
            TagAligner ta = new TagAligner(referenceFileS, this.bwaPath, mergedTagAnnotationFileS, alignmentDirS);
            ta.setThreads(numThreads);
            System.out.println("Aligning tags is complemeted in " + String.format("%.4f", Benchmark.getTimeSpanHours(start)) + " hours");
        }
        else if (mode.equals("cs")) {
            if (workingDirS == null) {
                this.printIntroductionAndUsage();
                System.exit(1); 
                return;
            }
            System.out.println("Start calling SNPs...");
            String tagLibraryDirS = new File (this.workingDirS, this.subDirS[1]).getAbsolutePath();
            String tagAnnotationFileS = new File(tagLibraryDirS, "tag.tas").getAbsolutePath();
            String alignmentDirS = new File (this.workingDirS, this.subDirS[2]).getAbsolutePath();
            String rawSNPFileS = new File(tagLibraryDirS, "rawSNP.bin").getAbsolutePath();
            String samFileS = new File (alignmentDirS, "tag.sam.gz").getAbsolutePath();
            TagAnnotations tas = new TagAnnotations(tagAnnotationFileS);
            tas.removeAllSNP();
            tas.callSNP(samFileS, this.minMappingQ, this.maxMappingLength, this.maxDivergence);
            tas.writeBinaryFile(tagAnnotationFileS);
            SNPCounts snpSCs = tas.getSNPCounts();
            snpSCs.writeBinaryFile(rawSNPFileS);
            System.out.println("Calling SNPs is complemeted in " + String.format("%.4f", Benchmark.getTimeSpanHours(start)) + " hours");
        }
        else if (mode.equals("rs")) {
            if (workingDirS == null) {
                this.printIntroductionAndUsage();
                System.exit(1); 
                return;
            }
            System.out.println("Start removing low read count SNPs...");
            String tagLibraryDirS = new File (this.workingDirS, this.subDirS[1]).getAbsolutePath();
            String rawSNPFileS = new File(tagLibraryDirS, "rawSNP.bin").getAbsolutePath();
            SNPCounts snpSCs = new SNPCounts(rawSNPFileS);
            snpSCs.writeBinaryFile(rawSNPFileS, this.minReadCountAlt);
            System.out.println("Removing low read count SNPs is complemeted in " + String.format("%.4f", Benchmark.getTimeSpanHours(start)) + " hours");
        }
        else if (mode.equals("ca")) {
            if (workingDirS == null) {
                this.printIntroductionAndUsage();
                System.exit(1); 
                return;
            }
            System.out.println("Start calling alleles...");
            String tagLibraryDirS = new File (this.workingDirS, this.subDirS[1]).getAbsolutePath();
            String tagAnnotationFileS = new File(tagLibraryDirS, "tag.tas").getAbsolutePath();
            String alignmentDirS = new File (this.workingDirS, this.subDirS[2]).getAbsolutePath();
            String rawSNPFileS = new File(tagLibraryDirS, "rawSNP.bin").getAbsolutePath();
            String samFileS = new File (alignmentDirS, "tag.sam.gz").getAbsolutePath();
            TagAnnotations tas = new TagAnnotations(tagAnnotationFileS);
            tas.removeAllAllele();
            SNPCounts sc = new SNPCounts (rawSNPFileS);
            tas.callAllele(samFileS, sc, this.minMappingQ, this.maxMappingLength);
            tas.writeBinaryFile(tagAnnotationFileS);
            System.out.println("Calling alleles is complemeted in " + String.format("%.4f", Benchmark.getTimeSpanHours(start)) + " hours");
        }
        else if (mode.equals("cg")) {
            if (workingDirS == null) {
                this.printIntroductionAndUsage();
                System.exit(1); 
                return;
            }
            System.out.println("Start calling genotypes...");
            String tagBySampleDirS = new File (this.workingDirS, this.subDirS[0]).getAbsolutePath();
            String tagLibraryDirS = new File (this.workingDirS, this.subDirS[1]).getAbsolutePath();
            String tagAnnotationFileS = new File(tagLibraryDirS, "tag.tas").getAbsolutePath();
            String rawSNPFileS = new File(tagLibraryDirS, "rawSNP.bin").getAbsolutePath();
            String genotypeDirS = new File (this.workingDirS, this.subDirS[3]).getAbsolutePath();
            TagAnnotations tas = new TagAnnotations(tagAnnotationFileS);
            SNPCounts sc = new SNPCounts (rawSNPFileS);
            GBSVCFBuilder builder = new GBSVCFBuilder(tas, sc);
            builder.setThreads(numThreads);
            builder.setTagIdentifyThreshold(tagIdentifyThreshold);
            builder.callGenotype(tagBySampleDirS, genotypeDirS);
            System.out.println("Calling genotype is complemeted in " + String.format("%.4f", Benchmark.getTimeSpanHours(start)) + " hours");
        }
        else if (mode.equals("fd")) {
            System.out.println("Start filtering database with validated genotype");  
            String tagLibraryDirS = new File (this.workingDirS, this.subDirS[1]).getAbsolutePath();
            String tagAnnotationFileS = new File(tagLibraryDirS, "tag.tas").getAbsolutePath();
            String rawSNPFileS = new File(tagLibraryDirS, "rawSNP.bin").getAbsolutePath();
            String filteredGenotypeDirS = new File (this.workingDirS, this.subDirS[4]).getAbsolutePath();
            String filteredTagAnnotationFileS = new File(tagLibraryDirS, "tag.filtered.tas").getAbsolutePath();
            String filteredSNPFileS = new File(tagLibraryDirS, "SNP.filtered.bin").getAbsolutePath();
            TagAnnotations tas = new TagAnnotations(tagAnnotationFileS);
            SNPCounts sc = new SNPCounts (rawSNPFileS);
            List<ChrPos> validatedSNPPosList = tas.filterTagAnnotationsWithValidatedGenotype(filteredGenotypeDirS);
            tas.writeBinaryFile(filteredTagAnnotationFileS);
            sc.selectSNPs(validatedSNPPosList);
            sc.writeBinaryFile(filteredSNPFileS);
            System.out.println("Filering database is complemeted in " + String.format("%.4f", Benchmark.getTimeSpanHours(start)) + " hours");
        }
        else if (mode.equals("rg")) {
            System.out.println("Start retrieving genotype from database");  
            String tagBySampleDirS = new File (this.workingDirS, this.subDirS[0]).getAbsolutePath();
            String tagLibraryDirS = new File (this.workingDirS, this.subDirS[1]).getAbsolutePath();
            String filteredTagAnnotationFileS = new File(tagLibraryDirS, "tag.filtered.tas").getAbsolutePath();
            String filteredSNPFileS = new File(tagLibraryDirS, "SNP.filtered.bin").getAbsolutePath();
            String genotypeDirS = new File (this.workingDirS, this.subDirS[5]).getAbsolutePath();
            TagAnnotations tas = new TagAnnotations(filteredTagAnnotationFileS);
            SNPCounts sc = new SNPCounts (filteredSNPFileS);
            GBSVCFBuilder builder = new GBSVCFBuilder(tas, sc);
            builder.setTagIdentifyThreshold(3);
            builder.callGenotype(tagBySampleDirS, genotypeDirS);
            System.out.println("Retrieving genotype is complemeted in " + String.format("%.4f", Benchmark.getTimeSpanHours(start)) + " hours");
        }
        else {
            this.printIntroductionAndUsage();
            System.exit(1); 
            return;
        }
    }
    
    private void creatDirectories () {
        File workingDir = new File(this.workingDirS);
        workingDir.mkdir();
        for (int i = 0; i < this.subDirS.length; i++) {
            File f = new File (this.workingDirS, subDirS[i]);
            f.mkdir();
        }
    }
    
    @Override
    public void createOptions () {
        options = new Options();
        options.addOption("m", true, "Analysis mode.");
        options.addOption("w", true, "Working directory, where sub-directories are created for analysis.");
        options.addOption("b", true, "The barcode file, where sample barcoding information is stored.");
        options.addOption("f", true, "The libraryFastqMap file, where corresponding fastq files can be found for each flowcell_lane_library-index combination.");
        options.addOption("ef", true, "Recognition sequence of restriction enzyme in R1, e.g GGATCC");
        options.addOption("er", true, "Recognition sequence of restriction enzyme in R2, e.g CCGG");
        options.addOption("t", true, "Number of threads. The default value is 32. The actual number of running threads is less than the number of cores regardless of the input value, but -1 means the number of all available cores");
        options.addOption("g", true, "The reference genome of the species. The indexing files should be included in the same directory of the reference genome.");
        options.addOption("bwa", true, "The path of bwa executable file, e.g /Users/Software/bwa-0.7.15/bwa");
        options.addOption("mc", true, "The minimum read count of tag in database. The default value is 3.");
        options.addOption("mca", true, "The minimum read count of an alternative allele. The default value is 10.");
        options.addOption("mq", true, "The minimum read mapping quality for SNP calling and allele calling. The default value is 30.");
        options.addOption("ml", true, "The maximum range of paired-end read mapping. The default value is 1000.");
        options.addOption("md", true, "The maximum divergence between a tag and the reference genome, which is a quality control in SNP calling. The default value is 7.");
        options.addOption("it", true, "The tag identify threshold. While searching the tag DB, query tag having more mismatch than the value is not considered as a match. The default value is 3.");
    }
    
    @Override
    public String createIntroduction () {
        StringBuilder sb = new StringBuilder();
        sb.append("\nThe program GRT.jar is designed to genotype seuqecing samples made from two-enzyme GBS systems. ");
        sb.append("By using a large set of diverse samples, it builds up a database of genetic variants of a species.\n");
        sb.append("The genotype of a tested sample can be retrieved from the database.\n\n");
        sb.append("Command line example:\n\n");
        sb.append("java -Xms20g -Xmx50g -GRT.jar " +
                  "-m pf "
                + "-w /Users/user1/Lib_GBS/pipeOutput/ "
                + "-b /Users/user1/Lib_GBS/source/20180601_GBSLibrarybarcode.txt "
                + "-f /Users/user1/Lib_GBS/source/LibraryFastqMap.txt "
                + "-ef GGATCC "
                + "-er CCGG\n");
        return sb.toString();
    }
    
    public static void main (String[] args) {
        new GRTGo (args);
    }

    @Override
    public void printIntroductionAndUsage() {
        System.out.println("Incorrect parameter input. Program stops.");
        System.out.println(introduction);
        optionFormat.printHelp("GRT.jar", options);
    }
    
}
