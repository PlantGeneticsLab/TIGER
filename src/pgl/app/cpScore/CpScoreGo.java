/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.app.cpScore;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import pgl.infra.utils.CLIInterface;

/**
 *
 * @author feilu
 */
public class CpScoreGo implements CLIInterface {
    Options options = new Options();
    HelpFormatter optionFormat = new HelpFormatter();
    String introduction = null;
    String mode = null;
    String kmerLengthS = null;
    String referenceGenomeFileS = null;
    String anotherGenomeFileS = null;
    String libFileS = null;
    String outputDirS = null;
    int kmerLength = 0;
    
    public CpScoreGo (String[] args) {
        introduction = this.createIntroduction();
        this.createOptions();
        this.retrieveParameters (args);
        this.runProfiler();
    }

    void runProfiler () {
        if (mode.equals("b")) {
            ReferenceKmerLib lib = new ReferenceKmerLib(kmerLength, referenceGenomeFileS);
            lib.writeBinaryFile(libFileS);
        }
        else if (mode.equals("p")) {
            new GenomeProfiler(libFileS, referenceGenomeFileS, anotherGenomeFileS, outputDirS);
        }
    }
    
    @Override
    public void retrieveParameters (String[] args) {
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine line = parser.parse(options, args);
            mode = line.getOptionValue("m");
            kmerLengthS = line.getOptionValue("k");
            referenceGenomeFileS = line.getOptionValue("r");
            anotherGenomeFileS = line.getOptionValue("a");
            libFileS = line.getOptionValue("l");
            outputDirS = line.getOptionValue("o");
        }
        catch(Exception e) {
            e.printStackTrace();
        }
        if (mode == null) {
            this.printIntroductionAndUsage();
            return;
        }
        else if (mode.equals("b")) {
            if (kmerLengthS == null) {
                this.printIntroductionAndUsage();
                return;
            }
            else if (kmerLengthS.equals("32") || kmerLengthS.equals("16")) {
                kmerLength = Integer.valueOf(kmerLengthS);
            }
            else {
                this.printIntroductionAndUsage();
                return;
            }
        }
        else if (mode.equals("p")) {
            if (kmerLengthS == null) {
                this.printIntroductionAndUsage();
                return;
            }
            else if (kmerLengthS.equals("32") || kmerLengthS.equals("16")) {
                kmerLength = Integer.valueOf(kmerLengthS);
            }
            else {
                this.printIntroductionAndUsage();
                return;
            }
        }
        else {
            this.printIntroductionAndUsage();
            return;
        }
    }
    
    @Override
    public void printIntroductionAndUsage () {
        System.out.println("Incorrect parameter input. Program quits.");
        System.out.println(introduction);
        optionFormat.printHelp("CpScoreProfiler.jar", options );
    }
    
    @Override
    public void createOptions () {
        options = new Options();
        options.addOption("m", true, "Analysis mode. Two modes are available, building reference kmer library (b option) and profiling CpScore (p option). e.g. -m b");
        options.addOption("k", true, "Kmer length. Only 32 and 16 are supported. e.g. -k 32");
        options.addOption("r", true, "Reference genome file. e.g -i maizeAGPV4.fa");
        options.addOption("a", true, "Another genome file from which kmers are counted. e.g -a CML247.fa");
        options.addOption("l", true, "Kmer library file. e.g -l maize_32mer.lib");
        options.addOption("o", true, "Output directory. e.g -o CML247_Cp");
    }
    
    public String createIntroduction () {
        StringBuilder sb = new StringBuilder();
        sb.append("\nThe program CpScoreProfiler.jar is designed to calculate base copy number of a given genome. ");
        sb.append("It has 2 analysis modes. The first is to build kmer library from the reference genome. The second is to calculate the CpScore from another non-reference genome.\n\n");
        sb.append("Command line example:\n\n");
        sb.append("\t1. Build kmer library from reference genome. The reference genome should be in Fasta format. The header must be the chromosome number. e.g. chromosome 1 is >1\n");
        sb.append("\t\tjava -jar CpScoreProfiler.jar -m b -k 32 -r maizeAGPV4.fa -l maize_32mer.lib\n");
        sb.append("\t2. Calculate CpScore from another non-reference genome\n");
        sb.append("\t\tjava -jar CpScoreProfiler.jar -m p -k 32 -r maizeAGPV4.fa -l maize_32mer.lib -a CML247.fa -o CML247_Cp\n");
        return sb.toString();
    }
    
    public static void main (String[] args) { 
        new CpScoreGo (args);
    }
}
