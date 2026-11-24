package pgl.app.popdep;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import pgl.AppAbstract;
import pgl.AppUtils;
import pgl.PGLAPPEntrance;
import pgl.infra.dna.FastaBit;
import pgl.infra.table.RowTable;
import pgl.infra.utils.*;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.InputStreamReader;
import java.util.*;

public class PopDep extends AppAbstract {
    /**
     * File path of taxa and there corresponding bams
     */
    String taxaBamFileS = null;
    /**
     * chromosome and it length, used to sampling site to estimate mode
     */
    String chrLengthFileS = null;
    /**
     * Current chromosome for depth profiling
     */
    short chromosome = Short.MIN_VALUE;
    /**
     * Current chromosome length
     */
    int chrLength = Integer.MIN_VALUE;
    /**
     * Path of samtools
     */
    String samPath = null;
    /** Minimum samtools mapping quality (MAPQ); passed to depth/mpileup as {@code -q}. */
    int minMapq = 0;
    /**
     * Number of threads
     */
    int threadNum = 16;
    String outfileS = null;
    String[] taxa = null;
    String[] references = null;
    HashMap<String, String[]> taxaBamPathsMap = new HashMap<>();
    /** Per-taxon genome-wide mean depth (column 2 of taxaBamMap), used to normalize to relative depth. */
    HashMap<String, Double> taxaCoverageMap = new HashMap<>();
    /** Coverage aligned to the sorted {@link #taxa} order. */
    double[] taxaCoverage = null;

    /**
     * Estimated range of mode
     */
    int maxDepthRange = 200;
    /**
     * Current pipeline step
     */
    int step = 0;
    /**
     * Sampling size for estimating mode
     */
    int samplingSize = 10000;
    /**
     * Window size to profile depth
     */
    int windowSize = 500_000;
    
    double[][] depth = null;

    public PopDep (String[] args) {
        this.creatAppOptions();
        this.retrieveAppParameters(args);
        this.profileDepth();
    }

    public void profileDepth () {
        int[][] windows = PArrayUtils.getSubsetsIndicesBySubsetSize(this.chrLength, windowSize);
        int[][] subIndices = PArrayUtils.getSubsetsIndicesBySubsetSize(taxa.length, this.threadNum);
        try {
            StringBuilder sb = new StringBuilder();
            BufferedWriter bw = IOUtils.getTextGzipWriter(this.outfileS);
            bw.write("Position\tDepth_Mean\tDepth_SD\tRelativeDepth_Mean\tRelativeDepth_SD");
            bw.newLine();
            for (int i = 0; i < windows.length; i++) {
                String[] commands = this.getSamCommands(windows[i][0], windows[i][1]);
                depth = new double[windows[i][1]-windows[i][0]][taxa.length];
                int startIndex = windows[i][0];
                int border = windows[i][1];
                for (int u = 0; u < subIndices.length; u++) {
                    List<Integer> indices = PArrayUtils.getIndexList(subIndices[u][0], subIndices[u][1]);
                    indices.parallelStream().forEach(j -> {
                        try {
                            Runtime rt = Runtime.getRuntime();
                            Process p = rt.exec(commands[j]);
                            BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
                            String temp = null;
                            List<String> l = new ArrayList<>();
                            int v = 0;
                            int pos = -1;
                            while ((temp = br.readLine()) != null) {
                                v = 0;
                                l = PStringUtils.fastSplit(temp);
                                for (int k = 0; k < l.size()-2; k+=2) {
                                    v+=Integer.parseInt(l.get(k+2));
                                }
                                pos = Integer.parseInt(l.get(1));
                                if (pos > border) break;
                                depth[pos-1-startIndex][j] = v;
                            }
                            br.close();
                            p.waitFor();
                        }
                        catch (Exception e) {
                            e.printStackTrace();
                        }
                    });
                }
                DescriptiveStatistics d = new DescriptiveStatistics();
                DescriptiveStatistics dRel = new DescriptiveStatistics();
                int blockSize = windows[i][1] - windows[i][0];
                for (int j = 0; j < blockSize; j++) {
                    sb.setLength(0);
                    d = new DescriptiveStatistics(depth[j]);
                    dRel.clear();
                    for (int t = 0; t < taxa.length; t++) {
                        dRel.addValue(depth[j][t] / taxaCoverage[t]);
                    }
                    sb.append(j+windows[i][0]+1).append("\t").append((float)d.getMean()).append("\t").append((float)d.getStandardDeviation())
                            .append("\t").append((float)dRel.getMean()).append("\t").append((float)dRel.getStandardDeviation());
                    bw.write(sb.toString());
                    bw.newLine();
                }
                sb.setLength(0);
                sb.append("Current position: ").append(windows[i][1]).append(" on chromosome ").append(this.chromosome);
                System.out.println(sb.toString());
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("PopDep on chromosome "+String.valueOf(this.chromosome) + " is finished.");
    }

    @Override
    public void creatAppOptions() {
        options.addOption("app", true, "App name.");
        options.addOption("a", true, "The taxaBamMap file: Taxa\\tCoverage\\tBam1[\\tBam2...] per line " +
                "(column 2 is the genome-wide mean depth of the taxon, used for relative depth). " +
                "Bam files must have a .bai index in the same folder.");
        options.addOption("b", true, "The chromosome which will be scanned.");
        options.addOption("c", true, "The length of the chromosome");
        options.addOption("d", true, "The path of samtools.");
        options.addOption("e", true, "Number of threads.");
        options.addOption("f", true, "The output file in gz format");
        options.addOption("g", true, "Minimum samtools mapping quality (MAPQ, depth/mpileup -q). Default 0.");
    }

    @Override
    public void retrieveAppParameters(String[] args) {
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine line = parser.parse(options, args);
            this.taxaBamFileS = line.getOptionValue("a");
            this.chromosome = Short.parseShort(line.getOptionValue("b"));
            this.chrLength = Integer.parseInt(line.getOptionValue("c"));
            this.samPath = line.getOptionValue("d");
            this.threadNum = Integer.parseInt(line.getOptionValue("e"));
            this.outfileS = line.getOptionValue("f");
            if (line.hasOption("g")) {
                this.minMapq = Integer.parseInt(line.getOptionValue("g"));
                if (this.minMapq < 0) {
                    throw new IllegalArgumentException("-g must be >= 0");
                }
            }

            BufferedReader br = IOUtils.getTextReader(this.taxaBamFileS);
            String temp = br.readLine();
            List<String> l;
            while ((temp = br.readLine()) != null) {
                if (temp.trim().isEmpty()) continue;
                l = PStringUtils.fastSplit(temp);
                String[] bams = new String[l.size()-2];
                for (int i = 0; i < bams.length; i++) {
                    bams[i] = l.get(i+2);
                }
                this.taxaBamPathsMap.put(l.get(0), bams);
                this.taxaCoverageMap.put(l.get(0), Double.parseDouble(l.get(1)));
            }
            Set<String> tSet= this.taxaBamPathsMap.keySet();
            this.taxa = tSet.toArray(new String[tSet.size()]);
            Arrays.sort(taxa);
            this.taxaCoverage = new double[taxa.length];
            for (int i = 0; i < taxa.length; i++) {
                this.taxaCoverage[i] = this.taxaCoverageMap.get(taxa[i]);
            }
            this.references = new String[taxa.length];
            br.close();
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
        System.out.println("Below are the commands of PopDep.");
        this.printUsage();
    }

    private String[] getSamCommands(int startIndex, int endIndex) {
        String[] commands = new String[taxa.length];
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < commands.length; i++) {
            sb.setLength(0);
            sb.append(this.samPath).append(" depth -Q 20 -q ").append(this.minMapq).append(" -r ")
                    .append(this.chromosome).append(":").append(startIndex+1).append("-").append(endIndex);
            String[] paths = this.taxaBamPathsMap.get(taxa[i]);
            for (int j = 0; j < paths.length; j++) {
                sb.append(" ").append(paths[j]);
            }
            commands[i] = sb.toString();
        }
        return commands;
    }

   
    private String[] getSamCommand (String siteFileS) {
        String[] commands = new String[taxa.length];
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < commands.length; i++) {
            sb.setLength(0);
            sb.append(this.samPath).append(" mpileup -A -B -Q 20 -q ").append(this.minMapq).append(" -f ")
                    .append(this.references[i]);
            String[] paths = this.taxaBamPathsMap.get(taxa[i]);
            for (int j = 0; j < paths.length; j++) {
                sb.append(" ").append(paths[j]);
            }
            sb.append(" -l ").append(siteFileS);
            commands[i] = sb.toString();
        }
        return commands;
    }

}
