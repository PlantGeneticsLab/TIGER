package pgl.app.popdep;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import pgl.AppAbstract;
import pgl.PGLAPPEntrance;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PArrayUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

/**
 * Population depth profiling with sequential BAM I/O.
 *
 * <p>Unlike {@link PopDep}, which re-scans each taxon's BAM once per genomic window
 * (default 500 kb), PopDepFull runs one {@code samtools depth} call per taxon for the
 * entire chromosome region. On HDD this avoids repeated seek-heavy passes over the
 * same large BAM files.</p>
 *
 * <p>Output format is identical to PopDep:
 * {@code Position\tDepth_Mean\tDepth_SD}</p>
 */
public class PopDepFull extends AppAbstract {
    String taxaBamFileS = null;
    short chromosome = Short.MIN_VALUE;
    int chrLength = Integer.MIN_VALUE;
    String samPath = null;
    /** Minimum samtools mapping quality (MAPQ); passed to depth as {@code -q}. */
    int minMapq = 0;
    int threadNum = 16;
    String outfileS = null;
    String[] taxa = null;
    HashMap<String, String[]> taxaBamPathsMap = new HashMap<>();
    /** Per-taxon genome-wide mean depth (column 2 of taxaBamMap), used to normalize to relative depth. */
    HashMap<String, Double> taxaCoverageMap = new HashMap<>();
    /** Coverage aligned to the sorted {@link #taxa} order. */
    double[] taxaCoverage = null;

    /** Chunk size for progress messages while writing output. */
    int progressWindowSize = 500_000;

    public PopDepFull(String[] args) {
        this.creatAppOptions();
        this.retrieveAppParameters(args);
        this.profileDepth();
    }

    public void profileDepth() {
        double[] sum = new double[this.chrLength];
        double[] sumSq = new double[this.chrLength];
        double[] relSum = new double[this.chrLength];
        double[] relSumSq = new double[this.chrLength];
        int nTaxa = this.taxa.length;
        String[] commands = this.getSamCommands(0, this.chrLength);
        int[][] subIndices = PArrayUtils.getSubsetsIndicesBySubsetSize(taxa.length, this.threadNum);

        try {
            for (int u = 0; u < subIndices.length; u++) {
                int batchStart = subIndices[u][0];
                int batchEnd = subIndices[u][1];
                int batchSize = batchEnd - batchStart;
                int[][] batchDepths = new int[batchSize][this.chrLength];

                java.util.stream.IntStream.range(batchStart, batchEnd).parallel().forEach(j -> {
                    try {
                        fillDepthFromSamtools(commands[j], batchDepths[j - batchStart], this.chrLength);
                    } catch (Exception e) {
                        throw new RuntimeException("Depth scan failed for taxon " + taxa[j], e);
                    }
                });

                for (int b = 0; b < batchSize; b++) {
                    int[] d = batchDepths[b];
                    double cov = this.taxaCoverage[batchStart + b];
                    for (int pos = 0; pos < this.chrLength; pos++) {
                        double v = d[pos];
                        sum[pos] += v;
                        sumSq[pos] += v * v;
                        double rv = v / cov;
                        relSum[pos] += rv;
                        relSumSq[pos] += rv * rv;
                    }
                }
                System.out.println("Finished taxa " + batchEnd + " / " + nTaxa + " on chromosome " + this.chromosome);
            }

            StringBuilder sb = new StringBuilder();
            BufferedWriter bw = IOUtils.getTextGzipWriter(this.outfileS);
            bw.write("Position\tDepth_Mean\tDepth_SD\tRelativeDepth_Mean\tRelativeDepth_SD");
            bw.newLine();
            int[][] progressWindows = PArrayUtils.getSubsetsIndicesBySubsetSize(this.chrLength, progressWindowSize);
            for (int i = 0; i < progressWindows.length; i++) {
                int start = progressWindows[i][0];
                int end = progressWindows[i][1];
                for (int pos = start; pos < end; pos++) {
                    sb.setLength(0);
                    double mean = sum[pos] / nTaxa;
                    double var = (sumSq[pos] - sum[pos] * sum[pos] / nTaxa) / (nTaxa - 1);
                    double sd = var > 0 ? Math.sqrt(var) : 0;
                    double relMean = relSum[pos] / nTaxa;
                    double relVar = (relSumSq[pos] - relSum[pos] * relSum[pos] / nTaxa) / (nTaxa - 1);
                    double relSd = relVar > 0 ? Math.sqrt(relVar) : 0;
                    sb.append(pos + 1).append("\t").append((float) mean).append("\t").append((float) sd)
                            .append("\t").append((float) relMean).append("\t").append((float) relSd);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                System.out.println("Current position: " + end + " on chromosome " + this.chromosome);
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("PopDepFull on chromosome " + this.chromosome + " is finished.");
    }

    /**
     * Run one samtools depth command and store per-position depth in {@code depths}.
     * Positions not reported by samtools remain 0 (uncovered).
     */
    private void fillDepthFromSamtools(String command, int[] depths, int expectedLength) throws Exception {
        Runtime rt = Runtime.getRuntime();
        Process p = rt.exec(command);
        try {
            p.getOutputStream().close();
        } catch (Exception ignored) {}

        Thread errDrainer = new Thread(() -> {
            try (BufferedReader ebr = new BufferedReader(new InputStreamReader(p.getErrorStream()))) {
                String errLine;
                while ((errLine = ebr.readLine()) != null) {
                    System.err.println("[samtools std err]: " + errLine);
                }
            } catch (Exception e) {
                throw new RuntimeException(e);
            }
        }, "samtools-stderr-drainer");
        errDrainer.setDaemon(true);
        errDrainer.start();

        try (BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()))) {
            String temp;
            List<String> l;
            while ((temp = br.readLine()) != null) {
                int v = 0;
                l = PStringUtils.fastSplit(temp);
                for (int k = 0; k < l.size() - 2; k += 2) {
                    v += Integer.parseInt(l.get(k + 2));
                }
                int pos = Integer.parseInt(l.get(1));
                if (pos < 1 || pos > expectedLength) continue;
                depths[pos - 1] = v;
            }
        }

        errDrainer.join();
        int exit = p.waitFor();
        if (exit != 0) {
            throw new RuntimeException("Command failed (exit=" + exit + "): " + command);
        }
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
        options.addOption("g", true, "Minimum samtools mapping quality (MAPQ, depth -q). Default 0.");
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
            br.readLine();
            List<String> l;
            String temp;
            while ((temp = br.readLine()) != null) {
                if (temp.trim().isEmpty()) continue;
                l = PStringUtils.fastSplit(temp);
                String[] bams = new String[l.size() - 2];
                for (int i = 0; i < bams.length; i++) {
                    bams[i] = l.get(i + 2);
                }
                this.taxaBamPathsMap.put(l.get(0), bams);
                this.taxaCoverageMap.put(l.get(0), Double.parseDouble(l.get(1)));
            }
            Set<String> tSet = this.taxaBamPathsMap.keySet();
            this.taxa = tSet.toArray(new String[tSet.size()]);
            Arrays.sort(taxa);
            this.taxaCoverage = new double[taxa.length];
            for (int i = 0; i < taxa.length; i++) {
                this.taxaCoverage[i] = this.taxaCoverageMap.get(taxa[i]);
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.out.println("\nThere are input errors in the command line. Program stops.");
            this.printInstructionAndUsage();
            System.exit(0);
        }
    }

    @Override
    public void printInstructionAndUsage() {
        System.out.println(PGLAPPEntrance.getTIGERIntroduction());
        System.out.println("Below are the commands of PopDepFull.");
        this.printUsage();
    }

    private String[] getSamCommands(int startIndex, int endIndex) {
        String[] commands = new String[taxa.length];
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < commands.length; i++) {
            sb.setLength(0);
            sb.append(this.samPath).append(" depth -Q 20 -q ").append(this.minMapq).append(" -r ")
                    .append(this.chromosome).append(":").append(startIndex + 1).append("-").append(endIndex);
            String[] paths = this.taxaBamPathsMap.get(taxa[i]);
            for (String path : paths) {
                sb.append(" ").append(path);
            }
            commands[i] = sb.toString();
        }
        return commands;
    }
}
