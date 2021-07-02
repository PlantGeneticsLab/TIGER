package pgl.app.popdep;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import pgl.AppUtils;
import pgl.infra.table.RowTable;
import pgl.infra.utils.Dyad;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PArrayUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.InputStreamReader;
import java.util.*;

public class PopDep {
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
    /**
     * Number of threads
     */
    int threadNum = 16;
    String outfileS = null;
    String[] taxa = null;
    String[] references = null;
    HashMap<String, String[]> taxaBamPathsMap = new HashMap<>();

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

    public PopDep (String parameterFileS) {
        this.parseParameters(parameterFileS);
        this.profileDepth();
    }

    public void profileDepth () {
        int[][] windows = PArrayUtils.getSubsetsIndicesBySubsetSize(this.chrLength, windowSize);
        int[][] subIndices = PArrayUtils.getSubsetsIndicesBySubsetSize(taxa.length, this.threadNum);
        try {
            StringBuilder sb = new StringBuilder();
            BufferedWriter bw = IOUtils.getTextGzipWriter(this.outfileS);
            bw.write("Position\tDepth_Mean\tDepth_SD");
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
                int blockSize = windows[i][1] - windows[i][0];
                for (int j = 0; j < blockSize; j++) {
                    sb.setLength(0);
                    d = new DescriptiveStatistics(depth[j]);
                    sb.append(j+windows[i][0]+1).append("\t").append((float)d.getMean()).append("\t").append((float)d.getStandardDeviation());
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

    private String[] getSamCommands(int startIndex, int endIndex) {
        String[] commands = new String[taxa.length];
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < commands.length; i++) {
            sb.setLength(0);
            sb.append(this.samPath).append(" depth -Q 20 -r ").append(this.chromosome).append(":").append(startIndex+1).append("-").append(endIndex);
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
            sb.append(this.samPath).append(" mpileup -A -B -Q 20 -f ").append(this.references[i]);
            String[] paths = this.taxaBamPathsMap.get(taxa[i]);
            for (int j = 0; j < paths.length; j++) {
                sb.append(" ").append(paths[j]);
            }
            sb.append(" -l ").append(siteFileS);
            commands[i] = sb.toString();
        }
        return commands;
    }

    public void parseParameters(String parameterFileS) {
        Dyad<List<String>, List<String>> d = AppUtils.getParameterList(parameterFileS);
        List<String> pLineList = d.getFirstElement();
        this.taxaBamFileS = pLineList.get(0);
        this.chromosome = Short.parseShort(pLineList.get(1));
        this.chrLength = Integer.parseInt(pLineList.get(2));
        this.samPath = pLineList.get(3);
        this.threadNum = Integer.parseInt(pLineList.get(4));
        this.outfileS = pLineList.get(5);
        try {
            BufferedReader br = IOUtils.getTextReader(this.taxaBamFileS);
            String temp = br.readLine();
            List<String> l = new ArrayList<>();
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String[] bams = new String[l.size()-1];
                for (int i = 0; i < bams.length; i++) {
                    bams[i] = l.get(i+1);
                }
                this.taxaBamPathsMap.put(l.get(0), bams);
            }
            Set<String> tSet= this.taxaBamPathsMap.keySet();
            this.taxa = tSet.toArray(new String[tSet.size()]);
            Arrays.sort(taxa);
            this.references = new String[taxa.length];
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

}
