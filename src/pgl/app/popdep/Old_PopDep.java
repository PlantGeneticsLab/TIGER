package pgl.app.popdep;

import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PArrayUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.InputStreamReader;
import java.util.*;

public class Old_PopDep {
    /**
     * File path of taxa and there corresponding bams
     */
    String taxaRefBamFileS = null;
    /**
     * File with info of each taxon and its depth mode
     */
    String taxaDepthModeFileS = null;
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
    double[] mode = null;
    double minMode = 5;
    double[] depWithMinMode = null;
    int[] taxaMinModeIndices = null;

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

    public Old_PopDep(String parameterFileS) {
        this.parseParameters(parameterFileS);
        if (this.step == 1) {
            this.mkTaxaDepthMode();
        }
        else if (this.step == 2) {
            this.profileDepth();
        }
    }

    public void profileDepth () {
        int[][] windows = PArrayUtils.getSubsetsIndicesBySubsetSize(this.chrLength, windowSize);
        int[][] subIndices = PArrayUtils.getSubsetsIndicesBySubsetSize(taxa.length, this.threadNum);
        try {
            StringBuilder sb = new StringBuilder();
            BufferedWriter bw = IOUtils.getTextGzipWriter(this.outfileS);
            bw.write("Position\tDepth_Mean\tDepth_SD\tDepth_Mean_Standardized\tDepth_SD_Standardized");
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
                    Arrays.fill(this.depWithMinMode, 0);
                    for (int k = 0; k < this.taxaMinModeIndices.length; k++) {
                        this.depWithMinMode[k] = depth[j][this.taxaMinModeIndices[k]]/mode[this.taxaMinModeIndices[k]];
                    }

//                    for (int k = 0; k < taxa.length; k++) {
//                        depth[j][k] = depth[j][k]/mode[k];
//                    }
                    d = new DescriptiveStatistics(this.depWithMinMode);
                    sb.append("\t").append((float)d.getMean()).append("\t").append((float)d.getStandardDeviation());
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

    private double getMean (double[] vs) {
        return Arrays.stream(vs).sum()/vs.length;
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

    public void mkTaxaDepthMode () {
        RowTable<String> t = new RowTable<>(this.chrLengthFileS);
        long totalLength = 0;
        short[] chrs = new short[t.getRowNumber()];
        int[] chrLengths = new int[t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            chrs[i] = Short.parseShort(t.getCell(i,0));
            chrLengths[i] = Integer.parseInt(t.getCell(i,1));
            totalLength+=chrLengths[i];
        }
        if (totalLength < samplingSize) totalLength = samplingSize;
        int[][] sites = new int[chrs.length][];
        for (int i = 0; i < sites.length; i++) {
            sites[i] = new int[(int)(samplingSize *((double)chrLengths[i]/totalLength))];
            double v = (double)chrLengths[i]/sites[i].length;
            for (int j = 0; j < sites[i].length; j++) {
                sites[i][j] = (int)(j*v+1);
            }
        }
        String siteFileS = new File(new File(this.taxaDepthModeFileS).getParent(), "sites.txt").getAbsolutePath();
        try {
            BufferedWriter bw = IOUtils.getTextWriter(siteFileS);
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < sites.length; i++) {
                for (int j = 0; j < sites[i].length; j++) {
                    sb.setLength(0);
                    sb.append(chrs[i]).append("\t").append(sites[i][j]);
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        String[] commands = this.getSamCommand(siteFileS);
        int[] modes = new int[commands.length];
        int[][] subIndices = PArrayUtils.getSubsetsIndicesBySubsetSize(commands.length, this.threadNum);
        for (int i = 0; i < subIndices.length; i++) {
            List<Integer> indices = PArrayUtils.getIndexList(subIndices[i][0], subIndices[i][1]);
            indices.parallelStream().forEach(j -> {
                int[] depth = new int[maxDepthRange+1];
                try {
                    Runtime rt = Runtime.getRuntime();
                    Process p = rt.exec(commands[j]);
                    BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
                    String temp = null;
                    List<String> l = new ArrayList<>();
                    int v = 0;
                    while ((temp = br.readLine()) != null) {
                        v = 0;
                        l = PStringUtils.fastSplit(temp);
                        for (int k = 0; k < l.size()-3; k+=3) {
                            v+=Integer.parseInt(l.get(k+3));
                        }
                        if (v > this.maxDepthRange) continue;
                        depth[v]++;
                    }
                    br.close();
                    p.waitFor();
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
                int maxValue = Integer.MIN_VALUE;
                int maxDepth = Integer.MAX_VALUE;
                for (int k = 0; k < depth.length; k++) {
                    if (depth[k] > maxValue) {
                        maxValue = depth[k];
                        maxDepth = k;
                    }
                }
                modes[j] = maxDepth;
            });
            System.out.println(String.valueOf(subIndices[i][1])+ " depth calculation finished in Step 1");
        }
        try {
            BufferedWriter bw = IOUtils.getTextWriter(this.taxaDepthModeFileS);
            StringBuilder sb = new StringBuilder();
            bw.write("Taxa\tMode");
            bw.newLine();
            for (int i = 0; i < taxa.length; i++) {
                sb.setLength(0);
                sb.append(taxa[i]).append("\t").append(modes[i]);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        new File(siteFileS).delete();
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
        List<String> pLineList = new ArrayList<>();
        try {
            BufferedReader br = IOUtils.getTextReader(parameterFileS);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) continue;
                if (temp.startsWith("@")) {
                    this.step = Integer.parseInt(temp.split("\t")[1]);
                    continue;
                }
                if (temp.isEmpty()) continue;
                pLineList.add(temp);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        if (step == 1) {
            this.parseStep1(pLineList);
        }
        else if (step == 2) {
            this.parseStep2(pLineList);
        }
    }

    private void parseStep1(List<String> pLineList) {
        this.taxaRefBamFileS = pLineList.get(0);
        this.taxaDepthModeFileS = pLineList.get(1);
        this.chrLengthFileS = pLineList.get(2);
        this.samPath = pLineList.get(3);
        this.threadNum = Integer.parseInt(pLineList.get(4));
        HashMap<String, String> taxaRefMap = new HashMap<>();
        try {
            BufferedReader br = IOUtils.getTextReader(this.taxaRefBamFileS);
            String temp = br.readLine();
            List<String> l = new ArrayList<>();
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String[] bams = new String[l.size()-2];
                for (int i = 0; i < bams.length; i++) {
                    bams[i] = l.get(i+2);
                }
                this.taxaBamPathsMap.put(l.get(0), bams);
                taxaRefMap.put(l.get(0), l.get(1));
            }
            Set<String> tSet= this.taxaBamPathsMap.keySet();
            this.taxa = tSet.toArray(new String[tSet.size()]);
            Arrays.sort(taxa);
            this.references = new String[taxa.length];
            for (int i = 0; i < taxa.length; i++) {
                this.references[i] = taxaRefMap.get(taxa[i]);
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void parseStep2(List<String> pLineList) {
        this.taxaRefBamFileS = pLineList.get(0);
        this.taxaDepthModeFileS = pLineList.get(1);
        this.chromosome = Short.parseShort(pLineList.get(2));
        this.chrLength = Integer.parseInt(pLineList.get(3));
        this.minMode = Double.parseDouble(pLineList.get(4));
        this.samPath = pLineList.get(5);
        this.threadNum = Integer.parseInt(pLineList.get(6));
        HashMap<String, String> taxaRefMap = new HashMap<>();
        this.outfileS = pLineList.get(7);
        try {
            BufferedReader br = IOUtils.getTextReader(this.taxaRefBamFileS);
            String temp = br.readLine();
            List<String> l = new ArrayList<>();
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String[] bams = new String[l.size()-2];
                for (int i = 0; i < bams.length; i++) {
                    bams[i] = l.get(i+2);
                }
                this.taxaBamPathsMap.put(l.get(0), bams);
                taxaRefMap.put(l.get(0), l.get(1));
            }
            Set<String> tSet= this.taxaBamPathsMap.keySet();
            this.taxa = tSet.toArray(new String[tSet.size()]);
            Arrays.sort(taxa);
            this.references = new String[taxa.length];
            for (int i = 0; i < taxa.length; i++) {
                this.references[i] = taxaRefMap.get(taxa[i]);
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        HashMap<String, Double> taxaModeMap = new HashMap<>();
        RowTable<String> t = new RowTable<>(this.taxaDepthModeFileS);
        TIntArrayList minModeIndexList = new TIntArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            double mode = Double.parseDouble(t.getCell(i,1));
            taxaModeMap.put(t.getCell(i,0), mode);
            if (mode < this.minMode) continue;
            minModeIndexList.add(i);
        }
        this.depWithMinMode = new double[minModeIndexList.size()];
        this.taxaMinModeIndices = minModeIndexList.toArray();
        this.mode = new double[this.taxa.length];
        for (int i = 0; i < taxa.length; i++) {
            mode[i] = taxaModeMap.get(taxa[i]);
        }
    }
}
