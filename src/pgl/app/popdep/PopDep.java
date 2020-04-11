package pgl.app.popdep;

import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PArrayUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.InputStreamReader;
import java.util.*;

public class PopDep {
    String taxaBamFileS = null;
    String taxaDepthModeFileS = null;
    boolean ifModeExist = false;
    short chromosome = Short.MIN_VALUE;
    int chrLength = Integer.MIN_VALUE;
    String samPath = null;
    int threadNum = 16;
    String outputDirS = null;
    String[] taxa = null;
    HashMap<String, String[]> taxaBamPathsMap = new HashMap<>();
    HashMap<String, Double> taxaModeMap = new HashMap<>();
    int maxDepthRange = 200;

    public PopDep (String parameterFileS) {
        this.parseParameters(parameterFileS);
        if (!this.ifModeExist) this.mkTaxaDepthMode();
    }

    public void mkTaxaDepthMode () {
        int size = 50000;
        int[] sites = null;
        if (size > this.chrLength) {
            size = chrLength;
            sites = new int[size];
            for (int i = 0; i < size; i++) {
                sites[i] = i + 1;
            }
        }
        else {
            double v = (double)this.chrLength/size;
            sites = new int[size];
            for (int i = 0; i < sites.length; i++) {
                sites[i] = (int)(i*v+1);
            }
        }
        String siteFileS = new File(this.outputDirS, "sites.txt").getAbsolutePath();
        try {
            BufferedWriter bw = IOUtils.getTextWriter(siteFileS);
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < sites.length; i++) {
                sb.setLength(0);
                sb.append(this.chromosome).append("\t").append(sites[i]);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        String[] commands = this.getSamCommand(siteFileS);
        int[] modes = new int[commands.length];
        List<Integer> indices = PArrayUtils.getIndexList(modes.length);
        indices.parallelStream().forEach(i -> {
            int[] depth = new int[maxDepthRange+1];
            try {
                Runtime rt = Runtime.getRuntime();
                Process p = rt.exec(commands[i]);
                BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
                String temp = null;
                List<String> l = new ArrayList<>();
                int v = 0;
                while ((temp = br.readLine()) != null) {
                    v = 0;
                    l = PStringUtils.fastSplit(temp);
                    for (int j = 0; j < l.size()-2; j+=2) {
                        v+=Integer.parseInt(l.get(j+2));
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
            for (int j = 0; j < depth.length; j++) {
                if (depth[j] > maxValue) {
                    maxValue = depth[j];
                    maxDepth = j;
                }
            }
            modes[i] = maxDepth;
        });
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
    }

    private String[] getSamCommand (String siteFileS) {
        String[] commands = new String[taxa.length];
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < commands.length; i++) {
            sb.setLength(0);
            sb.append(this.samPath).append(" depth -Q 20 -b ").append(siteFileS);
            String[] paths = this.taxaBamPathsMap.get(taxa[i]);
            for (int j = 0; j < paths.length; j++) {
                sb.append(" ").append(paths[j]);
            }
            commands[i] = sb.toString();
        }
        return commands;
    }

    public void parseParameters(String parameterFileS) {
        List<String> pLineList = new ArrayList<>();
        try {
            BufferedReader br = IOUtils.getTextReader(parameterFileS);
            String temp = null;
            boolean ifOut = false;
            if (!(temp = br.readLine()).equals("@App:\tPopDep")) ifOut = true;
            if (!(temp = br.readLine()).equals("@Author:\tFei Lu")) ifOut = true;
            if (!(temp = br.readLine()).equals("@Email:\tflu@genetics.ac.cn; dr.lufei@gmail.com")) ifOut = true;
            if (!(temp = br.readLine()).equals("@Homepage:\thttps://plantgeneticslab.weebly.com/")) ifOut = true;
            if (ifOut) {
                System.out.println("Thanks for using HapScanner.");
                System.out.println("Please keep the authorship in the parameter file. Program stops.");
                System.exit(0);
            }
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) continue;
                if (temp.isEmpty()) continue;
                pLineList.add(temp);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        this.taxaBamFileS = pLineList.get(0);
        this.taxaDepthModeFileS = pLineList.get(1);
        int v = Integer.parseInt(pLineList.get(2));
        if (v == 0) this.ifModeExist = false;
        else if (v == 1) this.ifModeExist = true;
        this.chromosome = Short.parseShort(pLineList.get(3));
        this.chrLength = Integer.parseInt(pLineList.get(4));
        this.samPath = pLineList.get(5);
        this.threadNum = Integer.parseInt(pLineList.get(6));
        this.outputDirS = pLineList.get(7);
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
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
