/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.app.grt;

import java.io.File;
import java.util.Arrays;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOUtils;

/**
 *
 * @author feilu
 */
class TagMerger {
    String inputDirS = null;
    String outputFileS = null;
    int collapseThreshold = (int)(Integer.MAX_VALUE * 0.75);
    int minReadCount = 3;
    
    public TagMerger (String inputDirS, String outputFileS) {
        this.mergeTagAnnotations(inputDirS, outputFileS);
    }
    
    public TagMerger (String inputDirS, String outputFileS, int minReadCount) {
        this.minReadCount = minReadCount;
        this.mergeTagAnnotations(inputDirS, outputFileS);
    }
    
    public void mergeTagAnnotations (String inputDirS, String outputFileS) {
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".tas");
        Arrays.sort(fs);
        System.out.println("Merging "+String.valueOf(fs.length)+" individual TagAnnotations files");
        TagAnnotations ta = new TagAnnotations(fs[0].getAbsolutePath());
        boolean ifCollapsed = false;
        int cnt = 0;
        for (int i = 1; i < fs.length; i++) {
            ifCollapsed = false;
            TagAnnotations ata = new TagAnnotations(fs[i].getAbsolutePath());
            ta.addTagAnnotations(ata);
            cnt++;
            if (cnt%100 == 0) {
                ta.collapseCounts(1);
                System.out.println(String.valueOf(cnt) + " TagAnnotations files have been merged");
                System.out.println("Memory used: " + String.format("%.4f", Benchmark.getUsedMemoryGb()) + " Gb\n");
            }
            if (ta.getMaxTagNumberAcrossGroups() < collapseThreshold) continue;
            ta.collapseCounts(minReadCount);
            ifCollapsed = true;
        }
        if (!ifCollapsed) {
            ta.collapseCounts(minReadCount);
        }
        System.out.println("A total of " + String.valueOf(fs.length) + " TagAnnotations files are merged");
        System.out.println(String.valueOf(ta.getTagNumber()) + " tags are devided into " + String.valueOf(ta.getGroupNumber()) + " tag groups");
        ta.writeBinaryFile(outputFileS);
    }
}
