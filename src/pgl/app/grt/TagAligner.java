/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.app.grt;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.InputStreamReader;
import pgl.infra.utils.IOUtils;

/**
 *
 * @author feilu
 */
class TagAligner {
    String referenceFileS = null;
    String bwaPath = null;
    String mergedTagCountFileS = null;
    String alignmentDirS = null;
    String r1FastqFileS = null;
    String r2FastqFileS = null;
    String samFileS = null;
    int numThreads = 32;
    
    public TagAligner (String referenceFileS, String bwaPath, String mergedTagCountFileS, String alignmentDirS) {
        this.referenceFileS = referenceFileS;
        this.bwaPath = bwaPath;
        this.mergedTagCountFileS = mergedTagCountFileS;
        this.alignmentDirS = alignmentDirS;
        this.bwaAlign();
    }
    
    public void setThreads (int numThreads) {
        this.numThreads = numThreads;
    }
    
    public void bwaAlign () {
        r1FastqFileS = new File (alignmentDirS, "tag_r1.fq.gz").getAbsolutePath();
        r2FastqFileS = new File (alignmentDirS, "tag_r2.fq.gz").getAbsolutePath();
        samFileS = new File (alignmentDirS, "tag.sam.gz").getAbsolutePath();
        String perlFileS = new File (alignmentDirS, "runBWA.pl").getAbsolutePath();
        TagAnnotations tc = new TagAnnotations(mergedTagCountFileS);
        tc.writeFastqFile(r1FastqFileS, r2FastqFileS);
        try {
            StringBuilder sb = new StringBuilder();
            int nThreads = Runtime.getRuntime().availableProcessors();
            if (nThreads > this.numThreads) nThreads = this.numThreads;
            sb.append(bwaPath).append(" mem -t ").append(nThreads).append(" ").append(this.referenceFileS).append(" ")
                    .append(r1FastqFileS).append(" ").append(r2FastqFileS).append(" | gzip > ").append(samFileS);      
            String cmd = sb.toString();
            System.out.println(cmd);
            try {
                BufferedWriter bw = IOUtils.getTextWriter(perlFileS);
                bw.write("system (\""+ cmd + "\");");
                bw.newLine();
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            Runtime run = Runtime.getRuntime();
            Process p = run.exec("perl "+perlFileS);
            BufferedReader br = new BufferedReader(new InputStreamReader(p.getErrorStream()));
            String temp = null;
            while ((temp = br.readLine()) != null) {
                System.out.println(temp);
            }
            p.waitFor();
            new File(perlFileS).delete();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
}
