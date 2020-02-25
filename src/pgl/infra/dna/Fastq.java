/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package pgl.infra.dna;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.zip.GZIPOutputStream;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;


/**
 * Hold Fastq single end file. 
 * <p>
 * Due to the large size of Fastq file, it is recommended to load a small chunk of the file.
 * <p>
 * The beginning of a Fastq file generally has lower quality, it is recommended to skip the first 100,000 reads
 * @author Fei Lu
 */
public class Fastq {
    private int minStartIndex = 100000;
    private int maxReadNum = 1000000;
    Read[] reads = null;
    int phredScale = Integer.MIN_VALUE;
    
    /**
     * Constructor to sample Fastq file, ignore those low-quality sequence at the beginning of a fastq file
     * @param fastqFileS
     * @param startIndex 100,000 by default
     * @param readNum The maximum is 1,000,000 by default
     */
    public Fastq (String fastqFileS, int startIndex, int readNum) {
        if (startIndex < minStartIndex) {
            startIndex = minStartIndex;
            System.out.println("Start index of read was set to " + String.valueOf(startIndex));
        }
        if (readNum > maxReadNum) {
            readNum = maxReadNum;
            System.out.println("Number of read was set to " + String.valueOf(readNum));
        }
        if (fastqFileS.endsWith(".gz")) {
            this.readFastq(fastqFileS, IOFileFormat.TextGzip, startIndex, readNum);
        }
        else {
            this.readFastq(fastqFileS, IOFileFormat.Text, startIndex, readNum);
        }
        
    }
    
    /**
     * Constructor to read in whole Fastq, fastq file should be small for test
     * @param fastqFileS 
     */
    public Fastq (String fastqFileS) {
        if (fastqFileS.endsWith("gz")) {
            this.readFastq(fastqFileS, IOFileFormat.TextGzip);
        }
        else {
            this.readFastq(fastqFileS, IOFileFormat.Text);
        }
    }
    
    public Fastq (Read[] reads, int phredScore) {
        this.reads = reads;
        this.phredScale = phredScore;
    }
    
    
    /**
     * Return phred score scale of the Fastq file, 33 or 64
     * @return 
     */
    public int getPhredScale () {
        return this.phredScale;
    }
    
    private void readFastq (String fastqFileS, IOFileFormat format) {
        System.out.println("Reading fastq file from " + fastqFileS);
        BufferedReader br = null;
        int cnt = 0;
        String temp;
        try {
            if (format == IOFileFormat.Text) {
                br = IOUtils.getTextReader(fastqFileS);
                this.setPhredScale(br, 10);
                br = IOUtils.getTextReader(fastqFileS);
            }
            else if (format == IOFileFormat.TextGzip) {
                br = IOUtils.getTextGzipReader(fastqFileS);
                this.setPhredScale(br, 10);
                br = IOUtils.getTextGzipReader(fastqFileS);
            }
            else {}
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        this.readFastq(br);
    }
    
    private void readFastq (String fastqFileS, IOFileFormat format, int startIndex, int readNum) {
        if (readNum <= 0) return;
        BufferedReader br = null;
        try {
            if (format == IOFileFormat.Text) {
                br = IOUtils.getTextReader(fastqFileS);
                this.setPhredScale(br, readNum);
                br = IOUtils.getTextReader(fastqFileS);
            }
            else if (format == IOFileFormat.TextGzip) {
                br = IOUtils.getTextGzipReader(fastqFileS);
                this.setPhredScale(br, readNum);
                br = IOUtils.getTextGzipReader(fastqFileS);
            }
            else {}
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("Reading fastq file from " + fastqFileS);
        this.readFastq(br, startIndex, readNum);
    }
    
    private void setPhredScale (BufferedReader br, int readNum) {
        try {
            int size = 10;
            if (readNum < 10) size = readNum;
            for (int i = 0; i < size; i++) {
                br.readLine();br.readLine();br.readLine();
                byte[] qualB = br.readLine().getBytes();
                for (int j = qualB.length-1; j > -1; j--) {
                    if (qualB[j] < 65) {
                        this.phredScale = 33;
                        return;
                    }
                }
            }
            this.phredScale = 64;
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private void readFastq (BufferedReader br) {
        int readCount = 0;
        try {
            Read r = null;
            String temp;
            ArrayList<Read> rList = new ArrayList();
            while ((temp = br.readLine()) != null) {
                r = new Read (temp, br.readLine(), br.readLine(), br.readLine(), this.getPhredScale());
                rList.add(r);
                readCount++;
            }
            this.reads = rList.toArray(new Read[rList.size()]);
            br.close();
            
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println(String.valueOf(readCount) + " reads imported");
    }
    
    private void readFastq (BufferedReader br, int startIndex, int readNum) {
        this.reads = new Read[readNum];
        try {
            Read r = null;
            String temp;
            int index = 0;
            int readCount = 0;
            while ((temp = br.readLine()) != null) {
                if (index >= startIndex) {
                    r = new Read (temp, br.readLine(), br.readLine(), br.readLine(), this.getPhredScale());
                    reads[readCount] = r;
                    readCount++;
                    if (readCount == readNum) break;
                    
                }
                index++;
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println(String.valueOf(readNum) + " reads imported");
    }
    
    /**
     * Write Fastq file
     * @param outputFileS
     * @param format 
     */
    public void writeFastq (String outputFileS, IOFileFormat format) {
        BufferedWriter bw = null;
        try {
            if (format == IOFileFormat.Text) {
                bw = new BufferedWriter(new FileWriter(outputFileS), 65536);
            }
            else if (format == IOFileFormat.TextGzip) { 
                bw = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outputFileS), 65536)), 65536);
            }
            else {}
            for (int i = 0; i < this.getReadNum(); i++) {
                bw.write(reads[i].getID());
                bw.newLine();
                bw.write(reads[i].getSequence());
                bw.newLine();
                bw.write(reads[i].getDescription());
                bw.newLine();
                bw.write(reads[i].getQualS(this.getPhredScale()));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("Fastq file written to " + outputFileS);
    }
    
    public void writeFasta (String outfileS) {
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            for (int i = 0; i < this.getReadNum(); i++) {
                bw.write(">"+this.reads[i].ID);
                bw.newLine();
                bw.write(this.reads[i].getSequence());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void writeFasta (String outfileS, boolean[] ifOut) {
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            for (int i = 0; i < this.getReadNum(); i++) {
                if(!ifOut[i]) continue;
                bw.write(">"+this.reads[i].ID);
                bw.newLine();
                bw.write(this.reads[i].getSequence());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    /**
     * Return Fastq read
     * @param index
     * @return 
     */
    public Read getRead (int index) {
        return reads[index];
    }
    
    /**
     * Return number of Fastq read
     * @return 
     */
    public int getReadNum () {
        if (reads == null) return 0;
        return reads.length;
    }
}
