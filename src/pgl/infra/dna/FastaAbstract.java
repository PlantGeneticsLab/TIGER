/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.dna;

import java.io.BufferedWriter;
import java.io.File;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;

/**
 * Providing functions of sorting, searching and collecting statistics.
 * @author feilu
 */
public abstract class FastaAbstract implements FastaInterface {
    FastaRecordInterface[] records = null;
    protected enum sortType {byName, byNameValue, byID, byLengthAscending, byLengthDescending}
    sortType sType = null;
    
    @Override
    public void writeFasta (String outfileS, boolean[] ifOut, IOFileFormat format) {
        int cnt = 0;
        try {
            BufferedWriter bw = null;
            if (format == IOFileFormat.Text) {
                bw = IOUtils.getTextWriter(outfileS);
            }
            else if (format == IOFileFormat.TextGzip) {
                bw = IOUtils.getTextGzipWriter(outfileS);
            }
            else {
                throw new UnsupportedOperationException("Invalid operation for output");
            }
            for (int i = 0; i < records.length; i++) {
                if (!ifOut[i]) continue;
                bw.write(">"+records[i].getDescription());
                bw.newLine();
                bw.write(PStringUtils.getMultiplelineString(60, records[i].getSequence()));
                bw.newLine();
                cnt++;
            }
            bw.flush();
            bw.close();
            System.out.println(cnt+ " sequences are written in " + outfileS);
        }
        catch (Exception e) {
            System.out.println("Error while writing "+ outfileS);
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    @Override
    public void writeFasta (String outfileS, int index, IOFileFormat format) {
        try {
            BufferedWriter bw = null;
            if (format == IOFileFormat.Text) {
                bw = IOUtils.getTextWriter(outfileS);
            }
            else if (format == IOFileFormat.TextGzip) {
                bw = IOUtils.getTextGzipWriter(outfileS);
            }
            else {
                throw new UnsupportedOperationException("Invalid operation for output");
            }
            bw.write(">"+records[index].getDescription());
            bw.newLine();
            bw.write(PStringUtils.getMultiplelineString(60, records[index].getSequence()));
            bw.newLine();
            bw.flush();
            bw.close();
            System.out.println("No." + index +" sequence is written in " + outfileS);
        }
        catch (Exception e) {
            System.out.println("Error while writing "+ outfileS);
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    @Override
    public void writeFasta (String outfileS, IOFileFormat format) {
        try {
            BufferedWriter bw = null;
            if (format == IOFileFormat.Text) {
                bw = IOUtils.getTextWriter(outfileS);
            }
            else if (format == IOFileFormat.TextGzip) {
                bw = IOUtils.getTextGzipWriter(outfileS);
            }
            else {
                throw new UnsupportedOperationException("Invalid operation for output");
            }
            for (int i = 0; i < records.length; i++) {
                bw.write(">"+records[i].getDescription());
                bw.newLine();
                bw.write(PStringUtils.getMultiplelineString(60, records[i].getSequence()));
                bw.newLine();
            }
            bw.flush();
            bw.close();
            System.out.println(records.length+ " sequences are written in " + outfileS);
        }
        catch (Exception e) {
            System.out.println("Error while writing "+ outfileS);
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    /**
     * Write fast file by chromosome, designed for large genome files.
     * @param outfileDirS
     * @param format
     * @param prefix
     * @param postfix 
     */
    public void writeFastaByChr (String outfileDirS, IOFileFormat format, String prefix, String postfix) {
        List<FastaRecordInterface> rList = Arrays.asList(records);
        rList.parallelStream().forEach(r -> {
            try {
                String outfileS = new File (outfileDirS, (prefix+r.getDescription()+postfix)).getAbsolutePath();
                BufferedWriter bw = null;
                if (format == IOFileFormat.Text) {
                    bw = IOUtils.getTextWriter(outfileS);
                }
                else if (format == IOFileFormat.TextGzip) {
                    bw = IOUtils.getTextGzipWriter(outfileS);
                }
                else {
                    throw new UnsupportedOperationException("Invalid operation for output");
                }
                bw.write(">"+r.getDescription());
                bw.newLine();
                bw.write(PStringUtils.getMultiplelineString(60, r.getSequence()));
                bw.newLine();
                bw.flush();
                bw.close();
                System.out.println(r.getDescription() +" sequence is written in " + outfileS);
            }
            catch (Exception e) {
                System.out.println("Error while writing "+ outfileDirS);
                e.printStackTrace();
                System.exit(1);
            }
        });
        
    }
    
    @Override
    public int getN50 () {
        if (sType != sortType.byLengthDescending) this.sortByLengthDescending();
        long sum = this.getTotalSeqLength();
        long halfSum = sum/2;
        int current = 0;
        for (int i = 0; i < this.getSeqNumber(); i++) {
            current+=this.getSeqLength(i);
            if (current > halfSum) return i+1;
        }
        return -1;
    }
    
    @Override
    public int getL50 () {
        if (sType != sortType.byLengthDescending) this.sortByLengthDescending();
        long sum = this.getTotalSeqLength();
        long halfSum = sum/2;
        int current = 0;
        for (int i = 0; i < this.getSeqNumber(); i++) {
            current+=this.getSeqLength(i);
            if (current > halfSum) return this.getSeqLength(i);
        }
        return -1;
    }
    
    @Override
    public long getTotalSeqLength () {
        long sum = 0;
        for (int i = 0; i < this.getSeqNumber(); i++) {
            sum+=this.getSeqLength(i);
        }
        return sum;
    }
    
    @Override
    public int getSeqNumber () {
        return records.length;
    }
    
    @Override
    public int getSeqLength (int index) {
        return records[index].getSequenceLength();
    }
    
    @Override
    public String[] getDescriptions() {
        String[] names = new String[this.getSeqNumber()];
        for (int i = 0; i < names.length; i++) names[i] = this.getDescription(i);
        return names;
    }
    
    @Override
    public String getDescription(int index) {
        return records[index].getDescription();
    }
    
    @Override
    public String getSeq (int index) {
        return records[index].getSequence();
    }
    
    @Override
    public FastaRecordInterface getFastaRecord(int index) {
        return this.records[index];
    }
    
    @Override
    public String getSeq (int index, int startIndex, int endIndex) {
        return records[index].getSequence(startIndex, endIndex);
    }
    
    @Override
    public void setDescription(String description, int index) {
        records[index].setDescription(description);
    }
    
    @Override
    public void sortByDescription() {
        Arrays.parallelSort(records, new sortByName());
        this.sType = sortType.byName;
    }
    

    @Override
    public void sortByDescriptionValue() {
        Arrays.parallelSort(records, new sortByName());
        this.sType = sortType.byNameValue;
    }
    
    @Override
    public void sortByID () {
        Arrays.parallelSort (records, new sortByID());
        this.sType = sortType.byID;
    }
    
    @Override
    public void sortByLengthAscending () {
        Arrays.parallelSort (records, new sortByLengthAscending());
        this.sType = sortType.byLengthAscending;
    }
    
    @Override
    public void sortByLengthDescending () {
        Arrays.parallelSort (records, new sortByLengthDescending());
        this.sType = sortType.byLengthDescending;
    }
    
    @Override
    public boolean isThereN () {
        boolean value;
        for (int i = 0; i < records.length; i++) {
             value = records[i].isThereN();
             if (value) return true;
        }
        return false;
    }
    
    protected class sortByID implements Comparator <FastaRecordInterface> {
        @Override
        public int compare(FastaRecordInterface o1, FastaRecordInterface o2) {
            return o1.getID() - o2.getID();
        }
    }
    
    protected class sortByName implements Comparator <FastaRecordInterface> {
        @Override
        public int compare (FastaRecordInterface o1, FastaRecordInterface o2) {
            return o1.getDescription().compareTo(o2.getDescription());
        }
    }
    
    protected class sortByNameValue implements Comparator <FastaRecordInterface> {
        @Override
        public int compare (FastaRecordInterface o1, FastaRecordInterface o2) {
            int n1 = Integer.parseInt(o1.getDescription());
            int n2 = Integer.parseInt(o2.getDescription());
            if (n1 < n2) return -1;
            else if (n1 > n2) return 1;
            return 0;
        }
    }
    
    protected class sortByLengthAscending implements Comparator <FastaRecordInterface> {
        @Override
        public int compare (FastaRecordInterface o1, FastaRecordInterface o2) {
            return o1.getSequenceLength() - o2.getSequenceLength();
        }
    }
    
    protected class sortByLengthDescending implements Comparator <FastaRecordInterface> {
        @Override
        public int compare (FastaRecordInterface o1, FastaRecordInterface o2) {
            return o2.getSequenceLength()-o1.getSequenceLength();
        }
    }
}
