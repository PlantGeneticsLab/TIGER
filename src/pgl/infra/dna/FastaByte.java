/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.dna;

import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;

/**
 * Holding Fasta format sequence, providing functions of sorting, searching and collecting statistics,
 * representing full IUPAC codes. (https://www.bioinformatics.org/sms/iupac.html).
 * <p>
 * Higher speed, more memory cost than {@link FastaBit}.
 * @author feilu
 */
public class FastaByte extends FastaAbstract {
    /**
     * Constructs a {@link FastaByte} from input file, the file should be either txt format or txt.gz format.
     * @param infileS 
     */
    public FastaByte (String infileS) {
        if (infileS.endsWith(".gz")) {
            this.readFasta(infileS, IOFileFormat.TextGzip);
        }
        else {
            this.readFasta(infileS, IOFileFormat.Text);
        }
    }
    
    /**
     * Constructs a {@link FastaByte} from input file, the file should be either txt format or txt.gz format.
     * @param infileS
     * @param format 
     */
    public FastaByte (String infileS, IOFileFormat format) {
        this.readFasta(infileS, format);
    }
    
    /**
     * Constructs a {@link FastaByte}.
     * @param descriptions
     * @param seqs
     * @param ids 
     */
    public FastaByte (String[] descriptions, String[] seqs, int[] ids) {
        records = new FastaRecordByte[descriptions.length];
        for (int i = 0; i < records.length; i++) {
            records[i] = new FastaRecordByte(descriptions[i], seqs[i], ids[i]);
        }
    }
    
    /**
     * Constructs a new {@link FastaByte} from a list of {@link FastaByte}.
     * By rebuilding the references
     * @param fArray 
     */
    public FastaByte (FastaByte[] fArray) {
        int size = 0;
        for (int i = 0; i < fArray.length; i++) {
            size += fArray[i].getSeqNumber();
        }
        records = new FastaRecordByte[size];
        int cnt = 0;
        for (int i = 0; i < fArray.length; i++) {
            for (int j = 0; j < fArray[i].getSeqNumber(); j++) {
                records[cnt] = fArray[i].records[j];
                cnt++;
            }
        }
        for (int i = 0; i < records.length; i++) {
            records[i].setID(i+1);
        }
        sType = sortType.byID;
    }

    /**
     * Construct an object from a {@link FastaByte}.
     * @param f
     */
    public FastaByte (FastaRecordByte f) {
        records = new FastaRecordByte[1];
        records[0] = f;
        records[0].setID(1);
        sType = sortType.byID;
    }

    /**
     * Construct an object from an array of {@link FastaByte}
     * @param fs
     */
    public FastaByte (FastaRecordByte[] fs) {
        records = new FastaRecordByte[fs.length];
        for (int i = 0; i < fs.length; i++) {
            records[i] = fs[i];
            records[i].setID(i+1);
        }
        sType = sortType.byID;
    }
    
    @Override
    public void readFasta (String infileS, IOFileFormat format) {
        System.out.println("Reading Fasta file...");
        List<FastaRecordByte> fl = new ArrayList<>();
        try {
            BufferedReader br = null;
            if (format == IOFileFormat.Text) {
                br = IOUtils.getTextReader(infileS);
            }
            else if (format == IOFileFormat.TextGzip) {
                br = IOUtils.getTextGzipReader(infileS);
            }
            else {
                throw new UnsupportedOperationException("Invalid input format for the Fasta file");
            }
            String temp = null, description = null, seq = null;
            StringBuilder sb = new StringBuilder();
            FastaRecordByte fr;
            boolean first = true;
            int cnt = 1;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith(">")) {
                    if (first == false) {
                        seq = sb.toString();
                        fr = new FastaRecordByte(description, seq, cnt);
                        fl.add(fr);
                        sb = new StringBuilder();
                        if (cnt%1000000 == 0) {
                            System.out.println("Read "+String.valueOf(cnt)+" sequences");
                        }
                        cnt++;
                    }
                    description = temp.substring(1, temp.length());
                    first = false;
                }
                else {
                    sb.append(temp);
                }
            }
            if (!description.equals("")) {
                seq = sb.toString();
                fr = new FastaRecordByte(description, seq, cnt);
                fl.add(fr);
            }
            records = fl.toArray(new FastaRecordByte[fl.size()]);
            sType = sortType.byID;
            System.out.println(records.length + " sequences in the file " + infileS);
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    @Override
    public int getIndexByDescription(String name) {
        if (this.sType != sortType.byName) {
            this.sortByDescription();
        }
        return Arrays.binarySearch(records, new FastaRecordByte(name,"A",-1), new sortByDescription());
    }

    /**
     * Return if the fasta has non-"ACGTN" base
     * @return
     */
    public boolean isThereNonACGTNBase () {
        boolean value;
        for (int i = 0; i < records.length; i++) {
             value = ((FastaRecordByte)records[i]).isThereNonACGTNBase();
             if (value) return true;
        }
        return false;
    }

    /**
     * Return a FastaRecord
     * @param index
     * @return
     */
    public FastaRecordByte getFastaRecordByte (int index) {
        return (FastaRecordByte)this.records[index];
    }

    public FastaRecordByte getFastaRecordByte (int index, int startPositionIndex, int endPositionIndex) {
        SequenceByte sb = (SequenceByte)this.records[index].getSequenceInterface(startPositionIndex, endPositionIndex);
        FastaRecordByte frb = new FastaRecordByte(this.getDescription(index), sb, this.records[index].getID());
        return frb;
    }
}
