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
 * Holding FastA format sequence, providing functions of sorting, searching and collecting statistics
 * Representing full IUPAC codes. (https://www.bioinformatics.org/sms/iupac.html).
 * Higher speed, more memory cost than {@link pgl.infra.dna.FastaBit}.
 * @author feilu
 */
public class FastaByte extends FastaAbstract {
    /**
     * Constructs a {@link pgl.infra.dna.FastaByte} from input file. The file should be either txt format or gz format.
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
     * Constructs a {@link pgl.infra.dna.FastaByte} from input file.
     * @param infileS
     * @param format 
     */
    public FastaByte (String infileS, IOFileFormat format) {
        this.readFasta(infileS, format);
    }
    
    /**
     * Constructs a {@link pgl.infra.dna.FastaByte}
     * @param names
     * @param seqs
     * @param ids 
     */
    public FastaByte (String[] names, String[] seqs, int[] ids) {
        records = new FastaRecordByte[names.length];
        for (int i = 0; i < records.length; i++) {
            records[i] = new FastaRecordByte(names[i], seqs[i], ids[i]);
        }
    }
    
    /**
     * Constructs a new {@link pgl.infra.dna.FastaByte} from a list of them
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
     * Construct an object from a FastaRecord
     * @param f
     */
    public FastaByte (FastaRecordByte f) {
        records = new FastaRecordByte[1];
        records[0] = f;
        records[0].setID(1);
        sType = sortType.byID;
    }

    /**
     * Construct an object from a group of FastaRecord
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
            String temp = null, name = null, seq = null;
            StringBuilder sb = new StringBuilder();
            FastaRecordByte fr;
            boolean first = true;
            int cnt = 1;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith(">")) {
                    if (first == false) {
                        seq = sb.toString();
                        fr = new FastaRecordByte(name, seq, cnt);
                        fl.add(fr);
                        sb = new StringBuilder();
                        if (cnt%1000000 == 0) {
                            System.out.println("Read "+String.valueOf(cnt)+" sequences");
                        }
                        cnt++;
                    }
                    name = temp.substring(1, temp.length());
                    first = false;
                }
                else {
                    sb.append(temp);
                }
            }
            if (!name.equals("")) {
                seq = sb.toString();
                fr = new FastaRecordByte(name, seq, cnt);
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
    public int getIndexByName (String name) {
        if (this.sType != sortType.byName) {
            this.sortByName();
        }
        return Arrays.binarySearch(records, new FastaRecordByte(name,"A",-1), new sortByName());
    }
     
    @Override
    public boolean isThereNonACGTNBase () {
        boolean value;
        for (int i = 0; i < records.length; i++) {
             value = records[i].isThereNonACGTNBase();
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
        FastaRecordByte frb = new FastaRecordByte(this.getName(index), sb, this.records[index].getID());
        return frb;
    }
}
