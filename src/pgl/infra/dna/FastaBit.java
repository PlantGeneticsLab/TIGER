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
 * Holding FastA format sequence, providing functions of sorting, searching and collecting statistics.
 * Bases are compressed to 3 bits, representing 'A', 'C', 'G', 'T', 'N'. Specially designed for large Fasta sequence.
 * @author feilu
 */
public class FastaBit extends FastaAbstract {
    /**
     * Constructs a {@link pgl.infra.dna.FastaBit} from input file. The file should be either txt format or gz format.
     * @param infileS 
     */
    public FastaBit (String infileS) {
        if (infileS.endsWith(".gz")) {
            this.readFasta(infileS, IOFileFormat.TextGzip);
        }
        else {
            this.readFasta(infileS, IOFileFormat.Text);
        }
    }
    
    /**
     * Constructs a {@link pgl.infra.dna.FastaBit} from input file.
     * @param infileS
     * @param format 
     */
    public FastaBit (String infileS, IOFileFormat format) {
        this.readFasta(infileS, format);
    }
    
    /**
     * Constructs a {@link pgl.infra.dna.FastaBit}
     * @param names
     * @param seqs
     * @param ids 
     */
    public FastaBit (String[] names, String[] seqs, int[] ids) {
        records = new FastaRecordBit[names.length];
        for (int i = 0; i < records.length; i++) {
            records[i] = new FastaRecordBit(names[i], seqs[i], ids[i]);
        }
    }
    
    /**
     * Constructs a new {@link pgl.infra.dna.FastaBit} from a list of them
     * By rebuilding the references
     * @param fArray 
     */
    public FastaBit (FastaBit[] fArray) {
        int size = 0;
        for (int i = 0; i < fArray.length; i++) {
            size += fArray[i].getSeqNumber();
        }
        records = new FastaRecordBit[size];
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
    public FastaBit (FastaRecordBit f) {
        records = new FastaRecordBit[1];
        records[0] = f;
        records[0].setID(1);
        sType = sortType.byID;
    }

    /**
     * Construct an object from a group of FastaRecord
     * @param fs
     */
    public FastaBit (FastaRecordBit[] fs) {
        records = new FastaRecordBit[fs.length];
        for (int i = 0; i < fs.length; i++) {
            records[i] = fs[i];
            records[i].setID(i+1);
        }
        sType = sortType.byID;
    }
    
    @Override
    public void readFasta (String infileS, IOFileFormat format) {
        System.out.println("Reading Fasta file...");
        List<FastaRecordBit> fl = new ArrayList<>();
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
            FastaRecordBit fr;
            boolean first = true;
            int cnt = 1;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith(">")) {
                    if (first == false) {
                        seq = sb.toString();
                        fr = new FastaRecordBit(name, seq, cnt);
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
                fr = new FastaRecordBit(name, seq, cnt);
                fl.add(fr);
            }
            records = fl.toArray(new FastaRecordBit[fl.size()]);
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
        return Arrays.binarySearch(records, new FastaRecordBit(name,"A",-1), new sortByName());
    }
    
    @Override
    public boolean isThereNonACGTNBase() {
        return false;
    }

    /**
     * Return an FastaRecord
     * @param index
     * @return
     */
    public FastaRecordBit getFastaRecordBit (int index) {
        return (FastaRecordBit)this.records[index];
    }

    public FastaRecordBit getFastaRecordBit (int index, int startPositionIndex, int endPositionIndex) {
        Sequence3Bit b3 = (Sequence3Bit)this.records[index].getSequenceInterface(startPositionIndex, endPositionIndex);
        FastaRecordBit frb = new FastaRecordBit(this.getName(index), b3, this.records[index].getID());
        return frb;
    }
}
