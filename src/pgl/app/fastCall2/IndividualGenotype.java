package pgl.app.fastCall2;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import pgl.infra.dna.BaseEncoder;
import pgl.infra.utils.IOUtils;

import java.io.DataInputStream;
import java.util.ArrayList;
import java.util.Collections;

class IndividualGenotype implements Comparable<IndividualGenotype> {
    String taxonName = null;
    short chrom = Short.MIN_VALUE;
    int binStart = Integer.MIN_VALUE;
    int binEnd = Integer.MIN_VALUE;
    IntArrayList codedAllelePackList = new IntArrayList();
    IntArrayList indelIndexList = new IntArrayList();
    ArrayList<long[]> indelSeqLList = new ArrayList<>();

    public IndividualGenotype (String fileS) {
        this.readFile(fileS);
    }

    private void readFile (String fileS) {
        try {
            DataInputStream dis = IOUtils.getBinaryGzipReader(fileS);
            this.taxonName = dis.readUTF();
            this.chrom = dis.readShort();
            this.binStart = dis.readInt();
            this.binEnd = dis.readInt();
            int currentRecord = 0;
            int counter = 0;
            int indelLength = 0;
            long[] indelSeqL = null;
            while ((currentRecord = dis.readInt()) != Integer.MIN_VALUE) {
                codedAllelePackList.add(currentRecord);
                indelLength = FastCall2.getIndelLength(currentRecord);
                if (indelLength != 0) {
                    if (indelLength < BaseEncoder.longChunkSize) {
                        indelSeqL = new long[1];
                        indelSeqL[0] = dis.readLong();
                    }
                    else {
                        indelSeqL = new long[2];
                        for (int i = 0; i < indelSeqL.length; i++) {
                            indelSeqL[i] = dis.readLong();
                        }
                    }
                    indelIndexList.add(counter);
                    indelSeqLList.add(indelSeqL);
                }
                counter++;
            }
            dis.close();
        }
        catch (Exception e) {
            System.out.println(fileS);
            e.printStackTrace();
        }

    }

    public String getTaxonName () {
        return this.taxonName;
    }

    public int getPositionNumber () {
        return codedAllelePackList.size();
    }

    public int getAllelePosition (int alleleIndex) {
        return FastCall2.getAllelePosition(codedAllelePackList.getInt(alleleIndex), binStart);
    }

    public int getAlleleCoding(int alleleIndex) {
        return FastCall2.getAlleleCoding(codedAllelePackList.getInt(alleleIndex));
    }

    public char getAlleleBase (int alleleIndex) {
        return FastCall2.getAlleleBase(codedAllelePackList.getInt(alleleIndex));
    }

    public byte getIndelLength (int alleleIndex) {
        return FastCall2.getIndelLength(codedAllelePackList.getInt(alleleIndex));
    }

    public int getCodedAllelePack(int alleleIndex) {
        return this.codedAllelePackList.getInt(alleleIndex);
    }

    public long[] getIndelSeqL (int alleleIndex) {
        int index = Collections.binarySearch(indelIndexList, alleleIndex);
        return indelSeqLList.get(index);
    }

    @Override
    public int compareTo(IndividualGenotype o) {
        return taxonName.compareTo(o.taxonName);
    }
}
