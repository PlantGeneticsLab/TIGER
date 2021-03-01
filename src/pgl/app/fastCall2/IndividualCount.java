package pgl.app.fastCall2;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import pgl.infra.utils.IOUtils;

import java.io.BufferedReader;
import java.io.DataInputStream;

class IndividualCount implements Comparable<IndividualCount> {
    String taxonName = null;
    short chrom = Short.MIN_VALUE;
    int binStart = Integer.MIN_VALUE;
    int binEnd = Integer.MIN_VALUE;
    byte[] alleleNum = null;
    short[][] alleleCounts = null;

    public IndividualCount (String infileS) {
        this.readBinaryFileS(infileS);
    }

    private void readBinaryFileS (String infileS) {
        try {
            DataInputStream dis = IOUtils.getBinaryGzipReader(infileS);
            this.taxonName = dis.readUTF();
            this.chrom = dis.readShort();
            this.binStart = dis.readInt();
            this.binEnd = dis.readInt();
            int positionNum = dis.readInt();
            alleleNum = new byte[positionNum];
            alleleCounts = new short[positionNum][];
            for (int i = 0; i < positionNum; i++) {
                alleleNum[i] = dis.readByte();
                if (alleleNum[i] < 0) continue;
                alleleCounts[i] = new short[alleleNum[i]];
                for (int j = 0; j < alleleNum[i]; j++) {
                    alleleCounts[i][j] = dis.readShort();
                }
            }
            dis.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    @Override
    public int compareTo(IndividualCount o) {
        return this.taxonName.compareTo(o.taxonName);
    }
}
