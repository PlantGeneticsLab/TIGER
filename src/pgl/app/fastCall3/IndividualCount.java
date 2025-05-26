package pgl.app.fastCall3;

import pgl.infra.utils.IOUtils;

import java.io.BufferedWriter;
import java.io.DataInputStream;

class IndividualCount implements Comparable<IndividualCount> {
    String taxonName = null;
    short chrom = Short.MIN_VALUE;
    int binStart = Integer.MIN_VALUE;
    int binEnd = Integer.MIN_VALUE;
    byte[] alleleNum = null;
    //set null if the site is missing
    short[][] alleleCounts = null;

    public IndividualCount(String infileS) {
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
            System.out.println(infileS);
            e.printStackTrace();
            System.exit(1);
        }
    }

    public void writeTextFile (String outfileS) {
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write(this.taxonName);
            bw.newLine();
            bw.write("Chromosome: "+String.valueOf(this.chrom));
            bw.newLine();
            bw.write("BinStart: "+String.valueOf(this.binStart));
            bw.newLine();
            bw.write("BinEnd: "+String.valueOf(this.binEnd));
            bw.newLine();
            bw.write("PositionNum: "+String.valueOf(this.alleleNum.length));
            bw.newLine();
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < this.alleleNum.length; i++) {
                bw.write(this.getAlleleCountInfo(sb, i));
                bw.newLine();
                sb.setLength(0);
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    private String getAlleleCountInfo (StringBuilder sb, int index) {
        sb.append(this.alleleNum[index]).append(":");
        for (int i = 0; i < this.alleleNum[index]; i++) {
            sb.append("\t").append(String.valueOf(this.alleleCounts[index][i]));
        }
        return sb.toString();
    }

    @Override
    public int compareTo(IndividualCount o) {
        return this.taxonName.compareTo(o.taxonName);
    }
}
