package pgl.app.fastCall2;

import pgl.infra.utils.IOUtils;

import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.util.ArrayList;

public class IndividualGenotype implements Comparable<IndividualGenotype> {
    String taxonName = null;
    short chrom = Short.MIN_VALUE;
    int binStart = Integer.MIN_VALUE;
    int binEnd = Integer.MIN_VALUE;
    ArrayList<int[]> allelePackList = new ArrayList<>();

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
            int currentInt = 0;
            while ((currentInt = dis.readInt()) != Integer.MAX_VALUE) {
                int [] currentRecord = new int[AllelePackage.getAllelePackSizeFromFirstInt(currentInt)];
                currentRecord[0] = currentInt;
                for (int i = 0; i < currentRecord.length-1; i++) {
                    currentRecord[i+1] = dis.readInt();
                }
                allelePackList.add(currentRecord);
            }
            dis.close();
        }
        catch (Exception e) {
            System.out.println(fileS);
            e.printStackTrace();
        }

    }

    public void writeTextFile (String outfileS) {
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write(this.getTaxonName());
            bw.newLine();
            bw.write("Chromosome: "+String.valueOf(this.chrom));
            bw.newLine();
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < this.getPositionNumber(); i++) {
                AllelePackage ap = new AllelePackage (this.allelePackList.get(i));
                bw.write(ap.getAlleleInfo(this.binStart, sb).toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public String getTaxonName () {
        return this.taxonName;
    }

    public int getPositionNumber () {
        return allelePackList.size();
    }

    public int getAlleleChromPosition(int alleleIndex) {
        return AllelePackage.getAlleleChromPosition(allelePackList.get(alleleIndex), binStart);
    }

    public int[] getAllelePack(int alleleIndex) {
        return this.allelePackList.get(alleleIndex);
    }

    @Override
    public int compareTo(IndividualGenotype o) {
        return taxonName.compareTo(o.taxonName);
    }
}
