package pgl.app.fastCall2;

import htsjdk.samtools.util.IOUtil;
import it.unimi.dsi.fastutil.bytes.ByteArrayList;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import pgl.infra.utils.IOUtils;

import java.io.DataInputStream;
import java.io.File;

class IndividualGenotype implements Comparable<IndividualGenotype> {
    String taxonName = null;
    short chrom = Short.MIN_VALUE;
    int binStart = Integer.MIN_VALUE;
    int binEnd = Integer.MIN_VALUE;
    int[] positions = null;
    byte[] codedAlleles = null;

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
            int currentPosition = 0;
            byte codedAllele = 0;
            int siteDepth = 0;
            int alleleDepth = 0;
            IntArrayList positionList = new IntArrayList();
            ByteArrayList alleleList = new ByteArrayList();
            while ((currentPosition = dis.readInt()) != Integer.MIN_VALUE) {
                positionList.add(currentPosition);
                alleleList.add(dis.readByte());
            }
            this.positions = positionList.toIntArray();
            this.codedAlleles = alleleList.toByteArray();
            dis.close();
        }
        catch (Exception e) {
            System.out.println(fileS);
            e.printStackTrace();
        }


    }

    @Override
    public int compareTo(IndividualGenotype o) {
        return taxonName.compareTo(o.taxonName);
    }
}
