package pgl.app.fastCall2;

import htsjdk.samtools.util.IOUtil;
import pgl.infra.utils.IOUtils;

import java.io.DataInputStream;
import java.io.File;

class IndividualGenotype implements Comparable<IndividualGenotype> {
    String taxonName = null;
    short chrom = Short.MIN_VALUE;
    int binStart = Integer.MIN_VALUE;
    int binEnd = Integer.MIN_VALUE;

    public IndividualGenotype (String fileS) {
        this.readFile(fileS);
    }

    private void readFile (String fileS) {
        try {
            DataInputStream dis = IOUtils.getBinaryGzipReader(fileS);
            this.taxonName = dis.readUTF();
            this.chrom = dis.readShort();
            this.binStart = dis.readInt();
            this.binEnd = dis.read();


            dis.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }


    }

    @Override
    public int compareTo(IndividualGenotype o) {
        return taxonName.compareTo(o.taxonName);
    }
}
