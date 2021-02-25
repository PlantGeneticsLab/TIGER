package pgl.app.fastCall2;

import htsjdk.samtools.util.IOUtil;
import it.unimi.dsi.fastutil.bytes.ByteArrayList;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import pgl.infra.dna.allele.AlleleEncoder;
import pgl.infra.utils.IOUtils;

import java.io.DataInputStream;
import java.io.File;

class IndividualGenotype implements Comparable<IndividualGenotype> {
    String taxonName = null;
    short chrom = Short.MIN_VALUE;
    int binStart = Integer.MIN_VALUE;
    int binEnd = Integer.MIN_VALUE;
    IntArrayList codedAlleleInfo = null;

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
            codedAlleleInfo = new IntArrayList();
            int currentRecord = 0;
            while ((currentRecord = dis.readInt()) != Integer.MIN_VALUE) {
                codedAlleleInfo.add(currentRecord);
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
        return codedAlleleInfo.size();
    }

    public int getAllelePosition (int alleleIndex) {
        return FastCall2.getAllelePosition(codedAlleleInfo.getInt(alleleIndex), binStart);
    }

    public int getAlleleByte (int alleleIndex) {
        return FastCall2.getAlleleByte(codedAlleleInfo.getInt(alleleIndex));
    }

    public char getAlleleBase (int alleleIndex) {
        return FastCall2.getAlleleBase(codedAlleleInfo.getInt(alleleIndex));
    }

    public byte getIndelLength (int alleleIndex) {
        return FastCall2.getIndelLength(codedAlleleInfo.getInt(alleleIndex));
    }

    public int getCodedAlleleInfo (int alleleIndex) {
        return this.codedAlleleInfo.getInt(alleleIndex);
    }

    @Override
    public int compareTo(IndividualGenotype o) {
        return taxonName.compareTo(o.taxonName);
    }
}
