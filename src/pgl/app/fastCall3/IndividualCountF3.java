package pgl.app.fastCall3;

import pgl.infra.utils.IOUtils;

import java.io.BufferedWriter;
import java.io.Closeable;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.IOException;

/**
 * A class for storing and managing allele count information for individual samples in the FastCall3 pipeline.
 * <p>
 * This class implements {@link Comparable} to enable sorting of individual counts by taxon name.
 * It provides functionality to read binary formatted allele count data and export it to a human-readable text format.
 *
 * <p>Key features:
 * <ul>
 *   <li>Stores allele counts for a single individual across multiple genomic positions</li>
 *   <li>Supports reading from binary gzipped files</li>
 *   <li>Provides methods for text-based output of allele count data</li>
 *   <li>Organizes data by genomic bins for efficient processing</li>
 * </ul>
 *
 * <p>The binary file format consists of:
 * <ol>
 *   <li>Taxon name (UTF string)</li>
 *   <li>Chromosome number (short)</li>
 *   <li>Bin start position (int)</li>
 *   <li>Bin end position (int)</li>
 *   <li>Number of positions (int)</li>
 *   <li>For each position:
 *     <ul>
 *       <li>Number of alleles (byte, negative if data is missing)</li>
 *       <li>Counts for each allele (short[])</li>
 *     </ul>
 *   </li>
 * </ol>
 *
 * @author Fei Lu
 * @version 3.0
 * @since 1.0
 * @see Comparable
 */
class IndividualCountF3 implements Comparable<IndividualCountF3> {
    String taxonName = null;
    short chrom = Short.MIN_VALUE;
    int binStart = Integer.MIN_VALUE;
    int binEnd = Integer.MIN_VALUE;
    byte[] alleleNum = null;
    //set null if the site is missing
    short[][] alleleCounts = null;

    public IndividualCountF3(String infileS) {
        this.readBinaryFileS(infileS);
    }

    IndividualCountF3(String taxonName, short chrom, int binStart, int binEnd, byte[] alleleNum, short[][] alleleCounts) {
        this.taxonName = taxonName;
        this.chrom = chrom;
        this.binStart = binStart;
        this.binEnd = binEnd;
        this.alleleNum = alleleNum;
        this.alleleCounts = alleleCounts;
    }

    void writeBinaryFileS(String outfileS) {
        try {
            DataOutputStream dos = IOUtils.getBinaryGzipWriter(outfileS);
            dos.writeUTF(this.taxonName);
            dos.writeShort(this.chrom);
            dos.writeInt(this.binStart);
            dos.writeInt(this.binEnd);
            dos.writeInt(this.alleleNum.length);
            for (int i = 0; i < this.alleleNum.length; i++) {
                if (this.alleleCounts[i] == null) {
                    dos.writeByte((byte) -1);
                    continue;
                }
                dos.writeByte((byte) this.alleleCounts[i].length);
                for (int j = 0; j < this.alleleCounts[i].length; j++) {
                    dos.writeShort(this.alleleCounts[i][j]);
                }
            }
            dos.flush();
            dos.close();
        }
        catch (Exception e) {
            System.err.println(outfileS);
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * Read a binary file storing the counts of each allele for each site of a single individual.
     * The file should be in the format of gzip compressed binary file.
     * The content of the file is as follows:
     * 1. taxonName (String)
     * 2. chrom (short)
     * 3. binStart (int)
     * 4. binEnd (int)
     * 5. positionNum (int)
     * 6. For each position:
     *      alleleNum (byte): the number of alleles, negative if missing
     *      alleleCounts (short[]): the counts of each allele
     * @param infileS the path of the input file
     */
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
            System.err.println(infileS);
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * Write the counts of each allele for each site of this individual to a text file.
     * The content of the file is as follows:
     * 1. taxonName
     * 2. Chromosome: chrom
     * 3. BinStart: binStart
     * 4. BinEnd: binEnd
     * 5. PositionNum: positionNum
     * 6. For each position:
     *      position\talleleNum\talleleCounts
     * @param outfileS the path of the output file
     */
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

    /**
     * Construct a string as follows:
     * alleleNum: alleleCount1\talleleCount2\t...
     * @param sb the StringBuilder to be used
     * @param index the index in alleleNum and alleleCounts
     * @return the constructed string
     */
    private String getAlleleCountInfo (StringBuilder sb, int index) {
        sb.append(this.alleleNum[index]).append(":");
        for (int i = 0; i < this.alleleNum[index]; i++) {
            sb.append("\t").append(String.valueOf(this.alleleCounts[index][i]));
        }
        return sb.toString();
    }

    /**
     * Compare two {@code IndividualCountF3} objects.
     * @param o the other object
     * @return a negative integer, zero, or a positive integer as this object
     *         is less than, equal to, or greater than the specified object.
     */
    @Override
    public int compareTo(IndividualCountF3 o) {
        return this.taxonName.compareTo(o.taxonName);
    }

    /**
     * Sequential reader for one {@code .iac.gz}. Sites are yielded in library order
     * so a bin can be genotyped without loading every taxon into memory.
     */
    static final class SiteStream implements Closeable {
        final String taxonName;
        final short chrom;
        final int binStart;
        final int binEnd;
        final int positionNum;
        private final DataInputStream dis;
        private int remaining;

        private SiteStream(String taxonName, short chrom, int binStart, int binEnd,
                           int positionNum, DataInputStream dis) {
            this.taxonName = taxonName;
            this.chrom = chrom;
            this.binStart = binStart;
            this.binEnd = binEnd;
            this.positionNum = positionNum;
            this.dis = dis;
            this.remaining = positionNum;
        }

        static SiteStream open(String infileS) throws IOException {
            File f = new File(infileS);
            if (!f.exists()) {
                return null;
            }
            DataInputStream dis = IOUtils.getBinaryGzipReader(infileS);
            if (dis == null) {
                throw new IOException("Cannot open individual count file: " + infileS);
            }
            String taxonName = dis.readUTF();
            short chrom = dis.readShort();
            int binStart = dis.readInt();
            int binEnd = dis.readInt();
            int positionNum = dis.readInt();
            return new SiteStream(taxonName, chrom, binStart, binEnd, positionNum, dis);
        }

        short[] readSite() throws IOException {
            if (this.remaining <= 0) {
                throw new IOException("No remaining sites in individual count stream for " + this.taxonName);
            }
            this.remaining--;
            byte n = this.dis.readByte();
            if (n < 0) {
                return null;
            }
            short[] counts = new short[n];
            for (int i = 0; i < n; i++) {
                counts[i] = this.dis.readShort();
            }
            return counts;
        }

        @Override
        public void close() throws IOException {
            this.dis.close();
        }
    }
}
