package pgl.app.fastCall3;

import pgl.infra.utils.IOUtils;

import java.io.BufferedWriter;
import java.io.DataInputStream;

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
            System.out.println(infileS);
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
}
