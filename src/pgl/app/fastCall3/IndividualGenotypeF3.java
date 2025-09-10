package pgl.app.fastCall3;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import pgl.infra.utils.IOUtils;
import java.io.BufferedWriter;
import java.io.DataInputStream;

/**
 * A class for storing and managing genotype information for individual samples in the FastCall3 pipeline.
 * <p>
 * This class implements {@link Comparable} to enable sorting of individual genotypes by taxon name.
 * It provides functionality to read binary formatted genotype data and export it to a human-readable text format.
 *
 * <p>Key features:
 * <ul>
 *   <li>Stores genotype calls for a single individual across multiple genomic positions</li>
 *   <li>Uses compressed allele pack format for efficient storage</li>
 *   <li>Supports reading from binary gzipped files</li>
 *   <li>Provides methods for text-based output of genotype data</li>
 *   <li>Organizes data by genomic bins for efficient processing</li>
 * </ul>
 *
 * <p>The binary file format consists of:
 * <ol>
 *   <li>Taxon name (UTF string)</li>
 *   <li>Chromosome number (short)</li>
 *   <li>Bin start position (int)</li>
 *   <li>Bin end position (int)</li>
 *   <li>Sequence of allele packs (int[]), terminated by Integer.MAX_VALUE</li>
 * </ol>
 *
 * @author Fei Lu
 * @version 3.0
 * @since 1.0
 * @see Comparable
 * @see AllelePackageF3
 */
public class IndividualGenotypeF3 implements Comparable<IndividualGenotypeF3> {
    String taxonName = null;
    short chrom = Short.MIN_VALUE;
    int binStart = Integer.MIN_VALUE;
    int binEnd = Integer.MIN_VALUE;
    IntArrayList allelePackList = new IntArrayList();

    public IndividualGenotypeF3(String fileS) {
        this.readFile(fileS);
    }

    /**
     * Reads a binary file and stores the data in this object.
     * The file should contain taxon name (UTF string), chromosome number (short), bin start position (int), bin end position (int), and a sequence of allele packs (int[]), terminated by Integer.MAX_VALUE.
     * @param fileS the path to the binary file
     * @since 1.0
     */
    private void readFile (String fileS) {
        try {
            DataInputStream dis = IOUtils.getBinaryGzipReader(fileS);
            this.taxonName = dis.readUTF();
            this.chrom = dis.readShort();
            this.binStart = dis.readInt();
            this.binEnd = dis.readInt();
            int currentInt = 0;
            while ((currentInt = dis.readInt()) != Integer.MAX_VALUE) {
                allelePackList.add(currentInt);
            }
            dis.close();
        }
        catch (Exception e) {
            System.out.println(fileS);
            e.printStackTrace();
        }

    }

    /**
     * Writes the genotype data in a human-readable text format to a file.
     * The output file will contain the taxon name, chromosome number, and a list of alleles at each position, with each allele represented as a tab-delimited string of: position, allele number, allele count.
     * @param outfileS the path to the output file
     * @since 1.0
     */
    public void writeTextFile (String outfileS) {
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write(this.getTaxonName());
            bw.newLine();
            bw.write("Chromosome: "+String.valueOf(this.chrom));
            bw.newLine();
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < this.getPositionNumber(); i++) {
                bw.write(AllelePackageF3.getAlleleInfo(this.allelePackList.getInt(i), this.binStart, sb).toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * Gets the taxon name for this individual.
     * @return the taxon name
     * @since 1.0
     */
    public String getTaxonName () {
        return this.taxonName;
    }

    /**
     * Returns the number of positions in the genomic bin associated with this individual.
     * @return the number of positions
     * @since 1.0
     */
    public int getPositionNumber () {
        return allelePackList.size();
    }

    /**
     * Gets the chromosomal position of the allele at the given index.
     * The position is given as a zero-based index relative to the start of the genomic bin.
     * @param alleleIndex the index of the allele
     * @return the chromosomal position of the allele
     * @since 1.0
     */
    public int getAlleleChromPosition(int alleleIndex) {
        return AllelePackageF3.getAlleleChromPosition(allelePackList.get(alleleIndex), binStart);
    }

    /**
     * Returns the allele pack integer at the given index.
     * The allele pack integer can be decoded into individual alleles using the methods in {@link AllelePackageF3}.
     * @param alleleIndex the index of the allele pack
     * @return the allele pack integer
     * @since 1.0
     */
    public int getAllelePack(int alleleIndex) {
        return this.allelePackList.getInt(alleleIndex);
    }

    /**
     * Compares this object with the specified object for order.  Returns a
     * negative integer, zero, or a positive integer as this object is less
     * than, equal to, or greater than the specified object.
     * <p>
     * The comparison is based on the taxon name.
     * @param o the object to be compared.
     * @return a negative integer, zero, or a positive integer as this object
     *         is less than, equal to, or greater than the specified object.
     * @since 1.0
     */
    @Override
    public int compareTo(IndividualGenotypeF3 o) {
        return taxonName.compareTo(o.taxonName);
    }
}
