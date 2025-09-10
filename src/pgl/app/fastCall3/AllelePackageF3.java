package pgl.app.fastCall3;

import pgl.infra.dna.BaseEncoder;
import pgl.infra.dna.allele.AlleleEncoder;

/**
 * A utility class for compressing and manipulating allele information in a packed integer format.
 * <p>
 * This class provides methods to work with allele data that is compressed into a 32-bit integer (pack) with the following structure:
 * - Bits 0-22 (23 bits): Position within a bin (max value: 8,388,607)
 * - Bits 23-25 (3 bits): Allele coding
 * - Bits 26-31 (6 bits): Indel length (max value: 63)
 *
 * @author Fei Lu
 * @version 3.0
 * @since 1.0
 */
class AllelePackageF3 {

    /**
     * Get the allele information given an allele pack, bin start and StringBuilder.
     * The format of the returned string is: chromPosition\talleleBase\tindelLength
     * @param allelePack
     * @param binStart
     * @param sb
     * @return
     */
    static StringBuilder getAlleleInfo (int allelePack, int binStart, StringBuilder sb) {
        sb.setLength(0);
        sb.append(getAlleleChromPosition(allelePack,binStart)).append("\t").append(getAlleleBase(allelePack));
        int indelLength = getIndelLength(allelePack);
        if (indelLength != 0) {
            sb.append("\t").append(indelLength);
        }
        return sb;
    }

    /**
     * Compares allele information to an allele pack
     * the first 23 bits are used to record position in a bin (maximum size = 2^23 - 1 = 8,388,607)
     * the next 3 bits are used to record alleles
     * the last 6 bits are used record the length of indel (max size = 2^6 -1 = 63)
     * Above three components are called a "pack"
     * @param binStart
     * @param position
     * @param alleleCoding
     * @param indelLength
     * @return compressed information of an allele in allele pack format
     */
    static int getAllelePack (int binStart, int position, byte alleleCoding, int indelLength) {
        int v = (position-binStart) << 9;
        v = v + (alleleCoding << 6);
        if (indelLength > FastCall3.maxIndelLength) indelLength = FastCall3.maxIndelLength;
        v = v+indelLength;
        return v;
    }

    /**
     * Return the int size of an allele pack from Indel length
     * @param indelLength
     * @return
     */
    static int getAllelePackSizeFromIndelLength(int indelLength) {
        int remainder = indelLength % BaseEncoder.intChunkSize;
        if (remainder == 0) {
            return indelLength/BaseEncoder.intChunkSize+1;
        }
        else {
            return indelLength/BaseEncoder.intChunkSize+2;
        }
    }


    /**
     * Return the chromosomal position of the allele
     * @param firstIntOfAllelePack
     * @param binStart
     * @return
     */
    static int getAlleleChromPosition(int firstIntOfAllelePack, int binStart) {
        int v = (firstIntOfAllelePack >>> 9) + binStart;
        return v;
    }

    /**
     * Return the allele coding
     * @param allelePack
     * @return
     */
    static byte getAlleleCoding(int allelePack) {
        int v = (511 & allelePack) >>> 6;
        return (byte)v;
    }

    /**
     * Returns the base of an allele given an allele pack
     * @param allelePack a compressed information of an allele
     * @return the base of the allele
     */
    static char getAlleleBase (int allelePack) {
        return AlleleEncoder.alleleCodingToBaseMap.get(getAlleleCoding(allelePack));
    }

    /**
     * Return the length of an indel given an allele pack
     * @param allelePack a compressed information of an allele
     * @return the length of the indel
     */
    static byte getIndelLength (int allelePack) {
        return (byte)(FastCall3.maxIndelLength & allelePack);
    }
}
