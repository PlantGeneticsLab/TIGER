package pgl.app.fastCall3;

import pgl.infra.dna.BaseEncoder;
import pgl.infra.dna.allele.AlleleEncoder;

/**
 * A wrapper class of allele pack, which enabled {@link java.util.Set} functionality of allele pack
 * Using array of int, an allele pack stores allele state, relative position of the allele in bin, length of Indel (if it exists), and sequence of the Indel
 */
class AllelePackageF3 implements Comparable<AllelePackageF3>{
    int[] allelePack = null;

    public AllelePackageF3(int[] allelePacks) {
        this.allelePack = allelePacks;
    }

    public int[] getAllelePack() {
        return this.allelePack;
    }

    public int getFirstIntOfAllelePack() {
        return this.allelePack[0];
    }

    public int getAllelePackSize () {
        return this.allelePack.length;
    }

    public byte getAlleleCoding () {
        return getAlleleCoding(getFirstIntOfAllelePack());
    }

    public char getAlleleBase () {
        return getAlleleBase(getFirstIntOfAllelePack());
    }

    public byte getIndelLength () {
        return getIndelLength(getFirstIntOfAllelePack());
    }

    public String getIndelSeq () {
        if (this.getAllelePackSize() == 1) return "";
        int[] des = new int[this.getAllelePackSize()-1];
        System.arraycopy(this.allelePack, 1, des, 0, des.length);
        return BaseEncoder.getSequenceFromInts(des).substring(0, this.getIndelLength());
    }

    public int getAlleleChromPosition (int binStart) {
        return getAlleleChromPosition (getFirstIntOfAllelePack(), binStart);
    }

    public StringBuilder getAlleleInfo (int binStart, StringBuilder sb) {
        sb.setLength(0);
        sb.append(this.getAlleleChromPosition(binStart)).append("\t").append(this.getAlleleBase());
        if (this.getIndelLength() != 0) {
            sb.append("\t").append(this.getIndelSeq());
        }
        return sb;
    }

    static StringBuilder getAlleleInfo (int allelePack, int binStart, StringBuilder sb) {
        sb.setLength(0);
        sb.append(getAlleleChromPosition(allelePack,binStart)).append("\t").append(getAlleleBase(allelePack));
        int indelLength = getIndelLength(allelePack);
        if (indelLength != 0) {
            sb.append("\t").append(indelLength);
        }
        return sb;
    }

    @Override
    public int compareTo(AllelePackageF3 o) {
        if (allelePack.length != o.allelePack.length) {
            return allelePack.length - o.allelePack.length;
        }
        else {
            for (int i = 0; i < allelePack.length; i++) {
              if (allelePack[i] < o.allelePack[i]) {
                  return -1;
              }
              else if (allelePack[i] > o.allelePack[i]) {
                  return 1;
              }
              else continue;
            }
        }
        return 0;
    }

    public int hashCode() {
        return allelePack[0];
    }

    public boolean equals(Object o) {
        if (o == this) {
            return true;
        }
        if (!(o instanceof AllelePackageF3)) {
            return false;
        }
        AllelePackageF3 other = (AllelePackageF3) o;
        if (allelePack.length != other.allelePack.length) return false;
        for (int i = 0; i < allelePack.length; i++) {
            if (allelePack[i] != other.allelePack[i]) return false;
        }
        return true;
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
     * Return the int size of an allele pack from the first int of the allele pack
     * @param firstIntOfAllelePack
     * @return
     */
    static int getAllelePackSizeFromFirstInt (int firstIntOfAllelePack) {
        int indelLength = getIndelLength(firstIntOfAllelePack);
        return getAllelePackSizeFromIndelLength(indelLength);
    }

    /**
     * Return the chromosomal position of the allele
     * @param allelePack
     * @param binStart
     * @return
     */
    static int getAlleleChromPosition(int[] allelePack, int binStart) {
        return (allelePack[0] >>> 9) + binStart;
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

    static char getAlleleBase (int allelePack) {
        return AlleleEncoder.alleleCodingToBaseMap.get(getAlleleCoding(allelePack));
    }

    static byte getIndelLength (int allelePack) {
        return (byte)(FastCall3.maxIndelLength & allelePack);
    }
}
