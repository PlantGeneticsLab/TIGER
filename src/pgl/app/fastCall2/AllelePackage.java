package pgl.app.fastCall2;

import pgl.infra.dna.BaseEncoder;
import pgl.infra.dna.allele.AlleleEncoder;

/**
 * A wrapper class of allele packs, which enabled {@link java.util.Set} functionality of allele packs
 * Using array of int, an allele pack stores allele state, relative position of the allele in bin, length of Indel (if it exists), and sequence of the Indel
 */
class AllelePackage implements Comparable<AllelePackage>{
    int[] allelePacks = null;

    public AllelePackage (int[] allelePacks) {
        this.allelePacks = allelePacks;
    }

    public int[] getAllelePacks() {
        return this.allelePacks;
    }

    @Override
    public int compareTo(AllelePackage o) {
        return allelePacks[0]-o.allelePacks[0];
    }

    public int hashCode() {
        return allelePacks[0];
    }

    public boolean equals(Object o) {
        if (o == this) {
            return true;
        }
        if (!(o instanceof AllelePackage)) {
            return false;
        }
        AllelePackage other = (AllelePackage) o;
        if (allelePacks.length != other.allelePacks.length) return false;
        for (int i = 0; i < allelePacks.length; i++) {
            if (allelePacks[i] != other.allelePacks[i]) return false;
        }
        return true;
    }

    /**
     * Compares allele information to an allele pack
     * @param binStart
     * @param position
     * @param alleleCoding
     * @param indelLength
     * @param indelSeq
     * @return compressed information of an allele in allele pack format
     */
    static int[] getAllelePack (int binStart, int position, byte alleleCoding, int indelLength, String indelSeq) {
        int v = (position-binStart) << 9;
        v = v + (alleleCoding << 6);
        if (indelLength > FastCall2.maxIndelLength) indelLength = FastCall2.maxIndelLength;
        v = v+indelLength;
        int[] allelePack = new int[getAllelePackSizeFromIndelLength(indelLength)];
        allelePack[0] = v;
        if (allelePack.length ==1) return allelePack;
        byte[] seqCodings = BaseEncoder.convertToBaseCodingArray(indelSeq.getBytes());
        for (int i = 0; i < allelePack.length-1; i++) {
            int remainder = indelLength % BaseEncoder.intChunkSize;
            if (remainder == 0) {
                allelePack[i+1] = BaseEncoder.getIntSeqFromSubByteArray(seqCodings, i*BaseEncoder.intChunkSize, (i+1)*BaseEncoder.intChunkSize);
            }
            else {
                allelePack[i+1] = BaseEncoder.getIntSeqFromSubByteArray(seqCodings, i*BaseEncoder.intChunkSize, i*BaseEncoder.intChunkSize+remainder);
            }

        }
        return allelePack;
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
     * Return allele coding and Indel length
     * @param alleleCoding
     * @param indelLength
     * @return
     */
    static short getAlleleCodingLength(byte alleleCoding, int indelLength) {
        int v = (alleleCoding << 6);
        if (indelLength > FastCall2.maxIndelLength) indelLength = FastCall2.maxIndelLength;
        return (short)(v + indelLength);
    }

    /**
     * Return allele coding and Indel length
     * @param firstIntOfAllelePack
     * @return
     */
    static short getAlleleCodingLength(int firstIntOfAllelePack) {
        return (short)(firstIntOfAllelePack&511);
    }

    /**
     * Return the allele coding
     * @param firstIntOfAllelePack
     * @return
     */
    static byte getAlleleCoding(int firstIntOfAllelePack) {
        int v = (511 & firstIntOfAllelePack) >>> 6;
        return (byte)v;
    }

    static char getAlleleBase (int firstIntOfAllelePack) {
        return AlleleEncoder.alleleCodingToBaseMap.get(getAlleleCoding(firstIntOfAllelePack));
    }

    static byte getIndelLength (int firstIntOfAllelePack) {
        return (byte)(FastCall2.maxIndelLength & firstIntOfAllelePack);
    }

    static byte getAlleleCodingFromAlleleCodingLength(short alleleCodingLength) {
        int v = (alleleCodingLength>>>6);
        return (byte)v;
    }

    static char getAlleleBaseFromAlleleCodingLength(short alleleCodingLength) {
        return AlleleEncoder.alleleCodingToBaseMap.get(getAlleleCodingFromAlleleCodingLength(alleleCodingLength));
    }

    static byte getIndelLengthFromAlleleCodingLength(short alleleCodingLength) {
        return (byte)(FastCall2.maxIndelLength & alleleCodingLength);
    }
}
