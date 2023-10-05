/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.dna;

import com.koloboke.collect.map.hash.HashByteByteMap;
import com.koloboke.collect.map.hash.HashByteByteMaps;
import java.util.Arrays;
import java.util.BitSet;

/**
 * Class holding DNA sequence, with 3 bits for each base. Bases are converted to upper case. Non-"ACGT" bases are automatically converted to "N".
 * <P>
 * Well packed into memory, but lower speed performance than {@link pgl.infra.dna.SequenceByte}
 * This is designed to pack large sequence data set, e.g. wheat genome.
 * @author feilu
 */
public class Sequence3Bit implements SequenceInterface {
    protected static final HashByteByteMap ascIIByteMap = buildAscIIByteMap();
    protected static final HashByteByteMap byteAscIIMap = buildByteAscIIMap();
    BitSet seqS = null;
    int sequenceLength;
    protected int wordSize = 3;
    
    public Sequence3Bit () {
        
    }
    
    /**
     * Constructs the object from DNA sequence {@link java.lang.String}
     * @param seq 
     */
    public Sequence3Bit (String seq) {
        this.initialize(seq); 
    }

    /**
     * Construct the object from a BitSet, with 3bits for each base
     * @param seqS
     */
    public Sequence3Bit (BitSet seqS, int sequenceLength) {
        this.seqS = seqS;
        this.sequenceLength = sequenceLength;
    }

    /**
     * Constructs the object from DNA sequence {@link java.lang.String}
     * @param seq
     * @param wordSize 3 or 2
     */
    protected Sequence3Bit (String seq, int wordSize) {
        if (wordSize > 3 || wordSize < 2) {
            System.out.println("Word size should be either 3 or 2, program stops");
            System.exit(1);
        }
        this.wordSize = wordSize;
        this.initialize(seq);
    }
    
    /**
     * Initialize the object by setting bit set
     * @param seq 
     */
    private void initialize (String seq) {
        this.sequenceLength = seq.length();
        byte[] seqByte = seq.toUpperCase().getBytes();
        seqS = new BitSet(seq.length()*wordSize);
        for (int i = 0; i < seq.length(); i++) {
            this.setCodedBase(i, seqByte[i]);
        }
    }
    
    /**
     * Set the base in bit set from a position index of DNA sequence
     * @param positionIndex index of DNA sequence
     * @param value AscII value of a base
     */
    protected void setCodedBase (int positionIndex, byte value) {
        byte v = ascIIByteMap.get(value);
        int startIndex = wordSize * positionIndex;
        for (int i = 0; i < wordSize; i++) {
            boolean b = (v & 1) == 1;
            v = (byte) (v >> 1);
            seqS.set(startIndex+i, b);
        }
    }
    
    /**
     * Return a map converting ascII value of a base to required value.
     * A(00000000), C(00000001), G(00000010), T(0000000011), N(00000100)
     * @return 
     */
    private static HashByteByteMap buildAscIIByteMap () {
        if (ascIIByteMap != null) return ascIIByteMap;
        int size = 128;
        byte[] key = new byte[size];
        byte[] value = new byte[size];
        for (int i = 0; i < key.length; i++) {
            key[i] = (byte)i;
            value[i] = (byte)4;
        }
        value[65] = 0;
        value[67] = 1;
        value[71] = 2;
        value[84] = 3;
        value[78] = 4;
        HashByteByteMap ascIIByteMap = HashByteByteMaps.newImmutableMap(key, value);
        return ascIIByteMap;
    } 
    
    private static HashByteByteMap buildByteAscIIMap () {
        if (byteAscIIMap != null) return byteAscIIMap;
        int size = 5;
        byte[] key = new byte[size];
        byte[] value = new byte[size];
        for (int i = 0; i < size; i++) {
            key[i] = (byte)i;
        }
        value[0] = 65;
        value[1] = 67;
        value[2] = 71;
        value[3] = 84;
        value[4] = 78;
        HashByteByteMap byteAscIIMap = HashByteByteMaps.newImmutableMap(key, value);
        return byteAscIIMap;
    }
    
    @Override
    public int getSequenceLength() {
        return this.sequenceLength;
    }

    @Override
    public double getProportionA() {
        int cnt = 0;
        for (int i = 0; i < this.sequenceLength; i++) {
            if (this.getBase(i) == 'A') cnt++;
        }
        return (double)cnt/this.getSequenceLength();
    }

    @Override
    public double getProportionT() {
        int cnt = 0;
        for (int i = 0; i < this.sequenceLength; i++) {
            if (this.getBase(i) == 'T') cnt++;
        }
        return (double)cnt/this.getSequenceLength();
    }

    @Override
    public double getProportionG() {
        int cnt = 0;
        for (int i = 0; i < this.sequenceLength; i++) {
            if (this.getBase(i) == 'G') cnt++;
        }
        return (double)cnt/this.getSequenceLength();
    }

    @Override
    public double getProportionC() {
        int cnt = 0;
        for (int i = 0; i < this.sequenceLength; i++) {
            if (this.getBase(i) == 'C') cnt++;
        }
        return (double)cnt/this.getSequenceLength();
    }

    @Override
    public double getGCContent() {
        int cnt = 0;
        for (int i = 0; i < this.sequenceLength; i++) {
            if (this.getBase(i) == 'C' || this.getBase(i) == 'G') cnt++;
        }
        return (double)cnt/this.getSequenceLength();
    }

    @Override
    public char getBase (int positionIndex) {
        return (char)getBaseAscII(positionIndex);
    }
    
    protected byte getBaseAscII (int positionIndex) {
        return byteAscIIMap.get(getCodedBase(positionIndex));
    }
    
    /**
     * Return base coding in byte from a position
     * @param positionIndex
     * @return 
     */
    protected byte getCodedBase (int positionIndex) {
        byte value = 0;
        int startIndex = (positionIndex+1) * wordSize - 1;
        for (int i = 0; i < wordSize; i++) {
            if (seqS.get(startIndex-i)) {
                value = (byte)((value << 1) + 1);
            }
            else {
                value = (byte)((value << 1) + 0);
            }
        }
        return value;
    }
    
    /**
     * Print base coding in bits of the DNA sequence.
     * Used in tests
     */
    protected void printBits () {
        StringBuilder s = new StringBuilder();
        for( int i = 0; i < seqS.length();  i++ ) {
            s.append(seqS.get( i ) == true ? 1: 0 );
        }

        System.out.println( s );
    }
    
    @Override
    public String getSequence() {
        return getSequence(0, sequenceLength);
    }

    @Override
    public String getSequence(int startIndex, int endIndex) {
        return new String(getSequenceAscII(startIndex, endIndex));
    }
    
    @Override
    public SequenceInterface getSequenceInterface(int startIndex, int endIndex) {
        return this.getSequence3Bit(startIndex, endIndex);
    }

    /**
     * Return a {@link Sequence3Bit} based on position
     * @param startIndex inclusive
     * @param endIndex exclusive
     * @return
     */
    public Sequence3Bit getSequence3Bit (int startIndex, int endIndex) {
        BitSet bs = this.seqS.get((startIndex)*this.wordSize, (endIndex)*this.wordSize+1);
        return new Sequence3Bit (bs, endIndex-startIndex);
    }
    
    /**
     * Return byte array of ascII code from a stretch of DNA sequence
     * @param startIndex
     * @param endIndex
     * @return 
     */
    protected byte[] getSequenceAscII (int startIndex, int endIndex) {
        int size = endIndex - startIndex;
        byte[] values = new byte[size];
        for (int i = 0; i < size; i++) {
            values[i] = getBaseAscII(startIndex+i);
        }
        return values;
    }

    @Override
    public String getReverseComplementarySeq() {
        return getReverseComplementarySeq (0, sequenceLength);
    }

    @Override
    public String getReverseComplementarySeq(int startIndex, int endIndex) {
        HashByteByteMap baseCompleByteMap = DNAUtils.getBaseCompleAscIIMap();
        byte[] reverseByte = new byte[endIndex - startIndex];
        for (int i = 0; i < reverseByte.length; i++) {
            reverseByte[i] = baseCompleByteMap.get(getBaseAscII(endIndex-i-1));
        }
        return new String(reverseByte);
    }

    @Override
    public boolean isThereNonACGTNBase() {
        byte[] baseByteWithN = DNAUtils.getBaseWithNAscIIArray();
        for (int i = 0; i < this.getSequenceLength(); i++) {
            int index = Arrays.binarySearch(baseByteWithN, this.getBaseAscII(i));
            if (index < 0) {
                return true;
            }
        }
        return false;
    }
    
    @Override
    public boolean isThereN () {
        for (int i = 0; i < this.getSequenceLength(); i++) {
            if (this.getBaseAscII(i) == 78) return true;
        }
        return false;
    }
}
