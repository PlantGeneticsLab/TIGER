/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.dna;

import com.koloboke.collect.map.hash.HashByteByteMap;
import java.util.Arrays;

/**
 * The class uses AscII value to store a DNA base.
 * <p>
 * Bases are converted to upper case, with full representation of Non-"ATGC" bases.
 * It supports standard IUPAC DNA coding (https://www.bioinformatics.org/sms/iupac.html).
 * Higher memory cost, but higher speed than {@link pgl.infra.dna.Sequence3Bit}
 * @author Fei Lu
 */
public class SequenceByte implements SequenceInterface {
    
    byte[] seqAscII = null;

    /**
     * Construct an object.
     */
    public SequenceByte () {
        
    }
    
    /**
     * Constructs an object from DNA sequence in {@link String}, the lower case bases are converted to upper case.
     * @param seq A string of DNA sequence
     */
    public SequenceByte (String seq) {
        seqAscII = seq.toUpperCase().getBytes();
    }

    /**
     * Constructs an object from an array of AscII value of DNA sequence
     * @param seqAscII
     */
    public SequenceByte (byte[] seqAscII) {
        this.seqAscII = seqAscII;
    }
    
    @Override
    public int getSequenceLength() {
        return seqAscII.length;
    }

    @Override
    public double getProportionA() {
        int cnt = 0;
        for (int i = 0; i < this.getSequenceLength(); i++) {
            if (this.seqAscII[i] == 65) cnt++;
        }
        return (double)cnt/this.seqAscII.length;
    }

    @Override
    public double getProportionT() {
        int cnt = 0;
        for (int i = 0; i < this.getSequenceLength(); i++) {
            if (this.seqAscII[i] == 84) cnt++;
        }
        return (double)cnt/this.seqAscII.length;
    }

    @Override
    public double getProportionG() {
        int cnt = 0;
        for (int i = 0; i < this.getSequenceLength(); i++) {
            if (this.seqAscII[i] == 71) cnt++;
        }
        return (double)cnt/this.seqAscII.length;
    }

    @Override
    public double getProportionC() {
        int cnt = 0;
        for (int i = 0; i < this.getSequenceLength(); i++) {
            if (this.seqAscII[i] == 67) cnt++;
        }
        return (double)cnt/this.seqAscII.length;
    }

    @Override
    public double getGCContent() {
        int cnt = 0;
        for (int i = 0; i < this.getSequenceLength(); i++) {
            if (this.seqAscII[i] == 67 || this.seqAscII[i] == 84) cnt++;
        }
        return (double)cnt/this.seqAscII.length;
    }

    @Override
    public char getBase(int positionIndex) {
        return (char)this.seqAscII[positionIndex];
    }
    
    /**
     * Return the ascII value of a base
     * @param positionIndex
     * @return 
     */
    public byte getBaseAscII(int positionIndex) {
        return this.seqAscII[positionIndex];
    }
    
    @Override
    public String getSequence() {
        return new String(this.seqAscII);
    }

    @Override
    public String getSequence(int startIndex, int endIndex) {
        return new String(this.seqAscII, startIndex, endIndex-startIndex);
    }
    
    @Override
    public SequenceInterface getSequenceInterface(int startIndex, int endIndex) {
        return this.getSequenceByte(startIndex, endIndex);
    }

    /**
     * Return a {@link SequenceByte} based on positions
     * @param startIndex inclusive
     * @param endIndex exclusive
     * @return
     */
    public SequenceByte getSequenceByte (int startIndex, int endIndex) {
        byte[] bs = new byte[endIndex-startIndex];
        System.arraycopy(this.seqAscII, startIndex, bs, 0, endIndex-startIndex);
        return new SequenceByte(bs);
    }
    
    @Override
    public String getReverseComplementarySeq() {
        return this.getReverseComplementarySeq(0, this.getSequenceLength());
    }

    @Override
    public String getReverseComplementarySeq(int startIndex, int endIndex) {
        HashByteByteMap baseCompleByteMap = SequenceUtils.getBaseCompleAscIIMap();
        byte[] reverseByte = new byte[endIndex - startIndex];
        for (int i = 0; i < reverseByte.length; i++) {
            reverseByte[i] = baseCompleByteMap.get(seqAscII[endIndex-i-1]);
        }
        return new String(reverseByte);
    }
    
    @Override
    public boolean isThereN () {
        for (int i = 0; i < this.getSequenceLength(); i++) {
            if (this.seqAscII[i] == 78) return true;
        }
        return false;
    }

    /**
     * Return if the sequence has non-“ACGTN” base, e.g ".".
     * @return
     */
    public boolean isThereNonACGTNBase () {
        byte[] baseByteWithN = SequenceUtils.getBaseWithNAscIIArray();
        for (int i = 0; i < this.getSequenceLength(); i++) {
            int index = Arrays.binarySearch(baseByteWithN, this.getBaseAscII(i));
            if (index < 0) {
                System.out.println(this.getBase(i));
                return true;
            }
        }
        return false;
    }
}
