/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.dna;

import com.koloboke.collect.map.hash.HashByteByteMap;
import java.util.Arrays;

/**
 * The class uses one byte to store a DNA base. 
 * <p>
 * Bases are converted to upper case, with full representation of Non-"ATGC" bases.
 * It supports standard IUPAC DNA coding (https://www.bioinformatics.org/sms/iupac.html).
 * Higher memory cost, but higher speed than {@link pgl.infra.dna.Sequence3Bit}
 * @author Fei Lu
 */
public class SequenceByte implements SequenceInterface {
    
    byte[] seqByte = null;
    
    /**
     * Constructs a {@code SequenceByte} from {@code String}. The lower case bases are converted to upper case.
     * @param seq 
     *        The name of a supported DNA sequence
     */
    public SequenceByte (String seq) {
        seqByte = seq.toUpperCase().getBytes();
    }
    
    @Override
    public int getSequenceLength() {
        return seqByte.length;
    }

    @Override
    public double getProportionA() {
        int cnt = 0;
        for (int i = 0; i < this.getSequenceLength(); i++) {
            if (this.seqByte[i] == 65) cnt++;
        }
        return (double)cnt/this.seqByte.length;
    }

    @Override
    public double getProportionT() {
        int cnt = 0;
        for (int i = 0; i < this.getSequenceLength(); i++) {
            if (this.seqByte[i] == 84) cnt++;
        }
        return (double)cnt/this.seqByte.length;
    }

    @Override
    public double getProportionG() {
        int cnt = 0;
        for (int i = 0; i < this.getSequenceLength(); i++) {
            if (this.seqByte[i] == 71) cnt++;
        }
        return (double)cnt/this.seqByte.length;
    }

    @Override
    public double getProportionC() {
        int cnt = 0;
        for (int i = 0; i < this.getSequenceLength(); i++) {
            if (this.seqByte[i] == 67) cnt++;
        }
        return (double)cnt/this.seqByte.length;
    }

    @Override
    public double getGCContent() {
        int cnt = 0;
        for (int i = 0; i < this.getSequenceLength(); i++) {
            if (this.seqByte[i] == 67 || this.seqByte[i] == 84) cnt++;
        }
        return (double)cnt/this.seqByte.length;
    }

    @Override
    public char getBase(int positionIndex) {
        return (char)this.seqByte[positionIndex];
    }
    
    /**
     * Return the ascII value of a base
     * @param positionIndex
     * @return 
     */
    public byte getBaseByte (int positionIndex) {
        return this.seqByte[positionIndex];
    }
    
    @Override
    public String getSequence() {
        return new String(this.seqByte);
    }

    @Override
    public String getSequence(int startIndex, int endIndex) {
        return new String(this.seqByte, startIndex, endIndex-startIndex);
    }

    @Override
    public String getReverseComplementarySeq() {
        return this.getReverseComplementarySeq(0, this.getSequenceLength());
    }

    @Override
    public String getReverseComplementarySeq(int startIndex, int endIndex) {
        HashByteByteMap baseCompleByteMap = DNAUtils.getBaseCompleAscIIMap();
        byte[] reverseByte = new byte[endIndex - startIndex];
        for (int i = 0; i < reverseByte.length; i++) {
            reverseByte[i] = baseCompleByteMap.get(seqByte[endIndex-i-1]);
        }
        return new String(reverseByte);
    }
    
    @Override
    public boolean isThereN () {
        for (int i = 0; i < this.getSequenceLength(); i++) {
            if (this.seqByte[i] == 78) return true;
        }
        return false;
    }
    
    @Override
    public boolean isThereNonACGTNBase () {
        byte[] baseByteWithN = DNAUtils.getBaseWithNAscIIArray();
        for (int i = 0; i < this.getSequenceLength(); i++) {
            int index = Arrays.binarySearch(baseByteWithN, this.getBaseByte(i));
            if (index < 0) {
                System.out.println(this.getBase(i));
                return true;
            }
        }
        return false;
    }
}
