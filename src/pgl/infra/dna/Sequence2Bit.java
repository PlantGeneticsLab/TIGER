/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.dna;

import java.util.Arrays;


/**
 * Class holding DNA sequence, with 2 bits for each base. Non-"ACGT" bases are not allowed.
 * This is designed to pack large sequence data set, e.g. wheat genome.
 * Well packed into memory, but lower speed performance than {@link pgl.infra.dna.SequenceByte}
 * @author feilu
 */
public class Sequence2Bit extends Sequence3Bit {
    private static final byte[] baseAscII = DNAUtils.getBaseAscIIArray();
    
    public Sequence2Bit(String seq) {
        super(seq, 2);
    }
    
    @Override
    protected void setCodedBase (int positionIndex, byte value) {
        int index = Arrays.binarySearch(baseAscII, value);
        if (index < 0) {
            System.out.println("Sequence contains non-ACGT bases, cannot be converted to 2 bits");
            System.exit(1);
        }
        byte v = ascIIByteMap.get(value);
        int startIndex = wordSize * positionIndex;
        for (int i = 0; i < wordSize; i++) {
            boolean b = (v & 1) == 1;
            v = (byte) (v >> 1);
            seqS.set(startIndex+i, b);
        }
    }
    
    @Override
    public boolean isThereN () {
        return false;
    }
}
