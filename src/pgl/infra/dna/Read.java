/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package pgl.infra.dna;

import java.util.Arrays;
import org.apache.commons.lang.ArrayUtils;


/**
 * Class holding Illumina read in Fastq format.
 * @author Fei Lu
 */
public class Read extends SequenceByte {
    /**sequence Identifier*/
    String ID;
    /**description*/
    String des;
    /**quality value string*/
    byte[] qualValue = null;
    byte[] sortedQualValue = null;
    
    public Read (String ID, String seq, String des, String qual, int phredScale) {
        super(seq);
        this.ID = ID;
        this.des = des;
        this.qualValue = qual.getBytes();
        this.calibrateQualValue(phredScale);
    }
    
    private void calibrateQualValue(int phredScale) {
        for (int i = 0; i < qualValue.length; i++) {
            qualValue[i] = (byte)(qualValue[i]-phredScale);
        }
    }
    
    /**
     * Return mean quality value of read
     * @return 
     */
    public byte returnMeanQuality () {
        double mean = 0;
        for (int i= 0; i < qualValue.length; i++) {
            mean+=qualValue[i];
        }
        return (byte)(mean/super.getSequenceLength());
    }
    
    /**
     * Return median quality of read
     * @return 
     */
    public byte getMedianQuality (int phredScale) {
        if (sortedQualValue == null) {
            System.arraycopy(qualValue, 0, sortedQualValue, 0, this.seqByte.length);
            Arrays.sort(sortedQualValue);
        }
        return sortedQualValue[sortedQualValue.length/2];
    }
    
    /**
     * Return quality value of base
     * @param index 0-based position index
     * @return 
     */
    public byte getBaseQuality (int index) {
        return qualValue[index];
    }
    
    /**
     * Return description of a read
     * @return 
     */
    public String getDescription () {
        return this.des;
    }
    
    /**
     * Return identifier of a read
     * @return 
     */
    public String getID () {
        return this.ID;
    }
    
    /**
     * Return reverse quality string
     * @param phredScale
     * @return 
     */
    public String getReverseQual (int phredScale) {
        byte[] qual = new byte[qualValue.length];
        System.arraycopy(qualValue, 0, qual, 0, qual.length);
        for (int i = 0; i < this.getSequenceLength(); i++) {
            qual[i]+=phredScale;
        }
        ArrayUtils.reverse(qual);
        return new String(qual);
    }
    
    /**
     * Return quality string
     * @param phredScale
     * @return 
     */
    public String getQualS (int phredScale) {
        byte[] qual = new byte[qualValue.length];
        System.arraycopy(qualValue, 0, qual, 0, qual.length);
        for (int i = 0; i < this.getSequenceLength(); i++) {
            qual[i]+=phredScale;
        }
        return new String(qual);
    }
    
    /**
     * Return quality string from startIndex to endIndex
     * @param startIndex
     * @param endIndex
     * @param phredScale
     * @return 
     */
    public String getQualS (int startIndex, int endIndex, int phredScale) {
        byte[] qual = new byte[endIndex - startIndex];
        System.arraycopy(qualValue, startIndex, qual, 0, endIndex - startIndex);
        for (int i = 0; i < this.getSequenceLength(); i++) {
            qual[i]+=phredScale;
        }
        return new String(qual);
    }
}
