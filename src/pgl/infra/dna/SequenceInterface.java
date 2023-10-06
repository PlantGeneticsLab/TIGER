/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.dna;

/**
 * The interface holds basic methods for DNA sequence.
 * @author Fei Lu
 */
public interface SequenceInterface {

    /**
     * Return the length of sequence.
     * @return 
     */
    public int getSequenceLength ();
    
    
    /**
     * Return the proportion of base A.
     * @return 
     */
    public double getProportionA ();
    
    /**
     * Return the proportion of base T.
     * @return 
     */
    public double getProportionT ();
    
    /**
     * Return the proportion of base G.
     * @return 
     */
    public double getProportionG ();
    
    /**
     * Return the proportion of base C.
     * @return 
     */
    public double getProportionC ();
    
    /**
     * Return the proportion of base G and C.
     * @return 
     */
    public double getGCContent ();
    
    /**
     * Return the base of a given position
     * @param positionIndex is 0-based position
     * @return 
     */
    public char getBase (int positionIndex);
    
    /**
     * Return the entire DNA sequence.
     * @return 
     */
    public String getSequence();
    
    /**
     * Return a sub DNA sequence
     * @param startIndex inclusive
     * @param endIndex exclusive
     * @return 
     */
    public String getSequence(int startIndex, int endIndex);

    /**
     * Return a sub DNA sequence
     * @param startIndex inclusive
     * @param endIndex exclusive
     * @return
     */
    public SequenceInterface getSequenceInterface(int startIndex, int endIndex);

    /**
     * Return the reverse complementary sequence
     * @return 
     */
    public String getReverseComplementarySeq ();
    
    /**
     * Return a sub reverse complementary sequence
     * @param startIndex inclusive
     * @param endIndex exclusive
     * @return 
     */
    public String getReverseComplementarySeq (int startIndex, int endIndex);

    /**
     * Return if the sequence has "N" base.
     * @return 
     */
    public boolean isThereN ();
    
}
