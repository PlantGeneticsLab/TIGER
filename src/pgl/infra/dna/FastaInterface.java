/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.dna;

import pgl.infra.utils.IOFileFormat;

/**
 * Interface that provides operations on Fasta DNA sequences.
 * @author feilu
 */
public interface FastaInterface {
    
    /**
     * Read Fasta file from either txt or txt.gz file.
     * @param infileS
     * @param format 
     */
    public void readFasta (String infileS, IOFileFormat format);
        
    
    /**
     * Write Fasta file from selected sequences.
     * @param outfileS
     * @param ifOut 
     * @param format 
     */
    public void writeFasta (String outfileS, boolean[] ifOut, IOFileFormat format); 
    
    /**
     * Write Fasta file from a specified sequence
     * @param outfileS
     * @param index 
     * @param format 
     */
    public void writeFasta (String outfileS, int index, IOFileFormat format);
    
    /**
     * Write Fasta file.
     * @param outfileS 
     * @param format 
     */
    public void writeFasta (String outfileS, IOFileFormat format);
    
    /**
     * Return N50 statistic.
     * @return 
     */
    public int getN50 ();
    
    /**
     * Return L50 statistic.
     * @return 
     */
    public int getL50 ();
    
    /**
     * Return total sequence length in bp.
     * @return 
     */
    public long getTotalSeqLength ();
    
    /**
     * Return number of sequences.
     * @return 
     */
    public int getSeqNumber ();
    
    /**
     * Return index from sequence name,
     * the sequences will be sorted by description first if it is not.
     * @param name
     * @return 
     */
    public int getIndexByDescription(String name);
    
    /**
     * Return sequence length in bp.
     * @param index
     * @return 
     */
    public int getSeqLength (int index);
    
    /**
     * Return all of the sequence descriptions.
     * @return 
     */
    public String[] getDescriptions () ;
    
    /**
     * Return sequence description.
     * @param index
     * @return 
     */
    public String getDescription (int index) ;
    
    /**
     * Return sequence.
     * @param index
     * @return 
     */
    public String getSeq (int index) ;
    
    /**
     * Return a stretch of sequence.
     * @param index
     * @param startIndex inclusive
     * @param endIndex exclusive
     * @return 
     */
    public String getSeq (int index, int startIndex, int endIndex) ;
    
    /**
     * Return a FastaRecord {@link FastaRecordInterface}
     * @param index
     * @return 
     */ 
    public FastaRecordInterface getFastaRecord (int index);
    
    /**
     * Set sequence description.
     * @param description
     * @param index 
     */
    public void setDescription(String description, int index) ;
    
    /**
     * Sort the sequences of Fasta by description.
     */
    public void sortByDescription () ;
    
    /**
     * Sort the sequences of Fasta by description as value, e.g. when the names are numbers, 1, 2, 3...
     */
    public void sortByDescriptionValue ();
    
    /**
     * Sort the sequences of Fasta by ID.
     */
    public void sortByID () ;
    
    /**
     * Sort the sequences of Fasta by length in ascending order
     */
    public void sortByLengthAscending () ;
    
    /**
     * Sort the sequences of Fasta by length in descending order
     */
    public void sortByLengthDescending () ;
    
    /**
     * Return if the fasta has N in it
     * @return 
     */
    public boolean isThereN () ;

}
