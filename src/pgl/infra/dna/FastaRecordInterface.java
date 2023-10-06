/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.dna;

/**
 * Interface of fasta record.
 * @author feilu
 */
public interface FastaRecordInterface extends SequenceInterface {
    /**
     * Return the description of the record.
     * @return 
     */
    public String getDescription();
    
    /**
     * Return the ID of the record.
     * @return 
     */
    public int getID ();
    
    /**
     * Set the name of the record.
     * @param description
     */
    public void setDescription(String description);
    
    /**
     * Set the ID of the record.
     * @param id 
     */
    public void setID (int id);
    
    
}
