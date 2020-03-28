/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.dna;

/**
 * Interface of FastaRecord
 * @author feilu
 */
public interface FastaRecordInterface extends SequenceInterface {
    /**
     * Return the name of the record
     * @return 
     */
    public String getName ();
    
    /**
     * Return the ID of the record
     * @return 
     */
    public int getID ();
    
    /**
     * Set the name of the record
     * @param newName 
     */
    public void setName (String newName);
    
    /**
     * Set the ID of the record
     * @param id 
     */
    public void setID (int id);
    
    
}
