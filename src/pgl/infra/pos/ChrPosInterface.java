/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.pos;

/**
 * Interface for methods of a genomic site
 * @author feilu
 */
public interface ChrPosInterface extends Comparable<ChrPosInterface> {

    /**
     * Return the chromosome of a genomic site
     * @return
     */
    public short getChromosome ();

    /**
     * Return the position of genomic site
     * @return
     */
    public int getPosition ();

    /**
     * Set the chromosome of genomic site
     * @param chr
     */
    public void setChromosome (short chr);

    /**
     * Set the position of a genomic site
     * @param position
     */
    public void setPosition (int position);
    
}
