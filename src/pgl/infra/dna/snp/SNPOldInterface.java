/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.dna.snp;

/**
 *
 * @author feilu
 */
public interface SNPOldInterface {
    /**
     * Remove duplicated alternative allele in alternative allele list and sort the sort the byte list
     */
    public void removeDuplicatedAltAlleles ();
    
    public void addAltAllele (char altAllele);
    
    public void addAltAlleleByte (byte alt);
    
    public void sortAltAlleles ();
    
    public int getAltAlleleIndex (byte alt);
    
    public byte getAltAlleleNumber ();
    
    public int getAltAlleleIndex (char altAllele);
    
    public byte getRefAlleleByte ();
    
    public char getRefAllele ();
    
    public byte getAltAlleleByte (int alleleIndex);
    
    public char getAltAllele (int alleleIndex);
    
}
