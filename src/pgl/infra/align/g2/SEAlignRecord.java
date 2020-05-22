/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.align.g2;

/**
 * Class holding information of single end read alignment from BWA-MEM
 * @author feilu
 */
public class SEAlignRecord implements Comparable<SEAlignRecord>{
    /**Query sequence*/
    String query = null;
    /**Chromosome*/
    String hit = "";
    /**Reference starting position of query sequence, inclusive*/
    int startPos = Integer.MIN_VALUE;
    /**Reference ending position of query sequence, exclusive*/
    int endPos = Integer.MIN_VALUE;
    /**Combination of bitwise FLAGs, 4 means unmapped by default*/
    short flag = 4;
    /**Mapping quality*/
    short mappingQuality = Short.MIN_VALUE;
    /**The length of matched range*/
    short alignMatchNumber = Short.MIN_VALUE;
    /**Mismatch number in the matched range*/
    short editDistance = Short.MIN_VALUE;
    
    /**
     * Construct an object
     */
    public SEAlignRecord () {

    }
    
    /**
     * Construct an object by all fields
     * @param query
     * @param hit
     * @param startPos
     * @param endPos
     * @param flag
     * @param mappingQuality
     * @param alignMatchNumber
     * @param editDistance 
     */
    public SEAlignRecord (String query, String hit, int startPos, int endPos, short flag, short mappingQuality, short alignMatchNumber, short editDistance) {
        this.query = query;
        this.hit = hit;
        this.startPos = startPos;
        this.endPos = endPos;
        this.flag = flag;
        this.mappingQuality = mappingQuality;
        this.alignMatchNumber = alignMatchNumber;
        this.editDistance = editDistance;
    }
    
    /**
     * Set query of an alignment record
     * @param query
     * @return 
     */
    public SEAlignRecord setQuery (String query) {
        this.query = query;
        return this;
    }
    
    /**
     * Set hit of an alignment record
     * @param hit
     * @return 
     */
    public SEAlignRecord setHit (String hit) {
        this.hit = hit;
        return this;
    }
    
    /**
     * Set start position of an alignment record
     * @param startPos
     * @return 
     */
    public SEAlignRecord setStartPos (int startPos) {
        this.startPos = startPos;
        return this;
    }
    
    /**
     * Build an object using set-builder method
     * @return 
     */
    public SEAlignRecord build () {
        return this;
    }
    
    /**
     * Return if the read is mapped
     * @return 
     */
    public boolean isMapped () {
        if (this.getHit().equals("")) return false;
        return true;
    }
    
    /**
     * Return if the alignment is supplementary
     * @return 
     */
    public boolean isSupplementaryAlignment () {
        return SAMUtils.isSupplementaryAlignment(this.flag);
    }
    
    /**
     * Return query of an alignment
     * @return
     */
    public String getQuery () {
        return this.query;
    }
    
    /**
     * Return hit of an alignment, e.g. chromosome name of the reference genome
     * @return NULL if not aligned
     */
    public String getHit () {
        return this.hit;
    }
    
    /**
     * Return the 1-based starting position of an alignment
     * @return Integer.MIN_VALUE if now aligned
     */
    public int getStartPos () {
        return this.startPos;
    }
    
    /**
     * Return the 1-based ending position of an alignment
     * @return Integer.MIN_VALUE if now aligned
     */
    public int getEndPos () {
        return this.endPos;
    }
    
    /**
     * Return the strand of an alignment, 1 representing plus, 0 representing minus
     * @return Byte.MIN_VALUE if now aligned
     */
    public byte getStrand () {
        if (!this.isMapped()) return Byte.MIN_VALUE;
        if (SAMUtils.isReverseAligned(this.flag)) return 0;
        return 1;
    }
    
    /**
     * Return the FLAG of the alignment
     * @return 
     */
    public short getFlag () {
        return this.flag;
    }
    
    /**
     * Return the mapping quality of an alignment
     * Note: 255 indicates the mapping quality is not available
     * @return Short.MIN_VALUE if now aligned
     */
    public short getMappingQuality () {
        return this.mappingQuality;
    }
    
    /**
     * Return the total number of alignment matched base ('M') in CIGAR string, including both sequence match and mismatch
     * @return Short.MIN_VALUE if now aligned
     */
    public short getAlignMatchNumber () {
        return this.alignMatchNumber;
    }
    
    /**
     * Return the length of the alignment
     * @return 
     */
    public int getAlignmentLength () {
        return this.endPos-this.startPos;
    }
    
    /**
     * Return number of mismatched base, relative to the reference genome
     * @return Short.MIN_VALUE if now aligned
     */
    public short getEditDistance () {
        return this.editDistance;
    }
    
    /**
     * Sort by alignment position
     * @param o
     * @return 
     */
    @Override
    public int compareTo(SEAlignRecord o) {
        int value = this.hit.compareTo(o.hit);
        if (value == 0) {
            return this.startPos-o.startPos;
        }
        else if (value < 0) return -1;
        return 1;        
    }
}
