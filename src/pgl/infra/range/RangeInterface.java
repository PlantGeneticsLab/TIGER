/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.range;

/**
 * Interface holding basic for a range.
 * The range is specially designed for the interval on chromosomes. Hence, it has an additional attribute of chromosome.
 * @author feilu
 */
public interface RangeInterface extends Comparable<RangeInterface> {
    
    /**
     * Return the chromosome of a range, 1-based
     * @return 
     */
    public int getRangeChromosome ();
    
    /**
     * Return the start position of a range, 1-based, inclusive
     * @return 
     */
    public int getRangeStart ();
    
    /**
     * Return the ending position of a range, 1-based, exclusive
     * @return 
     */
    public int getRangeEnd ();
    
    /**
     * Return the size of a range
     * @return 
     */
    public int getRangeSize ();
    
    /**
     * Set the chromosome of range
     * @param chr 
     */
    public void setRangeChromosome (int chr);
    
    /**
     * Set the start position of a range
     * @param start 
     */
    public void setRangeStart (int start);
    
    /**
     * Set the end position of a range
     * @param end 
     */
    public void setRangeEnd (int end);
    
    /**
     * Return the information string, which is "Chr\tStart\tEnd"
     * @return 
     */
    public String getInfoString ();
    
    /**
     * Test if two ranges have overlap
     * @param ri
     * @return 
     */
    public boolean isOverlap (RangeInterface ri);
    
    /**
     * Test if the range is within another one
     * @param ri
     * @return 
     */
    public boolean isWithin (RangeInterface ri);
    
    /**
     * Test if the range contains another one
     * @param ri
     * @return 
     */
    public boolean isContain (RangeInterface ri);
    
    /**
     * Test if the range contains a position
     * @param chr
     * @param pos
     * @return 
     */
    public boolean isContain (int chr, int pos);
    
    /**
     * Return the intersection of two ranges.
     * @param other
     * @return null if there is no overlap
     */
    public RangeInterface getIntersection (RangeInterface other);
    
    /**
     * Compare ranges by starting position
     * @param ri
     * @return 
     */
    @Override
    public int compareTo(RangeInterface ri);
    
}
