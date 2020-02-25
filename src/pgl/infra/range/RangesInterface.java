/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.range;

import java.util.List;

/**
 * Holding basic method of a range
 * @author feilu
 */
public interface RangesInterface {
    
    /**
     * Return a {@link Range} object
     * @param <T>
     * @param rangeIndex
     * @return 
     */
    public <T extends Range> T getRange (int rangeIndex);
    
    /**
     * Return the whole range list
     * @return 
     */
    public List<? extends Range> getRangeList();
    /**
     * Sort by range size
     */
    public void sortBySize ();
    
    /**
     * Sort by starting position of a range
     */
    public void sortByStartPosition ();
     
    /**
     * Insert a {@link Range} into the list
     * @param <T> {@link Range}
     * @param rangeIndex
     * @param r 
     * @return true if it is successful, otherwise false
     */
    public <T extends Range> boolean insertRange  (int rangeIndex, T r);
    
    /**
     * Remove a {@link Range} from the list
     * @param rangeIndex 
     */
    public void removeRange (int rangeIndex);
    
    /**
     * Set a {@link Range} in the list
     * @param <T> {@link Range}
     * @param rangeIndex
     * @param r 
     * @return true if it is successful, otherwise false
     */
    public <T extends Range> boolean setRange (int rangeIndex, T r);
    
    /**
     * Return total number of {@link Range} in the list
     * @return 
     */
    public int getRangeNumber ();
    
    /**
     * Return the starting index of range on a chromosome. 
     * @param chr
     * @return negative value if the chromosome does not have a range
     */
    public int getStartIndexOfChromosome (int chr);
    
    /**
     * Return the ending index of range on a chromosome. 
     * @param chr
     * @return negative value if the chromosome does not have a range
     */
    public int getEndIndexOfChromosome (int chr);
    
    /**
     * Return all the chromosomes in an array
     * @return 
     */
    public int[] getChromosomes ();
    
    /**
     * Return the total number of chromosomes
     * @return 
     */
    public int getChromosomeNumber ();
    
    /**
     * Return the starting position of a range
     * @param rangeIndex
     * @return 
     */
    public int getRangeStart (int rangeIndex);
    
    /**
     * Return the ending position of a range
     * @param rangeIndex
     * @return 
     */
    public int getRangeEnd (int rangeIndex);
    
    /**
     * Return the chromosome of a range
     * @param rangeIndex
     * @return 
     */
    public int getRangeChromosome (int rangeIndex);
    
    /**
     * Return a {@link Ranges} from a chromosome
     * @param <T>
     * @param chr
     * @return null if the chromosome does not have a range
     */
    public <T extends RangesInterface> T  getRangesByChromosome (int chr);
    
    /**
     * Return a {@link Ranges}, in which a position is contained
     * @param <T>
     * @param chr
     * @param pos
     * @return 
     */
    public <T extends RangesInterface> T  getRangesContainsPosition (int chr, int pos);
    
    /**
     * Return the indices of Range, where a position is contained
     * @param chr
     * @param pos
     * @return 
     */
    public int[] getRangesIndicesContainsPosition (int chr, int pos);
    
    /**
     * Return the first index of range, in which a position is contained
     * @param chr
     * @param pos
     * @return negative value if the {@link RangeInterface} does not contain the position
     */
    public int getFirstRangeIndex (int chr, int pos);
    
    /**
     * Return a {@link Ranges} with no overlap between individual range
     * @return 
     */
    public Ranges getNonOverlapRanges ();
    
    /**
     * Return a {@link Ranges} by adding another Ranges object
     * @param <T>
     * @param rs
     * @return 
     */
    public <T extends RangesInterface> T getMergedRanges(T rs);
    
    /**
     * Add another {@link Range} to the current list
     * @param <T> {@link Range}
     * @param r 
     * @return  
     */
    public <T extends Range> boolean addRange(T r);

    /**
     * Add {@link Ranges} to the current list
     * @param <T> {@link RangesInterface}
     * @param rs 
     * @return  true if it is successful, otherwise false
     */
    public <T extends RangesInterface> boolean addRanges(T rs);
    
}
