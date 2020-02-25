/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.range;

/**
 *
 * @author feilu
 */
public class RangeVal extends Range {
    public double value;
    
    /**
     * Construct a {@link RangeVal} object
     * @param chr
     * @param start
     * @param end
     * @param val 
     */
    public RangeVal(int chr, int start, int end, double val) {
        super(chr, start, end);
        this.value = val;
    }
    
    /**
     * Construct a {@link RangeVal} object
     * @param r
     * @param val 
     */
    public RangeVal (Range r, double val) {
        this(r.getRangeChromosome(), r.getRangeStart(), r.getRangeEnd(), val);
    }
    
    @Override
    public String getInfoString () {
       StringBuilder sb = new StringBuilder();
       sb.append(chr).append("\t").append(start).append("\t").append(end).append("\t").append(value);
       return sb.toString();
    }
    
    /**
     * Return the value of the range
     * @return 
     */
    public double getRangeValue () {
        return value;
    }
    
    /**
     * Set the value of the range
     * @param value 
     */
    public void setRangeValue (double value) {
        this.value = value;
    }
    
    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final RangeVal other = (RangeVal) obj;
        if (this.chr != other.chr) {
            return false;
        }
        if (this.start != other.start) {
            return false;
        }
        if (this.end != other.end) {
            return false;
        }
        if (this.value != other.value) {
            return false;
        }
        return true;
    }
}
