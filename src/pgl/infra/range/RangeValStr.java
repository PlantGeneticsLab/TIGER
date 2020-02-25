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
public class RangeValStr extends RangeVal {
    /**
     * The strand of this range, 1 means plus, 0 means minus
     */
    public byte str = -1;
    
    /**
     * Construct a {@link RangeValStr} object
     * @param chr
     * @param start
     * @param end
     * @param val
     * @param str 
     */
    public RangeValStr(int chr, int start, int end, double val, byte str) {
        super(chr, start, end, val);
        this.str = str;
    }
    
    /**
     * Construct a {@link RangeValStr} object
     * @param r
     * @param val
     * @param str 
     */
    public RangeValStr (Range r, double val, byte str) {
        this(r.getRangeChromosome(), r.getRangeStart(), r.getRangeEnd(), val, str);
    }
    
    @Override
    public String getInfoString () {
       StringBuilder sb = new StringBuilder();
       sb.append(chr).append("\t").append(start).append("\t").append(end).append("\t").append(value).append("\t").append(str);
       return sb.toString();
    }
    
    /**
     * Return the strand of the range, 1 means plus, 0 means minus
     * @return 
     */
    public byte getRangeStrand () {
        return str;
    }
    
    /**
     * Set the value of the range
     * @param str
     */
    public void setRangeStrand (byte str) {
        this.str = str;
    }
    
    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final RangeValStr other = (RangeValStr) obj;
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
        if (this.str != other.str) {
            return false;
        }
        return true;
    }
}
