/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.range;

/**
 * Defining major field in this class, which are chr, start, and end. Note the maximum number of chromosomes are 65536.
 * @author feilu
 */
public abstract class RangeAbstract implements RangeInterface {
    public short chr;
    public int start;
    /**Exclusive*/
    public int end;
    
    public RangeAbstract (int chr, int start, int end) {
        this.chr = (short)chr;
        this.start = start;
        this.end = end;
    }
    
    @Override
    public int getRangeChromosome () {
        return this.chr;
    }
    
    @Override
    public int getRangeStart () {
        return this.start;
    }
    
    @Override
    public int getRangeEnd () {
        return this.end;
    }
    
    @Override
    public int getRangeSize () {
        return this.getRangeEnd() - this.getRangeStart();
    }
    
    @Override
    public void setRangeChromosome (int chr) {
        this.chr = (short)chr;
    }
    
    @Override
    public void setRangeStart (int start) {
        this.start = start;
    }
    
    @Override
    public void setRangeEnd (int end) {
        this.end = end;
    }
    
    
    @Override
    public String getInfoString () {
       StringBuilder sb = new StringBuilder();
       sb.append(chr).append("\t").append(start).append("\t").append(end);
       return sb.toString();
    }
    
    @Override
    public boolean isOverlap (RangeInterface ri) {
        if (this.getRangeChromosome() != ri.getRangeChromosome()) return false;
        if (this.getRangeEnd() <= ri.getRangeStart() || ri.getRangeEnd() <= this.getRangeStart()) return false;
        return true;
    }
    
    @Override
    public boolean isWithin (RangeInterface ri) {
        if (this.getRangeChromosome() != ri.getRangeChromosome()) return false;
        if (this.getRangeStart() >= ri.getRangeStart() && this.getRangeEnd() <= ri.getRangeEnd()) return true;
        return false;
    }
    
    @Override
    public boolean isContain (RangeInterface ri) {
        if (this.getRangeChromosome() != ri.getRangeChromosome()) return false;
        if (this.getRangeStart() <= ri.getRangeStart() && this.getRangeEnd() >= ri.getRangeEnd()) return true;
        return false;
    }
    
    @Override
    public boolean isContain (int chr, int pos) {
        if (this.getRangeChromosome() != chr) return false;
        if (this.getRangeStart() <= pos && this.getRangeEnd() > pos) return true;
        return false;
    }
    
    public RangeInterface getIntersection (RangeAbstract other) {
        if (!this.isOverlap(other)) return null;
        if (this.start < other.start) return new Range(this.chr, other.start, this.end < other.end? this.end : other.end);
        return new Range(this.chr, this.start, this.end < other.end? this.end : other.end);
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final RangeAbstract other = (RangeAbstract) obj;
        if (this.chr != other.chr) {
            return false;
        }
        if (this.start != other.start) {
            return false;
        }
        if (this.end != other.end) {
            return false;
        }
        return true;
    }
    
    @Override
    public int hashCode() {
        return start;
    }
    
    @Override
    public int compareTo(RangeInterface ri) {
        if (chr == ri.getRangeChromosome()) {
            if (start == ri.getRangeStart()) {
                return 0;
//                if (end == ri.getRangeEnd()) return 0;
//                else if (end < ri.getRangeEnd()) return -1;
//                return 1;
            }
            else if (start < ri.getRangeStart()) return -1;
            return 1;
        }
        else if (chr < ri.getRangeChromosome()) return -1;
        return 1;
    }
}
