/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.pos;

/**
 * Class of a genomic position, with position comparator can be extended
 * @author feilu
 */
public class ChrPos implements ChrPosInterface {
    short chr;
    int pos;

    /**
     * Construct an object
     */
    public ChrPos () {
        
    }

    /**
     * Construct an object with chromosome and position
     * @param chr
     * @param pos
     */
    public ChrPos (short chr, int pos) {
        this.chr = chr;
        this.pos = pos;
    }

    /**
     * Return the current instance
     * @return
     */
    public ChrPos getChrPos () {
        return this;
    }
    
    @Override
    public short getChromosome () {
        return chr;
    }
    
    @Override
    public int getPosition() {
        return pos;
    }

    @Override
    public void setChromosome (short chr) {
        this.chr = chr;
    }

    @Override
    public void setPosition (int position) {
        this.pos = position;
    }

    /**
     * The comparator based on the site position on chromosome
     * @param o
     * @return
     */
    @Override
    public int compareTo(ChrPosInterface o) {
        if (this.getChromosome() == o.getChromosome()) {
            if (this.getPosition() == o.getPosition()) return 0;
            else if (this.getPosition() < o.getPosition()) return -1;
            return 1;
        }
        else if (this.getChromosome() < o.getChromosome()) return -1;
        return 1;
    }
}
