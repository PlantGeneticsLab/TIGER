/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.app.grt;

import pgl.infra.pos.ChrPos;

/**
 *
 * @author feilu
 */
public class AlleleInfo extends ChrPos {
    byte allele = Byte.MIN_VALUE;
    byte base = Byte.MIN_VALUE;
    byte end = Byte.MIN_VALUE; // paired-end, 1 or 2
    byte relaPos = Byte.MIN_VALUE;
    
    public AlleleInfo(short chr, int pos) {
        super(chr, pos);
    }
    
    public AlleleInfo (short chr, int pos, byte end) {
        super(chr, pos);
        this.end = end;
    }
    
    public AlleleInfo (short chr, int pos, byte allele, byte base, byte end, byte relaPos) {
        super(chr, pos);
        this.allele = allele;
        this.base = base;
        this.end = end;
        this.relaPos = relaPos;
    }
    
    public byte getAllele () {
        return allele;
    }
    
    public byte getBase () {
        return base;
    }
    
    public byte getEnd () {
        return end;
    }
    
    public byte getRelativePosition () {
        return relaPos;
    }
    
    public void setRelativePosition (byte relaPos) {
        this.relaPos = relaPos;
    }
    
    public void setAllele (byte allele) {
        this.allele = allele;
    }
    
    public void setBase (byte base) {
        this.base = base;
    }
}
