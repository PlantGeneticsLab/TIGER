/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.dna.snp;

import pgl.infra.dna.allele.AlleleEncoder;
import pgl.infra.pos.ChrPos;
import gnu.trove.list.array.TByteArrayList;
import gnu.trove.set.hash.TByteHashSet;

/**
 *
 * @author feilu
 */
public class SNPOld extends ChrPos implements SNPOldInterface {
    byte ref = Byte.MIN_VALUE;
    TByteArrayList alts = new TByteArrayList();
    
    public SNPOld (short chr, int pos, char refAllele, char altAllele) {
        super(chr, pos);
        this.ref = AlleleEncoder.alleleBaseToCodingMap.get(refAllele);
        this.addAltAllele(altAllele);
    }
    
    public SNPOld (short chr, int pos, byte ref, byte alt) {
        super(chr, pos);
        this.ref = ref;
        this.addAltAlleleByte(alt);
    }
    
    public SNPOld (short chr, int pos, byte ref, TByteArrayList alts) {
        super(chr, pos);
        this.ref = ref;
        this.alts = alts;
    }
    
    public TByteArrayList getAltAlleleList () {
        return this.alts;
    }
    
    @Override
    public void removeDuplicatedAltAlleles () {
        TByteHashSet s = new TByteHashSet(alts);
        alts = new TByteArrayList(s);
        alts.sort();
    }
    
    @Override
    public void addAltAllele (char altAllele) {
        this.addAltAlleleByte(AlleleEncoder.alleleBaseToCodingMap.get(altAllele));
    }
    
    @Override
    public void addAltAlleleByte (byte alt) {
        this.alts.add(alt);
        if (alts.size() == Byte.MAX_VALUE) this.removeDuplicatedAltAlleles();
    }
    
    @Override
    public void sortAltAlleles () {
        alts.sort();
    }
    
    @Override
    public int getAltAlleleIndex (byte alt) {
        return alts.binarySearch(alt);
    }
    
    @Override
    public byte getAltAlleleNumber () {
        return (byte)this.alts.size();
    }
    
    @Override
    public int getAltAlleleIndex (char altAllele) {
        return this.getAltAlleleIndex(AlleleEncoder.alleleBaseToCodingMap.get(altAllele));
    }
    
    @Override
    public byte getRefAlleleByte() {
        return this.ref;
    }

    @Override
    public char getRefAllele() {
        return AlleleEncoder.alleleCodingToBaseMap.get(this.getRefAlleleByte());
    }

    @Override
    public byte getAltAlleleByte(int altIndex) {
        return this.alts.get(altIndex);
    }

    @Override
    public char getAltAllele(int alleleIndex) {
        return AlleleEncoder.alleleCodingToBaseMap.get(this.getAltAlleleByte(alleleIndex));
    }
}
