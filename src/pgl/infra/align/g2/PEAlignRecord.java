/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.align.g2;


/**
 * Class holding information of paired end read alignment from BWA-MEM
 * @author feilu
 */
public class PEAlignRecord implements Comparable<PEAlignRecord> {
    SEAlignRecord r1 = null;
    SEAlignRecord r2 = null;

    /**
     * Construct an object
     */
    public PEAlignRecord () {
        
    }

    /**
     * Construct an object from two {@link SEAlignRecord} objects
     * @param r1
     * @param r2
     */
    public PEAlignRecord (SEAlignRecord r1, SEAlignRecord r2) {
        this.r1 = r1;
        this.r2 = r2;
    }

    /**
     * Set r1 {@link SEAlignRecord}
     * @param r1
     * @return
     */
    public PEAlignRecord setR1AlignmentRecord (SEAlignRecord r1) {
        this.r1 = r1;
        return this;
    }

    /**
     * Set r2 {@link SEAlignRecord}
     * @param r2
     * @return
     */
    public PEAlignRecord setR2AlignmentRecord (SEAlignRecord r2) {
        this.r2 = r2;
        return this;
    }

    /**
     * Build an object using set-builder method
     * @return
     */
    public PEAlignRecord build () {
        return this;
    }

    /**
     * Return insert size of a PE read
     * @return Short.MIN_VALUE if the read is not properly aligned
     */
    public short getInsertSize () {
        if (!SAMUtils.isReadsMappedInProperPair(r1.getFlag())) return Short.MIN_VALUE;
        int length = 0;
        if (r1.getStartPos() < r2.getStartPos()) return (short)(r2.getEndPos()-r1.getStartPos());
        return (short)(r1.getEndPos()-r2.getStartPos());
    }

    /**
     * Return the r1 alignment record
     * @return
     */
    public SEAlignRecord getR1AlignmentRecord () {
        return r1;
    }

    /**
     * Return the r2 alignment record
     * @return
     */
    public SEAlignRecord getR2AlignmentRecord () {
        return r2;
    }

    @Override
    public int compareTo(PEAlignRecord o) {
        return this.getR1AlignmentRecord().compareTo(o.getR1AlignmentRecord());
    }
}
