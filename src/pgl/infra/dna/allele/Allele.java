package pgl.infra.dna.allele;

/**
 * Class holding allele and its basic information
 * Basic information are stored in bits
 * @author feilu
 */
public class Allele {
    byte baseVal = -1;
    byte feature = 0;

    /**
     * Construct an empty object of {@link Allele}
     */
    public Allele () {

    }

    /**
     * Construct an object of {@link Allele}
     * @param c
     */
    public Allele (char c) {
        this.baseVal = AlleleEncoder.getAlleleByteFromBase(c);
    }

    /**
     * Construct an object of {@link Allele}
     * @param alleleByte
     */
    public Allele (byte alleleByte) {
        this.baseVal = alleleByte;
    }

    /**
     * Return the base of the allele
     * @return
     */
    public char getAlleleBase () {
        return AlleleEncoder.getAlleleBaseFromByte(this.baseVal);
    }

    /**
     * Return the byte code of the allele, see {@link AlleleEncoder}
     * @return
     */
    public byte getAlleleByte () {
        return this.baseVal;
    }

    /**
     * Return the code of allele feature
     * @return
     */
    public byte getAlleleFeature () {
        return this.feature;
    }

    /**
     * Set allele type, see {@link AlleleType}
     * @param at
     */
    public void setAlleleType (AlleleType at) {
        this.feature = (byte)(feature | at.getFeature());
    }

    /**
     * Set the allele feature at once, see {@link AlleleType}
     * @param feature
     */
    public void setAlleleFeature (byte feature) {
        this.feature = feature;
    }

    /**
     * Remove allele type, see {@link AlleleType}
     * @param at
     */
    public void removeAlleleType (AlleleType at) {
        this.feature = (byte) (feature & (~at.getFeature()));
    }

    /**
     * Reset all allele types to false, see {@link AlleleType}
     */
    public void resetAlleleTypeToDefault () {
        this.feature = 0;
    }

    /**
     * Return the allele is a specific type, see {@link AlleleType}
     * @param at
     * @return
     */
    public boolean isAlleleTypeOf (AlleleType at) {
        if ((this.feature & at.getFeature()) == 0) return false;
        return true;
    }
}
