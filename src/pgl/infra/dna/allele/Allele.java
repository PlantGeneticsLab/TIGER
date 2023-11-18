package pgl.infra.dna.allele;

/**
 * Class holding allele and its basic information.
 * <p>
 * Basic information are stored in 8 bits as a feature.
 * @author feilu
 */
public class Allele {
    byte alleleCoding = -1;
    byte feature = 0;

    /**
     * Construct an object of {@link Allele}.
     * @param c
     */
    public Allele (char c) {
        this.alleleCoding = AlleleEncoder.getAlleleCodingFromBase(c);
        this.checkBaseCoding();
    }

    /**
     * Construct an object of {@link Allele}.
     * @param alleleCoding
     */
    public Allele (byte alleleCoding) {
        this.alleleCoding = alleleCoding;
        this.checkBaseCoding();
    }

    /**
     * Double check allele coding to avoid incorrect input.
     */
    private void checkBaseCoding() {
        if (this.alleleCoding < AlleleEncoder.alleleCodings[0] || this.alleleCoding > AlleleEncoder.alleleCodings[AlleleEncoder.alleleCodings.length-1]) {
            System.out.println("Base coding was incorrectly set in Allele. Program quits.");
            System.exit(1);
        }
    }

    /**
     * Return the base of the allele.
     * @return
     */
    public char getAlleleBase () {
        return AlleleEncoder.getAlleleBaseFromCoding(this.alleleCoding);
    }

    /**
     * Return the coding of the allele, see {@link AlleleEncoder}.
     * @return
     */
    public byte getAlleleCoding() {
        return this.alleleCoding;
    }

    /**
     * Return the byte code of allele feature.
     * @return
     */
    public byte getAlleleFeature () {
        return this.feature;
    }

    /**
     * Set allele type, see {@link AlleleType}.
     * @param at
     */
    public void setAlleleType (AlleleType at) {
        this.feature = (byte)(feature | at.getFeature());
    }

    /**
     * Set the allele feature at once, see {@link AlleleType}.
     * @param feature
     */
    public void setAlleleFeature (byte feature) {
        this.feature = feature;
    }

    /**
     * Remove allele type, see {@link AlleleType}.
     * @param at
     */
    public void removeAlleleType (AlleleType at) {
        this.feature = (byte) (feature & (~at.getFeature()));
    }

    /**
     * Reset all allele types to false (remove all allele types), see {@link AlleleType}.
     */
    public void resetAlleleTypeToDefault () {
        this.feature = 0;
    }

    /**
     * Return if the allele is a specific type, see {@link AlleleType}.
     * @param at
     * @return
     */
    public boolean isAlleleTypeOf (AlleleType at) {
        if ((this.feature & at.getFeature()) == 0) return false;
        return true;
    }
}
