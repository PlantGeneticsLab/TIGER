package pgl.infra.dna.allele;

/**
 * Hard coded standard of 8 allele types, representing ref, alt, 2nd alt, major, minor, 2nd minor, ancestral, derived.
 * @author feilu 
 */
public enum AlleleType {
    /**
     * Reference allele
     */
    Reference (0, (byte)0b1),
    /**
     * Alternative allele
     */
    Alternative (1, (byte)0b10),
    /**
     * The second alternative allele, if there is one
     */
    Alternative2 (2, (byte)0b100),
    /**
     * Major allele
     */
    Major (3, (byte)0b1000),
    /**
     * Minor allele
     */
    Minor (4, (byte)0b10000),
    /**
     * The second minor allele, if there is one
     */
    Minor2 (5, (byte)0b100000),
    /**
     * Ancestral allele
     */
    Ancestral (6, (byte)0b1000000),
    /**
     * Derived allele
     */
    Derived (7, (byte)0b10000000);

    private final int index;
    private final byte feature;

    AlleleType (int index, byte feature) {
        this.index = index;
        this.feature = feature;
    }

    /**
     * Return the number of allele type.
     * @return
     */
    public int getAlleleTypeNumber () {
        return AlleleType.values().length;
    }

    /**
     * Return feature value of a specific type.
     * @return
     */
    public byte getFeature () {
        return feature;
    }

    /**
     * Return the index of a specific type.
     * @return
     */
    public int getIndex () {
        return index;
    }
}
