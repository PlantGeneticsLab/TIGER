package pgl.infra.dna.genot;

/**
 * The file format standards of genotype table
 */
public enum GenoIOFormat {
    /**
     * VCF format
     */
    VCF,
    /**
     * VCF format compressed in gz
     */
    VCF_GZ,
    /**
     * Binary format
     */
    Binary,
    /**
     * Binary format compressed in gz
     */
    Binary_GZ,
    /**
     * HDF5 format
     */
    HDF5;
}
