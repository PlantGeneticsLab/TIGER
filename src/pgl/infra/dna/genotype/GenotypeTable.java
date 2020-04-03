package pgl.infra.dna.genotype;

import java.nio.ByteBuffer;

/**
 * Interface for methods of a genotype table
 * @author feilu
 */
public interface GenotypeTable {

    /**
     * Return the number of taxa
     * @return
     */
    public int getTaxaNumber ();

    /**
     * Return the number of all sites
     * @return
     */
    public int getSiteNumber ();

    /**
     * Return the taxon name
     * @param taxonIndex
     * @return
     */
    public String getTaxonName (int taxonIndex);

    /**
     * Return the taxa names
     * @return
     */
    public String[] getTaxaNames ();

    /**
     * Return the chromosome of a specific site
     * @param siteIndex
     * @return
     */
    public short getChromosome (int siteIndex);

    /**
     * Return the position of a specific site
     * @param siteIndex
     * @return
     */
    public int getPosition (int siteIndex);

    /**
     * Sort by the chromosome position of site, which means by row
     */
    public void sortBySite ();

    /**
     * Sort by taxa names, which means by column
     */
    public void sortByTaxa ();

    /**
     * Return the index of a specific taxon
     * @param taxon
     * @return negative value if the taxon does not exist
     */
    public int getTaxonIndex (String taxon);

    /**
     * Return the index of a specific site
     * @param chromosome
     * @param position
     * @return negative value if the site does not exist
     */
    public int getSiteIndex (short chromosome, int position);

    /**
     * Return the byte value of a specific genotype
     * @param siteIndex
     * @param taxonIndex
     * @return
     */
    public byte getGenotypeByte(int siteIndex, int taxonIndex);

    /**
     * Return if a specific genotype is heterozygous
     * @param siteIndex
     * @param taxonIndex
     * @return
     */
    public boolean isHeterozygous (int siteIndex, int taxonIndex);

    /**
     * Return if a specific genotype is homozygous
     * @param siteIndex
     * @param taxonIndex
     * @return
     */
    public boolean isHomozygous (int siteIndex, int taxonIndex);

    /**
     * Return if a specific genotype is missing
     * @param siteIndex
     * @param taxonIndex
     * @return
     */
    public boolean isMissing (int siteIndex, int taxonIndex);

    /**
     * Return the total number of missing genotype at a specific site
     * @param siteIndex
     * @return
     */
    public int getMissingNumberBySite (int siteIndex);

    /**
     * Return the total number of missing genotype in a specific taxon
     * @param taxonIndex
     * @return
     */
    public int getMissingNumberByTaxon (int taxonIndex);

    /**
     * Return the total number of non-missing genotype at a specific site
     * @param siteIndex
     * @return
     */
    public int getNonMissingNumberBySite (int siteIndex);

    /**
     * Return the total number of non-missing in a specific taxon
     * @param taxonIndex
     * @return
     */
    public int getNonMissingNumberByTaxon (int taxonIndex);

    /**
     * Return the total number of homozygous genotype at a specific site
     * @param siteIndex
     * @return
     */
    public int getHomozygoteNumberBySite (int siteIndex);

    /**
     * Return the total number of homozygous genotype in a specific taxon
     * @param taxonIndex
     * @return
     */
    public int getHomozygoteNumberByTaxon (int taxonIndex);

    /**
     * Return the total number of heterozygous genotype at a specific site
     * @param siteIndex
     * @return
     */
    public int getHeterozygoteNumberBySite (int siteIndex);

    /**
     * Return the total number of heterzygous genotype in a specific taxon
     * @param taxonIndex
     * @return
     */
    public int getHeterozygoteNumberByTaxon (int taxonIndex);

    /**
     * Return the heterozygosity of a specific taxon
     * @param taxonIndex
     * @return
     */
    public float getTaxonHeterozygosity (int taxonIndex);

    /**
     * Return the fraction of heterozygous genotype at a specific site
     * @param siteIndex
     * @return
     */
    public double getSiteHeterozygoteFraction (int siteIndex);

    /**
     * Return byte value of the minor allele
     * @param siteIndex
     * @return
     */
    public byte getMinorAlleleByte(int siteIndex);

    /**
     * Return base of the minor allele
     * @param siteIndex
     * @return
     */
    public char getMinorAlleleBase(int siteIndex);

    /**
     * Return the minor allele frequency
     * @param siteIndex
     * @return
     */
    public double getMinorAlleleFrequency (int siteIndex);

    /**
     * Return the byte value of the major allele
     * @param siteIndex
     * @return
     */
    public byte getMajorAlleleByte(int siteIndex);

    /**
     * Return base of the major allele
     * @param siteIndex
     * @return
     */
    public char getMajorAlleleBase(int siteIndex);

    /**
     * Return the major allele frequency
     * @param siteIndex
     * @return
     */
    public float getMajorAlleleFrequency (int siteIndex);

    /**
     * Return the byte value of the reference allele
     * @param siteIndex
     * @return
     */
    public byte getReferenceAlleleByte(int siteIndex);

    /**
     * Return the base of the reference allele
     * @param siteIndex
     * @return
     */
    public char getReferenceAlleleBase(int siteIndex);

    /**
     * Return the reference allele frequency
     * @param siteIndex
     * @return
     */
    public double getReferenceAlleleFrequency (int siteIndex);

    /**
     * Return the byte value of the alternative allele
     * @param siteIndex
     * @return
     */
    public byte getAlternativeAlleleByte(int siteIndex);

    /**
     * Return base of the alternative allele
     * @param siteIndex
     * @return
     */
    public char getAlternativeAlleleBase(int siteIndex);

    /**
     * Return the alternative allele frequency
     * @param siteIndex
     * @return
     */
    public double getAlternativeAlleleFrequency (int siteIndex);

    /**
     * Return the start index of a chromosome, inclusive
     * @param chromosome
     * @return -1 if the chromosome does not exist
     */
    public int getStartIndexOfChromosome (short chromosome);

    /**
     * Return the end index of a chromosome, exclusive
     * @param chromosome
     * @return -1 if chromosome does not exist
     */
    public int getEndIndexOfChromosome (short chromosome);

    /**
     * Return the genetic divergence of two taxa, Dxy is defined as 1 - IBS (Identify by state)
     * @param taxonIndex1
     * @param taxonIndex2
     * @return Double.NaN if no shared non-missing sites exist
     */
    public double getDxy (int taxonIndex1, int taxonIndex2);

    /**
     * Return the genetic divergence of two taxa in a specified region, Dxy is defined as 1 - IBS (Identify by state)
     * @param taxonIndex1
     * @param taxonIndex2
     * @param startSiteIndex inclusive
     * @param endSiteIndex exclusive
     * @return Double.NaN if no shared non-missing sites exist
     */
    public double getDxy (int taxonIndex1, int taxonIndex2, int startSiteIndex, int endSiteIndex);

    /**
     * Return the genetic divergence of two taxa based on a list of sites, Dxy is defined as 1 - IBS (Identify by state)
     * @param taxonIndex1
     * @param taxonIndex2
     * @param siteIndices
     * @return Double.NaN if no shared non-missing sites exist
     */
    public double getDxy (int taxonIndex1, int taxonIndex2, int[] siteIndices);

    /**
     * Return a matrix of genetic divergence of all taxa based on all sites
     * @return
     */
    public double[][] getDxyMatrix ();

    /**
     * Return a matrix of genetic divergence of all taxa based on 10K evenly distributed sites
     * When total site number > 10K, size = 10K; or size = total site number
     * @return
     */
    public double[][] getDxyMatrixFast10K ();

    /**
     * Return a matrix of genetic divergence of all taxa based selected sites
     * @param siteIndices
     * @return
     */
    public double[][] getDxyMatrix (int[] siteIndices);
    /**
     * Return a genotype table by sub-setting the current table by sites
     * @param siteIndices
     * @return
     */
    public GenotypeTable getSubGenotypeTableBySite (int[] siteIndices);

    /**
     * Return a genotype table by sub-setting the current table by taxa
     * @param taxaIndices
     * @return
     */
    public GenotypeTable getSubGenotypeTableByTaxa (int[] taxaIndices);

    /**
     * Return a VCF genotype record at a specific site
     * @param siteIndex
     * @return
     */
    public String getUnphasedVCFRecord(int siteIndex);

    /**
     * Return a binary genotype output at a specific site, see {@link ByteBuffer}
     * @param index
     * @return
     */
    ByteBuffer getBinaryOutput(int index, ByteBuffer bb);
}
