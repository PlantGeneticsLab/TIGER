package pgl.infra.dna.genot;

import pgl.infra.dna.allele.AlleleEncoder;
import pgl.infra.dna.allele.AlleleType;
import pgl.infra.dna.snp.BiSNP;
import pgl.infra.utils.PStringUtils;

import java.nio.ByteBuffer;
import java.util.BitSet;
import java.util.List;

/**
 * Class holding SNP and genotype information per site,
 * <p>
 * Supports only bi-allelic SNPs, 3rd+ allele will be ignored. Allele depth is ignored.
 *
 * @author feilu
 */
public class SiteGenotype extends BiSNP {
    //Bit set of the 1st homologous chromosome, 1 is alt, 0 is ref
    BitSet phase1 = null;
    //Bit set of the 2nd homologous chromosome, 1 is alt, 0 is ref
    BitSet phase2 = null;
    //Bit set of the missing genotype at a specific site, 1 is missing
    BitSet missing = null;
    //Taxa number coding, max taxa capacity is 65536
    short taxaNumber = Short.MIN_VALUE;
    //Minor allele frequency
    float maf = Float.MIN_VALUE;

    /**
     * Is used to build VCF record, not thread safe
     */
    static StringBuilder vsb = new StringBuilder();

    /**
     * Construct an object
     */
    public SiteGenotype() {
        
    }

    /**
     * Construct an object by initializing all fields
     * @param chr
     * @param pos
     * @param refBase
     * @param altBase
     * @param info
     * @param phase1
     * @param phase2
     * @param missing
     * @param taxaNumber
     */
    public SiteGenotype(short chr, int pos, char refBase, char altBase, String info, BitSet phase1, BitSet phase2, BitSet missing, int taxaNumber) {
        super(chr, pos, refBase, altBase, info);
        this.phase1 = phase1;
        this.phase2 = phase2;
        this.missing = missing;
        this.taxaNumber = (short)(taxaNumber+Short.MIN_VALUE);
    }

    /**
     * Sort by taxa based on a given order
     * @param orderIndices
     */
    public void sortByTaxa (int[] orderIndices) {
        BitSet p1 = new BitSet(this.getTaxaNumber());
        BitSet p2 = new BitSet(this.getTaxaNumber());
        BitSet m = new BitSet(this.getTaxaNumber());
        for (int i = 0; i < this.getTaxaNumber(); i++) {
            if (isPhase1Alternative(orderIndices[i])) p1.set(i);
            if (isPhase2Alternative(orderIndices[i])) p2.set(i);
            if (isMissing(orderIndices[i])) m.set(i);
        }
        this.phase1 = p1; p1 = null;
        this.phase2 = p2; p2 = null;
        this.missing = m; m = null;
    }

    /**
     * Return a site genotype {@link SiteGenotype} by selecting taxa
     * @param taxaIndices
     * @return
     */
    public SiteGenotype getSubGenotypeByTaxa (int[] taxaIndices) {
        BitSet p1 = new BitSet(taxaIndices.length);
        BitSet p2 = new BitSet(taxaIndices.length);
        BitSet m = new BitSet(taxaIndices.length);
        for (int i = 0; i < taxaIndices.length; i++) {
            if (isPhase1Alternative(taxaIndices[i])) p1.set(i);
            if (isPhase2Alternative(taxaIndices[i])) p2.set(i);
            if (isMissing(taxaIndices[i])) m.set(i);
        }
        SiteGenotype sgb = new SiteGenotype(this.getChromosome(),this.getPosition(),AlleleEncoder.getAlleleBaseFromCoding(this.getReferenceAlleleByte()),
                                AlleleEncoder.getAlleleBaseFromCoding(this.getAlternativeAlleleByte()), null, p1, p2, m, taxaIndices.length);
        return sgb;
    }

    /**
     * Return the number of taxa
     * @return
     */
    public int getTaxaNumber () {
        return this.taxaNumber-Short.MIN_VALUE;
    }

    /**
     * Return the byte value of a specific genotype
     * @param taxonIndex
     * @return
     */
    public byte getGenotypeByte (int taxonIndex) {
        if (isMissing(taxonIndex)) return AlleleEncoder.genotypeMissingCoding;
        byte ref = this.getReferenceAlleleByte();
        byte alt = this.getAlternativeAlleleByte();
        byte b1 = AlleleEncoder.alleleMissingCoding;
        byte b2 = AlleleEncoder.alleleMissingCoding;
        if (isPhase1Alternative(taxonIndex)) b1 = alt;
        else b1 = ref;
        if (isPhase2Alternative(taxonIndex)) b2 = alt;
        else b2 = ref;
        return AlleleEncoder.getGenotypeCoding(b1, b2);
    }

    /**
     * Return if a specific genotype is missing
     * @param taxonIndex
     * @return
     */
    public boolean isMissing (int taxonIndex) {
        if (missing.get(taxonIndex)) return true;
        return false;
    }

    /**
     * Return if a specific genotype is heterozygous
     * @param taxonIndex
     * @return
     */
    public boolean isHeterozygous (int taxonIndex) {
        if (this.isMissing(taxonIndex)) return false;
        if (isPhase1Alternative(taxonIndex) == isPhase2Alternative(taxonIndex)) return false;
        return true;
    }

    /**
     * Return if specific genotype is homozygous
     * @param taxonIndex
     * @return
     */
    public boolean isHomozygous (int taxonIndex) {
        if (this.isMissing(taxonIndex)) return false;
        if (isPhase1Alternative(taxonIndex) == isPhase2Alternative(taxonIndex)) return true;
        return false;
    }

    /**
     * Return the number of missing genotype
     * @return
     */
    public int getMissingNumber () {
        return missing.cardinality();
    }

    /**
     * Return the number of non-missing genotype
     * @return
     */
    public int getNonMissingNumber () {
        return this.getTaxaNumber()-this.getMissingNumber();
    }

    /**
     * Return if the allele on phase 1 chromosome is alternative
     * @param taxonIndex
     * @return
     */
    public boolean isPhase1Alternative (int taxonIndex) {
        if (this.isMissing(taxonIndex)) return false;
        return phase1.get(taxonIndex);
    }

    /**
     * Return if the allele on phase 2 chromosome is alternative
     * @param taxonIndex
     * @return
     */
    public boolean isPhase2Alternative (int taxonIndex) {
        if (this.isMissing(taxonIndex)) return false;
        return phase2.get(taxonIndex);
    }

    /**
     * Return if the allele on phase 1 chromosome is reference
     * @param taxonIndex
     * @return
     */
    public boolean isPhase1Reference (int taxonIndex) {
        if (this.isMissing(taxonIndex)) return false;
        if (this.isPhase1Alternative(taxonIndex)) return false;
        return true;
    }

    /**
     * Return if the allele on phase 1 chromosome is reference
     * @param taxonIndex
     * @return
     */
    public boolean isPhase2Reference (int taxonIndex) {
        if (this.isMissing(taxonIndex)) return false;
        if (this.isPhase2Alternative(taxonIndex)) return false;
        return true;
    }

    /**
     * Return the number of heterozygous genotype across all taxa
     * @return
     */
    public int getHeterozygoteNumber () {
        BitSet phase1C = phase1.get(0, phase1.length());
        phase1C.xor(phase2);
        return phase1C.cardinality();
    }

    /**
     * Return the number of homozygous genotype across all taxa
     * @return
     */
    public int getHomozygoteNumber () {
        return this.getNonMissingNumber()-this.getHeterozygoteNumber();
    }

    /**
     * Return the number of alternative allele across all taxa
     * @return
     */
    public int getAlternativeAlleleNumber () {
        return this.phase1.cardinality()+this.phase2.cardinality();
    }

    /**
     * Return the alternative allele frequency
     * @return
     */
    public float getAlternativeAlleleFrequency () {
        if (this.isAlternativeAlleleTypeOf(AlleleType.Minor)) return this.maf;
        if (this.isAlternativeAlleleTypeOf(AlleleType.Major)) return (float)(1-this.maf);
        return (float)((double)this.getAlternativeAlleleNumber()/(this.getNonMissingNumber()*2));
    }

    /**
     * Return the number of reference allele across all taxa
     * @return
     */
    public int getReferenceAlleleNumber () {
        return this.getNonMissingNumber()*2-this.getAlternativeAlleleNumber();
    }

    /**
     * Return the reference allele frequency
     * @return
     */
    public float getReferenceAlleleFrequency () {
        if (this.isReferenceAlleleTypeOf(AlleleType.Minor)) return this.maf;
        if (this.isReferenceAlleleTypeOf(AlleleType.Major)) return (float)(1-this.maf);
        return (float)((double)this.getReferenceAlleleNumber()/(this.getNonMissingNumber()*2));
    }

    /**
     * Return the byte value of minor allele
     * @return
     */
    public byte getMinorAlleleByte () {
        if (this.reference.isAlleleTypeOf(AlleleType.Minor)) return this.getReferenceAlleleByte();
        if (this.alternative.isAlleleTypeOf(AlleleType.Minor)) return this.getAlternativeAlleleByte();
        this.getMinorAlleleFrequency();
        return this.getMinorAlleleByte();
    }

    /**
     * Return the number of minor allele across all taxa
     * @return
     */
    public int getMinorAlleleNumber () {
        if (this.reference.isAlleleTypeOf(AlleleType.Minor)) return this.getReferenceAlleleNumber();
        if (this.alternative.isAlleleTypeOf(AlleleType.Minor)) return this.getAlternativeAlleleNumber();
        this.getMinorAlleleFrequency();
        return this.getMinorAlleleNumber();
    }

    /**
     * Return the minor allele frequency
     * @return
     */
    public float getMinorAlleleFrequency () {
        if (this.maf != Float.MIN_VALUE) return this.maf;
        float altFre = this.getAlternativeAlleleFrequency();
        if (altFre < 0.5) {
            this.setAlternativeAlleleType(AlleleType.Minor);
            this.setReferenceAlleleType(AlleleType.Major);
        }
        else {
            this.setReferenceAlleleType(AlleleType.Minor);
            this.setAlternativeAlleleType(AlleleType.Major);
            altFre = (float)(1 - altFre);
        }
        this.maf = altFre;
        return altFre;
    }

    /**
     * Return the byte value of major allele
     * @return
     */
    public byte getMajorAlleleByte () {
        if (this.reference.isAlleleTypeOf(AlleleType.Major)) return this.getReferenceAlleleByte();
        if (this.alternative.isAlleleTypeOf(AlleleType.Major)) return this.getAlternativeAlleleByte();
        this.getMinorAlleleFrequency();
        return this.getMajorAlleleByte();
    }

    /**
     * Return the number of major allele across all taxa
     * @return
     */
    public int getMajorAlleleNumber () {
        if (this.reference.isAlleleTypeOf(AlleleType.Major)) return this.getReferenceAlleleNumber();
        if (this.alternative.isAlleleTypeOf(AlleleType.Major)) return this.getAlternativeAlleleNumber();
        this.getMinorAlleleFrequency();
        return this.getMajorAlleleNumber();
    }

    /**
     * Return the major allele frequency
     * @return
     */
    public float getMajorAlleleFrequency () {
        if (this.maf != Float.MIN_VALUE) return (float)(1-this.maf);
        return (float)(1-this.getMinorAlleleFrequency());
    }

    /**
     * Return unphased genotype record of the current site in VCF format
     * @return
     */
    public String getUnphasedVCFOutput() {
        return this.getVCFOutput('/');
    }

    private String getVCFOutput(char delimiter) {
        vsb.setLength(0);
        vsb.append(this.getChromosome()).append("\t").append(this.getPosition()).append("\t").append(this.getChromosome()).append("-").append(this.getPosition()).append("\t");
        vsb.append(this.getReferenceAlleleBase()).append("\t").append(this.getAlternativeAlleleBase()).append("\t.\t.\t");
        if (info == null) vsb.append(".");
        else vsb.append(info);
        vsb.append("\t").append("GT");
        for (int i = 0; i < this.getTaxaNumber(); i++) {
            if (isMissing(i)) vsb.append("\t").append(".").append(delimiter).append(".");
            else {
                vsb.append("\t");
                if (isPhase1Alternative(i)) vsb.append("1");
                else vsb.append("0");
                vsb.append(delimiter);
                if (isPhase2Alternative(i)) vsb.append("1");
                else vsb.append("0");
            }
        }
        return vsb.toString();
    }

    /**
     * Return the binary data block of genotype from the current site
     * @param bb
     * @return
     */
    public ByteBuffer getBinaryOutput (ByteBuffer bb) {
        bb.putShort(this.getChromosome());
        bb.putInt(this.getPosition());
        bb.put(AlleleEncoder.getGenotypeCoding(this.getReferenceAlleleByte(), this.getAlternativeAlleleByte()));
        bb.put(this.getReferenceAlleleFeature());
        bb.put(this.getAlternativeAlleleFeature());
        int size = (bb.capacity()-GenotypeExport.getByteSizeOfSNPInBinary())/3;
        byte[] b = new byte[size];
        byte[] ob = phase1.toByteArray();
        System.arraycopy(ob,0,b,0,ob.length);
        bb.put(b);
        b = new byte[size];
        ob = phase2.toByteArray();
        System.arraycopy(ob,0,b,0,ob.length);
        bb.put(b);
        b = new byte[size];
        ob = missing.toByteArray();
        System.arraycopy(ob,0,b,0,ob.length);
        bb.put(b);
        return bb;
    }

    /**
     * Static method build an object from a binary data block of genotype of the current site
     * @param bb
     * @param taxaNumber
     * @return
     */
    public static SiteGenotype buildFromBinaryLine (ByteBuffer bb, int taxaNumber) {
        //short chr, int pos, char refBase, char altBase, String info, BitSet phase1, BitSet phase2, BitSet missing, int taxaNumber
        bb.flip();
        short chr = bb.getShort();
        int pos = bb.getInt();
        byte geno = bb.get();
        char refBase = AlleleEncoder.getAlleleBase1FromGenotypeCoding(geno);
        char altBase = AlleleEncoder.getAlleleBase2FromGenotypeCoding(geno);
        byte refFeature = bb.get();
        byte altFeature = bb.get();
        int size = (bb.capacity()-GenotypeExport.getByteSizeOfSNPInBinary())/3;
        byte[] ba = new byte[size];
        bb.get(ba);
        BitSet p1 = BitSet.valueOf(ba);
        ba = new byte[size];
        bb.get(ba);
        BitSet p2 = BitSet.valueOf(ba);
        ba = new byte[size];
        bb.get(ba);
        BitSet m = BitSet.valueOf(ba);
        SiteGenotype sgb = new SiteGenotype(chr, pos, refBase, altBase, null, p1, p2, m, taxaNumber);
        sgb.setReferenceAlleleFeature(refFeature);
        sgb.setAlternativeAlleleFeature(altFeature);
        bb.clear();
        return sgb;
    }

    /**
     * Build and return an object of {@link SiteGenotype} from line of VCF format
     * @param line
     * @return
     */
    public static SiteGenotype buildFromVCFLine (String line) {
        List<String> l = PStringUtils.fastSplit(line);
        List<String> ll = null;
        String current = null;
        short chr = Short.parseShort(l.get(0));
        int pos = Integer.parseInt(l.get(1));
        char refBase = l.get(3).charAt(0);
        char altBase = l.get(4).charAt(0);
        String info = l.get(7);
        int taxaNumber = l.size()-9;
        BitSet phase1 = new BitSet(taxaNumber);
        BitSet phase2 = new BitSet(taxaNumber);
        BitSet missingP = new BitSet(taxaNumber);
        byte[] values = null;
        for (int i = 0; i < taxaNumber; i++) {
            current = l.get(i+9);
            if (current.startsWith(".")) {
                missingP.set(i);
                continue;
            }
            ll = PStringUtils.fastSplit(current, ":");
            values = ll.get(0).getBytes();
            if (values[0] == 49) {
                phase1.set(i);
            }
            if (values[2] == 49) {
                phase2.set(i);
            }
        }
        //SiteGenotypeBit sgb = new SiteGenotypeBit(chr, pos, refBase, altBase, info, majorP, minorP, missingP, taxaNumber);
        SiteGenotype sgb = new SiteGenotype(chr, pos, refBase, altBase, null, phase1, phase2, missingP, taxaNumber);
        return sgb;
    }
}
