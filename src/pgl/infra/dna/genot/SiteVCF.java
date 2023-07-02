package pgl.infra.dna.genot;

import pgl.infra.dna.allele.AlleleEncoder;
import pgl.infra.utils.PStringUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * Class holding SNP and VCF genotype information per site.
 * <p>
 * It is designed to be used in VCF filtering, taxon insertion/deletion/merging
 * @author feilu
 */
public class SiteVCF {
    short chr = Short.MIN_VALUE;
    int pos = Integer.MIN_VALUE;
    byte[] alleles = null;
    //if the site is missing, the depth is set to Short.Min_Value
    List<short[]> alleleDepthList = null;
    //if the site is missing, the genotype is set to Byte.Min_Value
    List<byte[]> genoList = null;

    public SiteVCF (String vcfRecord) {
        this.parseVCFRecord(vcfRecord);
    }

    public void parseVCFRecord (String vcfRecord) {
        List<String> l = PStringUtils.fastSplit(vcfRecord);
        List<String> ll = null;
        List<String> lll = null;
        this.chr = Short.parseShort(l.get(0));
        this.pos = Integer.parseInt(l.get(1));
        ll = PStringUtils.fastSplit(l.get(4), ",");
        int alleleNumber = 1+ll.size();
        alleles = new byte[alleleNumber];
        alleles[0] = AlleleEncoder.getAlleleCodingFromBase(l.get(3).charAt(0));
        for (int i = 0; i < ll.size(); i++) {
            alleles[i+1] = AlleleEncoder.getAlleleCodingFromBase(ll.get(i).charAt(0));
        }
        this.alleleDepthList = new ArrayList<>();
        this.genoList = new ArrayList<>();
        for (int i = 0; i < l.size()-9; i++) {
            byte[] geno = new byte[2];
            short[] depth = new short[alleles.length];
            if (!l.get(i+9).startsWith(".")) {
                ll = PStringUtils.fastSplit(l.get(i+9),":");
                if (ll.get(0).contains("/")) {
                    lll = PStringUtils.fastSplit(ll.get(0), "/");
                }
                else {
                    lll = PStringUtils.fastSplit(ll.get(0), "|");
                }
                for (int j = 0; j < geno.length; j++) {
                    geno[j] = Byte.parseByte(lll.get(j));
                }
                lll = PStringUtils.fastSplit(ll.get(1), ",");

                for (int j = 0; j < alleles.length; j++) {
                    depth[j] = Short.parseShort(lll.get(j));
                }
            }
            else {
                for (int j = 0; j < geno.length; j++) {
                    geno[j] = Byte.MIN_VALUE;
                }
                for (int j = 0; j < alleles.length; j++) {
                    depth[j] = Short.MIN_VALUE;
                }
            }
            this.genoList.add(geno);
            this.alleleDepthList.add(depth);
        }
    }

    /**
     * Return taxa number of the VCF record
     * @return
     */
    public int getTaxaNumber () {
        return this.genoList.size();
    }

    /**
     * Return the number of alleles, including Ref and Alts
     * @return
     */
    public int getAlleleNumber () {
        return this.alleles.length;
    }

    /**
     * Return the genotype as alleles in the listed order, e.g. 0 is the ref allele, 1 is first alt allele, 2 is the second alt allele.
     * Missing is denoted as Byte.Min_Value
     * @param taxonIndex
     * @return
     */
    public byte[] getGenotype (int taxonIndex) {
        return this.genoList.get(taxonIndex);
    }

    /**
     * Return the read depth of alleles
     * Missing denoted as Short.Min_Value
     * @param taxonIndex
     * @return
     */
    public short[] getAlleleDepth (int taxonIndex) {
        return this.alleleDepthList.get(taxonIndex);
    }

    /**
     * Return the number of missing genotype across all taxa
     * @return
     */
    public int getMissingTaxaNumber () {
        int cnt = 0;
        for (int i = 0; i < this.getTaxaNumber(); i++) {
            if (this.genoList.get(i)[0] < 0) {
                cnt++;
            }
        }
        return cnt;
    }

    /**
     * Return the number of non-missing genotype across all taxa
     * @return
     */
    public int getNonMissingTaxaNumber () {
        return this.getTaxaNumber()-this.getMissingTaxaNumber();
    }

    /**
     * Return the total read depth of each allele in the listed order
     * @return
     */
    public int[] getTotalAlleleDepth () {
        int[] totalCounts = new int[this.getAlleleNumber()];
        for (int i = 0; i < this.getTaxaNumber(); i++) {
            for (int j = 0; j < this.alleleDepthList.get(i).length; j++) {
                if (this.alleleDepthList.get(i)[j] < 0) continue;
                totalCounts[j]+=this.alleleDepthList.get(i)[j];
            }
        }
        return totalCounts;
    }

    /**
     * Return the total depth of all alleles
     * @return
     */
    public int getTotalDepth () {
        int[] alleleCounts = this.getTotalAlleleDepth();
        int cnt = 0;
        for (int i = 0; i < alleleCounts.length; i++) {
            cnt+=alleleCounts[i];
        }
        return cnt;
    }

    /**
     * Return the number of presence of each allele across all taxa
     * @return
     */
    public int[] getTotalAllelePresence () {
        int[] totalCounts = new int[this.getAlleleNumber()];
        for (int i = 0; i < this.getTaxaNumber(); i++) {
            for (int j = 0; j < totalCounts.length; j++) {
                for (int k = 0; k < 2; k++) {
                    if (this.getGenotype(i)[k] == j) {
                        totalCounts[j]++;
                        break;
                    }
                }
            }
        }
        return totalCounts;
    }

    /**
     * Return the number of genotype, e.g. AA,AB,BB or AA,AB,AC,BB,BC,CC if 2 alt alleles
     * @return
     */
    public int[][] getGenotypeCount () {
        int[][] counts = new int[this.getAlleleNumber()][this.getAlleleNumber()];
        for (int i = 0; i < this.getTaxaNumber(); i++) {
            if (this.getGenotype(i)[0] < 0) continue;
            counts[this.getGenotype(i)[0]][this.getGenotype(i)[1]]++;
        }
        return counts;
    }

    /**
     * Return the number of heterozygotes across all taxa
     * @return
     */
    public int getHeterozygotesNumber () {
        int cnt = 0;
        for (int i = 0; i < this.getTaxaNumber(); i++) {
            if (this.genoList.get(i)[0] < 0) continue;
            if (this.genoList.get(i)[0] == this.genoList.get(i)[1]) continue;
            cnt++;
        }
        return cnt;
    }

    /**
     * Return the minor allele frequency
     * @return
     */
    public float getMinorAlleleFrequency () {
        int cnt = 0;
        for (int i = 0; i < this.getTaxaNumber(); i++) {
            for (int j = 0; j < this.genoList.get(i).length; j++) {
                if (this.genoList.get(i)[j] == 1) cnt++;
            }
        }
        float maf = (float)((double)cnt/this.getNonMissingTaxaNumber()/2);
        if (maf > 0.5) return 1-maf;
        return maf;
    }

    /**
     * Delete the genotype of certain taxon
     * @param taxonIndex
     */
    public void deleteTaxon (int taxonIndex) {
        this.genoList.remove(taxonIndex);
        this.alleleDepthList.remove(taxonIndex);
    }

    /**
     * Set genotype of certain taxon
     * @param taxonIndex
     * @param geno
     */
    public void setGenotype (int taxonIndex, byte[] geno) {
        this.genoList.set(taxonIndex, geno);
    }

    /**
     * Set allele depth of certain taxon
     * @param taxonIndex
     * @param alleleDepth
     */
    public void setAlleleDepth (int taxonIndex, short[] alleleDepth) {
        this.alleleDepthList.set(taxonIndex, alleleDepth);
    }

    /**
     * Return a VCF record with combined error rate of 0.05 by default
     * @return
     */
    public String getVCFRecord () {
        return this.getVCFRecord(0.05);
    }

    /**
     * Return a VCF record
     * @param combinedErrorRate error rate used to calculate the likelihood of genotype
     * @return
     */
    public String getVCFRecord (double combinedErrorRate) {
        StringBuilder sb = new StringBuilder ();
        sb.append(this.chr).append("\t").append(this.pos).append("\t").append(this.chr).append("-").append(this.pos)
                .append("\t").append(AlleleEncoder.getAlleleBaseFromCoding(this.alleles[0])).append("\t");
        for (int i = 0; i < alleles.length-1; i++) {
            sb.append(AlleleEncoder.getAlleleBaseFromCoding(this.alleles[i+1])).append(",");
        }
        sb.deleteCharAt(sb.length()-1);
        sb.append("\t.\t.\t");
        sb.append("DP=").append(this.getTotalDepth()).append(";NZ=").append(this.getNonMissingTaxaNumber()).append(";AD=");
        int[] counts = this.getTotalAlleleDepth();
        for (int i = 0; i < counts.length; i++) {
            sb.append(counts[i]).append(",");
        }
        sb.deleteCharAt(sb.length()-1);
        sb.append(";AC=");
        counts = this.getTotalAllelePresence();
        for (int i = 0; i < counts.length; i++) {
            sb.append(counts[i]).append(",");
        }
        sb.deleteCharAt(sb.length()-1);
        sb.append(";GN=");
        int[][] cnts = this.getGenotypeCount();
        int hetCnt = 0;
        for (int i = 0; i < cnts.length; i++) {
            for (int j = i; j < cnts.length; j++) {
                sb.append(cnts[i][j]).append(",");
                if (i == j) continue;
                hetCnt+=cnts[i][j];
            }
        }
        sb.deleteCharAt(sb.length()-1);
        sb.append(";HT=").append(hetCnt);
        sb.append(";MAF=").append(this.getMinorAlleleFrequency()).append("\tGT:AD:GL");
        for (int i = 0; i < this.getTaxaNumber(); i++) {
            if (this.getAlleleDepth(i)[0] < 0) {
                sb.append("\t./.");
            }
            else {
                sb.append("\t").append(VCFUtils.getGenotypeByShort(this.getAlleleDepth(i), combinedErrorRate));
            }

        }
        return sb.toString();
    }
}
