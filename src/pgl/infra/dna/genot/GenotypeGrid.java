package pgl.infra.dna.genot;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import pgl.PGLConstraints;
import pgl.infra.dna.allele.AlleleEncoder;
import pgl.infra.dna.allele.AlleleType;
import pgl.infra.dna.snp.BiSNP;
import pgl.infra.pos.ChrPos;
import pgl.infra.utils.*;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;
import java.util.concurrent.*;

/**
 * Class holding genotype data. Genotypes are presented by both site and taxon in bitset.
 * <p>
 * Supports only bi-allelic SNPs, 3rd+ allele will be ignored. Allele depth is ignored.
 *
 * @author feilu
 */
public class GenotypeGrid implements GenotypeTable, Swapper, IntComparator {
    /**
     * The direction of the grid, either by site (Row of VCF) or by taxon (Column of VCF)
     */
    public static enum GridDirection {BySite, ByTaxon};

    /**
     * The taxa in the genotype table
     */
    String[] taxa = null;
    /**
     * SNPs of the genotype table
     */
    BiSNP[] snps = null;

    /**
     * Minor allele frequency of all SNPs
     */
    float[] mafs = null;

    /**
     * Bit genotype by site, the first dimension is site; the second dimension is phase 1, phase 2, and missing
     */
    BitSet[][] genoSite = null;

    /**
     * Bit genotype by taxon, the first dimension is taxon; the second dimension is phase 1, phase 2, and missing
     */
    BitSet[][] genoTaxon = null;


    /**
     * Construct an object by reading a file
     * @param infileS
     * @param format
     */
    public GenotypeGrid (String infileS, GenoIOFormat format) {
        if (format == GenoIOFormat.VCF) {
            this.buildFromVCF(infileS);
        }
        else if (format == GenoIOFormat.VCF_GZ) {
            this.buildFromVCF(infileS);
        }
        else if (format == GenoIOFormat.Binary) {
            this.buildFromBinary(infileS);
        }
        else if (format == GenoIOFormat.Binary_GZ) {
            this.buildFromBinary(infileS);
        }
        else if (format == GenoIOFormat.HDF5) {
            throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        }
        this.sortByTaxa();
    }

    public GenotypeGrid (BitSet[][] bArray, GridDirection gd, String[] taxa, BiSNP[] snps) {
        this.taxa = taxa;
        this.snps = snps;
        if (gd == GridDirection.BySite) {
            genoSite = bArray;
            this.transposeSiteToTaxon();
        }
        else if (gd == GridDirection.ByTaxon) {
            genoTaxon = bArray;
            this.transposeTaxonToSite();
        }
        this.mafs = new float[this.getSiteNumber()];
        Arrays.fill(mafs, Float.MIN_VALUE);
    }

    @Override
    public int getTaxaNumber() {
        return this.taxa.length;
    }

    @Override
    public int getSiteNumber() {
        return this.snps.length;
    }

    @Override
    public String getTaxonName(int taxonIndex) {
        return this.taxa[taxonIndex];
    }

    @Override
    public String[] getTaxaNames() {
//        String nTaxa[]  = new String[this.getTaxaNumber()];
//        System.arraycopy(taxa, 0, nTaxa, 0, nTaxa.length);
//        return nTaxa;
        return this.taxa;
    }

    @Override
    public short getChromosome(int siteIndex) {
        return this.snps[siteIndex].getChromosome();
    }

    @Override
    public int getPosition(int siteIndex) {
        return this.snps[siteIndex].getPosition();
    }

    @Override
    public void sortBySite() {
        System.out.println("Start sorting genotype table by site");
        long start = System.nanoTime();
        GenericSorting.quickSort(0, this.getSiteNumber(), this, this);
        StringBuilder sb = new StringBuilder();
        sb.append("Sorting finished in ").append(Benchmark.getTimeSpanSeconds(start)).append(" seconds.");
        System.out.println(sb.toString());
        this.transposeSiteToTaxon();
    }

    @Override
    public void sortByTaxa() {
        System.out.println("Start sorting genotype table by taxon");
        long start = System.nanoTime();
        int[] indices = PArrayUtils.getIndicesByAscendingValue(this.taxa);
        Arrays.sort(this.taxa);
        BitSet[][] nGenoTaxon = new BitSet[this.getTaxaNumber()][];
        for (int i = 0; i < this.getTaxaNumber(); i++) {
            nGenoTaxon[i] = this.genoTaxon[indices[i]];
        }
        this.genoTaxon = nGenoTaxon;
        StringBuilder sb = new StringBuilder();
        sb.append("Sorting finished in ").append(Benchmark.getTimeSpanSeconds(start)).append(" seconds.");
        System.out.println(sb.toString());
        this.transposeTaxonToSite();
    }

    @Override
    public int getTaxonIndex(String taxon) {
        return Arrays.binarySearch(taxa, taxon);
    }

    @Override
    public int getSiteIndex(short chromosome, int position) {
        ChrPos query = new ChrPos (chromosome, position);
        int index = Arrays.binarySearch(this.snps, query);
        return index;
    }

    @Override
    public byte getGenotypeByte(int siteIndex, int taxonIndex) {
        if (isMissing(siteIndex, taxonIndex)) return AlleleEncoder.genotypeMissingCoding;
        byte ref = this.getReferenceAlleleByte(siteIndex);
        byte alt = this.getAlternativeAlleleByte(siteIndex);
        byte b1 = AlleleEncoder.alleleMissingCoding;
        byte b2 = AlleleEncoder.alleleMissingCoding;
        if (isPhase1Alternative(siteIndex, taxonIndex)) b1 = alt;
        else b1 = ref;
        if (isPhase2Alternative(siteIndex, taxonIndex)) b2 = alt;
        else b2 = ref;
        return AlleleEncoder.getGenotypeCoding(b1, b2);
    }

    @Override
    public boolean isHeterozygous(int siteIndex, int taxonIndex) {
        if (this.isMissing(siteIndex,taxonIndex)) return false;
        if (this.isPhase1Alternative(siteIndex, taxonIndex) == this.isPhase2Alternative(siteIndex, taxonIndex)) return false;
        return true;
    }

    @Override
    public boolean isHomozygous(int siteIndex, int taxonIndex) {
        if (this.isMissing(siteIndex,taxonIndex)) return false;
        if (this.isPhase1Alternative(siteIndex, taxonIndex) == this.isPhase2Alternative(siteIndex, taxonIndex)) return true;
        return false;
    }

    @Override
    public boolean isMissing(int siteIndex, int taxonIndex) {
        if (this.genoSite[siteIndex][2].get(taxonIndex)) return true;
        return false;
    }

    @Override
    public boolean isPhase1Alternative(int siteIndex, int taxonIndex) {
        if (this.genoSite[siteIndex][0].get(taxonIndex)) return true;
        return false;
    }

    @Override
    public boolean isPhase2Alternative(int siteIndex, int taxonIndex) {
        if (this.genoSite[siteIndex][1].get(taxonIndex)) return true;
        return false;
    }

    @Override
    public boolean isPhase1Reference(int siteIndex, int taxonIndex) {
        if (this.isMissing(siteIndex, taxonIndex)) return false;
        if (this.isPhase1Alternative(siteIndex, taxonIndex)) return false;
        return true;
    }

    @Override
    public boolean isPhase2Reference(int siteIndex, int taxonIndex) {
        if (this.isMissing(siteIndex, taxonIndex)) return false;
        if (this.isPhase2Alternative(siteIndex, taxonIndex)) return false;
        return true;
    }

    @Override
    public boolean isAlternativeAlleleTypeOf(AlleleType at, int siteIndex) {
        return snps[siteIndex].isAlternativeAlleleTypeOf(at);
    }

    @Override
    public boolean isReferenceAlleleTypeOf(AlleleType at, int siteIndex) {
        return snps[siteIndex].isReferenceAlleleTypeOf(at);
    }

    @Override
    public int getMissingNumberBySite(int siteIndex) {
        return this.genoSite[siteIndex][2].cardinality();
    }

    @Override
    public int getMissingNumberByTaxon(int taxonIndex) {
        return this.genoTaxon[taxonIndex][2].cardinality();
    }

    @Override
    public int getNonMissingNumberBySite(int siteIndex) {
        return this.getTaxaNumber() - this.getMissingNumberBySite(siteIndex);
    }

    @Override
    public int getNonMissingNumberByTaxon(int taxonIndex) {
        return this.getSiteNumber() - this.getMissingNumberByTaxon(taxonIndex);
    }

    @Override
    public int getHomozygoteNumberBySite(int siteIndex) {
        return this.getNonMissingNumberBySite(siteIndex) - this.getHeterozygoteNumberBySite(siteIndex);
    }

    @Override
    public int getHomozygoteNumberByTaxon(int taxonIndex) {
        return this.getNonMissingNumberByTaxon(taxonIndex) - this.getHeterozygoteNumberByTaxon(taxonIndex);
    }

    @Override
    public int getHeterozygoteNumberBySite(int siteIndex) {
        BitSet phase1C = this.genoSite[siteIndex][0].get(0, this.getTaxaNumber());
        phase1C.xor(this.genoSite[siteIndex][1]);
        return phase1C.cardinality();
    }

    @Override
    public int getHeterozygoteNumberByTaxon(int taxonIndex) {
        BitSet phase1C = this.genoTaxon[taxonIndex][0].get(0, this.getSiteNumber());
        phase1C.xor(this.genoTaxon[taxonIndex][1]);
        return phase1C.cardinality();
    }

    @Override
    public int getAlternativeAlleleNumberBySite(int siteIndex) {
        int cnt = 0;
        for (int i = 0; i < 2; i++) {
            cnt+=this.genoSite[siteIndex][i].cardinality();
        }
        return cnt;
    }

    @Override
    public int getAlternativeAlleleOccurrenceBySite (int siteIndex) {
        BitSet bs = (BitSet)this.genoSite[siteIndex][0].clone();
        bs.or(this.genoSite[siteIndex][1]);
        return bs.cardinality();
    }

    @Override
    public float getHeterozygousProportionByTaxon(int taxonIndex) {
        return (float)((double)this.getHeterozygoteNumberByTaxon(taxonIndex)/this.getNonMissingNumberByTaxon(taxonIndex));
    }

    @Override
    public float getHeterozygousProportionBySite(int siteIndex) {
        return (float)((double)this.getHeterozygoteNumberBySite(siteIndex)/this.getNonMissingNumberBySite(siteIndex));
    }

    @Override
    public byte getMinorAlleleByte(int siteIndex) {
        if (this.snps[siteIndex].reference.isAlleleTypeOf(AlleleType.Minor)) return this.getReferenceAlleleByte(siteIndex);
        if (this.snps[siteIndex].alternative.isAlleleTypeOf(AlleleType.Minor)) return this.getAlternativeAlleleByte(siteIndex);
        this.getMinorAlleleFrequency(siteIndex);
        return this.getMinorAlleleByte(siteIndex);
    }

    @Override
    public char getMinorAlleleBase(int siteIndex) {
        if (this.snps[siteIndex].reference.isAlleleTypeOf(AlleleType.Minor)) return this.getReferenceAlleleBase(siteIndex);
        if (this.snps[siteIndex].alternative.isAlleleTypeOf(AlleleType.Minor)) return this.getAlternativeAlleleBase(siteIndex);
        this.getMinorAlleleFrequency(siteIndex);
        return this.getMinorAlleleBase(siteIndex);
    }

    @Override
    public float getMinorAlleleFrequency(int siteIndex) {
        if (this.mafs[siteIndex] != Float.MIN_VALUE) return this.mafs[siteIndex];
        float altFre = this.getAlternativeAlleleFrequency(siteIndex);
        if (altFre < 0.5) {
            this.setAlternativeAlleleType(AlleleType.Minor, siteIndex);
            this.setReferenceAlleleType(AlleleType.Major, siteIndex);
        }
        else {
            this.setReferenceAlleleType(AlleleType.Minor, siteIndex);
            this.setAlternativeAlleleType(AlleleType.Major, siteIndex);
            altFre = (float)(1 - altFre);
        }
        mafs[siteIndex] = altFre;
        return altFre;
    }

    @Override
    public byte getMajorAlleleByte(int siteIndex) {
        if (this.snps[siteIndex].reference.isAlleleTypeOf(AlleleType.Major)) return this.getReferenceAlleleByte(siteIndex);
        if (this.snps[siteIndex].alternative.isAlleleTypeOf(AlleleType.Major)) return this.getAlternativeAlleleByte(siteIndex);
        this.getMinorAlleleFrequency(siteIndex);
        return this.getMajorAlleleByte(siteIndex);
    }

    @Override
    public char getMajorAlleleBase(int siteIndex) {
        if (this.snps[siteIndex].reference.isAlleleTypeOf(AlleleType.Major)) return this.getReferenceAlleleBase(siteIndex);
        if (this.snps[siteIndex].alternative.isAlleleTypeOf(AlleleType.Major)) return this.getAlternativeAlleleBase(siteIndex);
        this.getMinorAlleleFrequency(siteIndex);
        return this.getMajorAlleleBase(siteIndex);
    }

    @Override
    public float getMajorAlleleFrequency(int siteIndex) {
        return 1-this.getMinorAlleleFrequency(siteIndex);
    }

    @Override
    public byte getReferenceAlleleByte(int siteIndex) {
        return snps[siteIndex].getReferenceAlleleByte();
    }

    @Override
    public char getReferenceAlleleBase(int siteIndex) {
        return snps[siteIndex].getReferenceAlleleBase();
    }

    @Override
    public float getReferenceAlleleFrequency(int siteIndex) {
        if (this.isReferenceAlleleTypeOf(AlleleType.Minor, siteIndex)) return this.mafs[siteIndex];
        if (this.isReferenceAlleleTypeOf(AlleleType.Major, siteIndex)) return 1-this.mafs[siteIndex];
        this.getMinorAlleleFrequency(siteIndex);
        return this.getReferenceAlleleFrequency(siteIndex);
    }

    @Override
    public byte getAlternativeAlleleByte(int siteIndex) {
        return snps[siteIndex].getAlternativeAlleleByte();
    }

    @Override
    public char getAlternativeAlleleBase(int siteIndex) {
        return snps[siteIndex].getAlternativeAlleleBase();
    }

    @Override
    public float getAlternativeAlleleFrequency(int siteIndex) {
        if (this.isAlternativeAlleleTypeOf(AlleleType.Minor, siteIndex)) return mafs[siteIndex];
        if (this.isAlternativeAlleleTypeOf(AlleleType.Major, siteIndex)) return (float)(1-mafs[siteIndex]);
        return (float)((double)this.getAlternativeAlleleNumberBySite(siteIndex)/(this.getNonMissingNumberBySite(siteIndex)*2));
    }

    @Override
    public int getStartIndexOfChromosome(short chromosome) {
        int index = this.getSiteIndex(chromosome, Integer.MIN_VALUE);
        if (index < 0) {
            index = -index - 1;
            if (index < this.getSiteNumber() && this.getChromosome(index) == chromosome) return index;
            return -1;
        }
        else {
            while (index > 0 && this.getChromosome(index-1) == chromosome) {
                index--;
            }
            return index;
        }
    }

    @Override
    public int getEndIndexOfChromosome(short chromosome) {
        int index = this.getSiteIndex(chromosome, Integer.MAX_VALUE);
        if (index < 0) {
            index = -index - 2;
            if (this.getChromosome(index) == chromosome) return index+1;
            else return -1;
        }
        else {
            while ((index+1) < this.getSiteNumber() && this.getChromosome(index+1) == chromosome) {
                index++;
            }
            return index+1;
        }
    }

    @Override
    public float getIBSDistance(int taxonIndex1, int taxonIndex2) {
        return this.getIBSDistance(taxonIndex1, taxonIndex2, 0, this.getSiteNumber());
    }

    @Override
    public float getIBSDistance(int taxonIndex1, int taxonIndex2, int startSiteIndex, int endSiteIndex) {
        int length = endSiteIndex - startSiteIndex;
        BitSet[] t1 = new BitSet[2];
        BitSet[] t2 = new BitSet[2];
        BitSet bothMissing = this.genoTaxon[taxonIndex1][2].get(startSiteIndex, endSiteIndex);
        bothMissing.or(this.genoTaxon[taxonIndex2][2].get(startSiteIndex, endSiteIndex));
        for (int i = 0; i < t1.length; i++) {
            t1[i] = this.genoTaxon[taxonIndex1][i].get(startSiteIndex, endSiteIndex);
            t1[i].or(bothMissing);
            t2[i] = this.genoTaxon[taxonIndex2][i].get(startSiteIndex, endSiteIndex);
            t2[i].or(bothMissing);
        }
        int cnt = 0;
        for (int i = 0; i < t1.length; i++) {
            for (int j = 0; j < t2.length; j++) {
                BitSet bs = (BitSet)t1[i].clone();
                bs.xor(t2[j]);
                cnt+=bs.cardinality();
            }
        }
        float dxy = (float) ((double)cnt/4/(endSiteIndex-startSiteIndex-bothMissing.cardinality()));
//        float dxy = (float) ((double)cnt/4/(endSiteIndex-startSiteIndex)); //incorrect calculation which is used in Tassel, didn't count missing data
        return dxy;
    }

    @Override
    public float getIBSDistance(int taxonIndex1, int taxonIndex2, int[] siteIndices) {
        BitSet selected = new BitSet(this.getSiteNumber());
        for (int i = 0; i < siteIndices.length; i++) {
            selected.set(siteIndices[i]);
        }
        int cnt = 0;
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                BitSet bs1 = this.genoTaxon[taxonIndex1][i].get(0, this.getSiteNumber());
                bs1.and(selected);
                BitSet bs2 = this.genoTaxon[taxonIndex2][i].get(0, this.getSiteNumber());
                bs2.and(selected);
                bs1.xor(bs2);;
                cnt+=bs1.cardinality();
            }
        }
        BitSet bs1 = this.genoTaxon[taxonIndex1][2].get(0, this.getSiteNumber());
        bs1.and(selected);
        BitSet bs2 = this.genoTaxon[taxonIndex2][2].get(0, this.getSiteNumber());
        bs2.and(selected);
        bs1.or(bs2);
        float dxy = (float) ((double)cnt/4/(siteIndices.length-bs1.cardinality()));
//        float dxy = (float) ((double)cnt/4/(endSiteIndex-startSiteIndex)); //incorrect calculation which is used in Tassel, didn't count missing data
        return dxy;
    }

    @Override
    public float[][] getIBSDistanceMatrix() {
        return this.getIBSDistanceMatrix(0, this.getSiteNumber());
    }

    @Override
    public float[][] getIBSDistanceMatrix(int startIndex, int endIndex) {
        float[][] matrix = new float[this.getTaxaNumber()][this.getTaxaNumber()];
        List<Integer> indexList = new ArrayList<>();
        for (int i = 0; i < matrix.length-1; i++) {
            indexList.add(i);
        }
        indexList.parallelStream().forEach(i -> {
            for (int j = i + 1; j < this.getTaxaNumber(); j++) {
                matrix[i][j] = this.getIBSDistance(i, j, startIndex, endIndex);
                matrix[j][i] = matrix[i][j];
            }
        });
        return matrix;
    }

    @Override
    public float[][] getIBSDistanceMatrix(int[] siteIndices) {
        float[][] matrix = new float[this.getTaxaNumber()][this.getTaxaNumber()];
        List<Integer> indexList = new ArrayList<>();
        for (int i = 0; i < matrix.length-1; i++) {
            indexList.add(i);
        }
        indexList.parallelStream().forEach(i -> {
            for (int j = i + 1; j < this.getTaxaNumber(); j++) {
                matrix[i][j] = this.getIBSDistance(i, j, siteIndices);
                matrix[j][i] = matrix[i][j];
            }
        });
        return matrix;
    }

    @Override
    public GenotypeTable getSubGenotypeTableBySite(int[] siteIndices) {
        BitSet[][] bArray = new BitSet[siteIndices.length][3];
        for (int i = 0; i < siteIndices.length; i++) {
            for (int j = 0; j < bArray[0].length; j++) {
                bArray[i][j] = this.genoSite[siteIndices[i]][j];
            }
        }
        BiSNP[] nsnps = new BiSNP[siteIndices.length];
        for (int i = 0; i < siteIndices.length; i++) {
            nsnps[i] = this.snps[siteIndices[i]];
        }
        return new GenotypeGrid(bArray, GridDirection.BySite, this.taxa, nsnps);
    }

    @Override
    public GenotypeTable getSubGenotypeTableByTaxa(int[] taxaIndices) {
        BitSet[][] bArray = new BitSet[taxaIndices.length][3];
        for (int i = 0; i < taxaIndices.length; i++) {
            for (int j = 0; j < bArray[0].length; j++) {
                bArray[i][j] = this.genoTaxon[taxaIndices[i]][j];
            }
        }
        String[] nTaxa = new String[taxaIndices.length];
        for (int i = 0; i < taxaIndices.length; i++) {
            nTaxa[i] = this.getTaxonName(taxaIndices[i]);
        }
        BiSNP[] nsnps = new BiSNP[this.getSiteNumber()];
        for (int i = 0; i < this.getSiteNumber(); i++) {
            nsnps[i] = this.snps[i].replicateWithoutFeature();
        }
        return new GenotypeGrid(bArray, GridDirection.ByTaxon, nTaxa, nsnps);
    }

    @Override
    public String getUnphasedVCFRecord(int siteIndex) {
        StringBuilder vsb = new StringBuilder();
        char delimiter = '/';
        vsb.setLength(0);
        vsb.append(this.getChromosome(siteIndex)).append("\t").append(this.getPosition(siteIndex)).append("\t").append(this.getChromosome(siteIndex)).append("-").append(this.getPosition(siteIndex)).append("\t");
        vsb.append(this.getReferenceAlleleBase(siteIndex)).append("\t").append(this.getAlternativeAlleleBase(siteIndex)).append("\t.\t.\t");
        if (this.snps[siteIndex].getSNPInfo() == null) vsb.append(".");
        else vsb.append(this.snps[siteIndex].getSNPInfo());
        vsb.append("\t").append("GT");
        for (int i = 0; i < this.getTaxaNumber(); i++) {
            if (isMissing(siteIndex, i)) vsb.append("\t").append(".").append(delimiter).append(".");
            else {
                vsb.append("\t");
                if (isPhase1Alternative(siteIndex, i)) vsb.append("1");
                else vsb.append("0");
                vsb.append(delimiter);
                if (isPhase2Alternative(siteIndex, i)) vsb.append("1");
                else vsb.append("0");
            }
        }
        return vsb.toString();
    }

    @Override
    public ByteBuffer getBinaryOutput(int siteIndex, ByteBuffer bb) {
        bb.putShort(this.getChromosome(siteIndex));
        bb.putInt(this.getPosition(siteIndex));
        bb.put(AlleleEncoder.getGenotypeCoding(this.getReferenceAlleleByte(siteIndex), this.getAlternativeAlleleByte(siteIndex)));
        bb.put(this.snps[siteIndex].getReferenceAlleleFeature());
        bb.put(this.snps[siteIndex].getAlternativeAlleleFeature());
        int size = (bb.capacity()-GenotypeExport.getByteSizeOfSNPInBinary())/3;
        for (int i = 0; i < this.genoSite[siteIndex].length; i++) {
            byte[] b = new byte[size];
            byte[] ob = this.genoSite[siteIndex][i].toByteArray();
            System.arraycopy(ob,0,b,0,ob.length);
            bb.put(b);
        }
        return bb;
    }

    @Override
    public void setAlternativeAlleleType(AlleleType at, int siteIndex) {
        this.snps[siteIndex].setAlternativeAlleleType(at);
    }

    @Override
    public void setReferenceAlleleType(AlleleType at, int siteIndex) {
        this.snps[siteIndex].setReferenceAlleleType(at);
    }

    /**
     * Reader of a binary genotype file
     * @param infileS
     */
    private void buildFromBinary (String infileS) {
        try{
            DataInputStream dis = null;
            if (infileS.endsWith(".gz")) {
                dis = IOUtils.getBinaryGzipReader(infileS);
            }
            else {
                dis = IOUtils.getBinaryReader(infileS);
            }
            int siteNumber = dis.readInt();
            short taxaNumber = (short)dis.readInt();
            this.taxa = new String[taxaNumber];
            for (int i = 0; i < taxaNumber; i++) {
                this.taxa[i] = dis.readUTF();
            }
            ExecutorService pool = Executors.newFixedThreadPool(PGLConstraints.parallelLevel);
            List<Future<GenoSiteBlockBinary>> resultList = new ArrayList<>();
            byte[][] input = new byte[GenoSiteBlockBinary.blockSize][GenotypeExport.getByteSizeOfSiteInBinary(taxaNumber)];
            int siteCount = 0;
            int startIndex = 0;
            int actBlockSize = 0;
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < siteNumber; i++) {
                dis.read(input[actBlockSize]);
                siteCount++;
                actBlockSize++;
                if (siteCount% GenoSiteBlockBinary.blockSize == 0) {
                    GenoSiteBlockBinary gsb = new GenoSiteBlockBinary(input, startIndex, actBlockSize);
                    Future<GenoSiteBlockBinary> result = pool.submit(gsb);
                    resultList.add(result);
                    startIndex+=actBlockSize;
                    actBlockSize = 0;
                    input = new byte[GenoSiteBlockBinary.blockSize][GenotypeExport.getByteSizeOfSiteInBinary(taxaNumber)];
                }
                if (siteCount%1000000 == 0) {
                    sb.setLength(0);
                    sb.append("Read in ").append(siteCount).append(" SNPs from ").append(infileS);
                    System.out.println(sb.toString());
                }
            }
            dis.close();
            if (actBlockSize != 0) {
                GenoSiteBlockBinary gsb = new GenoSiteBlockBinary(input, startIndex, actBlockSize);
                Future<GenoSiteBlockBinary> result = pool.submit(gsb);
                resultList.add(result);
            }
            pool.shutdown();
            pool.awaitTermination(Long.MAX_VALUE, TimeUnit.MICROSECONDS);
            this.genoSite = new BitSet[siteCount][];
            this.snps = new BiSNP[siteCount];
            for (int i = 0; i < resultList.size(); i++) {
                GenoSiteBlockBinary block = resultList.get(i).get();
                for (int j = 0; j < block.actBlockSize; j++) {
                    this.snps[block.getStartIndex()+j] = block.getSNPBlock()[j];
                    this.genoSite[block.getStartIndex()+j] = block.getGenoSiteBlock()[j];
                }
            }
            sb.setLength(0);
            sb.append("A total of ").append(this.getSiteNumber()).append(" SNPs are in ").append(infileS).append("\n");
            sb.append("Genotype table is successfully built");
            System.out.println(sb.toString());
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        this.transposeSiteToTaxon();
        mafs = new float[this.getSiteNumber()];
        Arrays.fill(mafs, Float.MIN_VALUE);
    }

    private void buildFromVCF (String infileS) {
        try {
            List<String> vcfAnnotationList = new ArrayList<>();
            BufferedReader br = null;
            if (infileS.endsWith(".gz")) {
                br = IOUtils.getTextGzipReader(infileS);
            }
            else {
                br = IOUtils.getTextReader(infileS);
            }
            String temp = null;
            while ((temp = br.readLine()).startsWith("##")) {
                vcfAnnotationList.add(temp);
            }
            List<String> l = new ArrayList<>();
            l = PStringUtils.fastSplit(temp);
            this.taxa = new String[l.size()-9];
            for (int i = 9; i < l.size(); i++) {
                taxa[i-9] = l.get(i);
            }
            ExecutorService pool = Executors.newFixedThreadPool(PGLConstraints.parallelLevel);
            List<Future<GenoSiteBlockVCF>> resultList = new ArrayList<>();
            int siteCount = 0;
            int startIndex = 0;
            List<String> lines = new ArrayList ();
            StringBuilder sb = new StringBuilder();
            while ((temp = br.readLine()) != null) {
                lines.add(temp);
                if (lines.size()%GenoSiteBlockVCF.blockSize == 0) {
                    GenoSiteBlockVCF gsb = new GenoSiteBlockVCF(lines, startIndex);
                    Future<GenoSiteBlockVCF> result = pool.submit(gsb);
                    resultList.add(result);
                    startIndex+=lines.size();
                    lines = new ArrayList<>();
                }
                siteCount++;
                if (siteCount%1000000 == 0) {
                    sb.setLength(0);
                    sb.append("Read in ").append(siteCount).append(" SNPs from ").append(infileS);
                    System.out.println(sb.toString());
                }
            }
            br.close();
            if (lines.size() != 0) {
                GenoSiteBlockVCF gsb = new GenoSiteBlockVCF(lines, startIndex);
                Future<GenoSiteBlockVCF> result = pool.submit(gsb);
                resultList.add(result);
            }
            pool.shutdown();
            pool.awaitTermination(Long.MAX_VALUE, TimeUnit.MICROSECONDS);
            this.genoSite = new BitSet[siteCount][];
            this.snps = new BiSNP[siteCount];
            for (int i = 0; i < resultList.size(); i++) {
                GenoSiteBlockVCF block = resultList.get(i).get();
                for (int j = 0; j < block.actBlockSize; j++) {
                    this.snps[block.getStartIndex()+j] = block.getSNPBlock()[j];
                    this.genoSite[block.getStartIndex()+j] = block.getGenoSiteBlock()[j];
                }
            }
            sb.setLength(0);
            sb.append("A total of ").append(this.getSiteNumber()).append(" SNPs are in ").append(infileS).append("\n");
            sb.append("Genotype table is successfully built");
            System.out.println(sb.toString());
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        this.transposeSiteToTaxon();
        mafs = new float[this.getSiteNumber()];
        Arrays.fill(mafs, Float.MIN_VALUE);
    }

    private void transposeSiteToTaxon () {
        long start = System.nanoTime();
        genoTaxon = new BitSet[this.getTaxaNumber()][3];
        Arrays.stream(genoTaxon).parallel().forEach(e -> {
            for (int i = 0; i < genoTaxon[0].length; i++) {
                e[i] = new BitSet(this.getSiteNumber());
            }
        });
        List<Integer> tIndexList = PArrayUtils.getIndexList(this.getTaxaNumber());
        tIndexList.parallelStream().forEach(i -> {
            for (int j = 0; j < genoTaxon[0].length; j++) {
                for (int k = 0; k < this.getSiteNumber(); k++) {
                    if (genoSite[k][j].get(i)) genoTaxon[i][j].set(k);
                }
            }
        });
        StringBuilder sb = new StringBuilder();
        sb.append("Transpose genoSite to genoTaxon takes ").append(Benchmark.getTimeSpanSeconds(start)).append(" seconds.");
        System.out.println(sb.toString());
    }

    public void transposeTaxonToSite () {
        long start = System.nanoTime();
        genoSite = new BitSet[this.getSiteNumber()][3];
        Arrays.stream(genoSite).parallel().forEach(e -> {
            for (int i = 0; i < genoSite[0].length; i++) {
                e[i] = new BitSet(this.getTaxaNumber());
            }
        });
        List<Integer> tIndexList = PArrayUtils.getIndexList(this.getTaxaNumber());
        tIndexList.parallelStream().forEach(i -> {
            for (int j = 0; j < genoTaxon[0].length; j++) {
                for (int k = 0; k < this.getSiteNumber(); k++) {
                    if (genoTaxon[i][j].get(k)) genoSite[k][j].set(i);
                }
            }
        });
        StringBuilder sb = new StringBuilder();
        sb.append("Transpose genoTaxon to genoSite takes ").append(Benchmark.getTimeSpanSeconds(start)).append(" seconds.");
        System.out.println(sb.toString());
    }

    @Override
    public void swap(int index1, int index2) {
        BiSNP tempB = null;
        tempB = snps[index1];
        snps[index1] = snps[index2];
        snps[index2] = tempB;
        float tempM;
        tempM = mafs[index1];
        mafs[index1] = mafs[index2];
        mafs[index2] = tempM;
        BitSet[] tempS = null;
        tempS = genoSite[index1];
        genoSite[index1] = genoSite[index2];
        genoSite[index2] = tempS;
    }

    @Override
    public int compare(int index1, int index2) {
        return snps[index1].compareTo(snps[index2]);
    }
}
