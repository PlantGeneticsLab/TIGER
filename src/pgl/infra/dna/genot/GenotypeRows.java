package pgl.infra.dna.genot;

import pgl.PGLConstraints;
import pgl.infra.dna.allele.AlleleEncoder;
import pgl.infra.dna.allele.AlleleType;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PArrayUtils;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import pgl.infra.pos.ChrPos;

/**
 * The bit version implementation of a genotype table, see {@link GenotypeTable}
 * <p>
 * The class implements bitset genotype by site, which means it is smaller than {@link GenotypeGrid} in memory,
 * but slower in taxon based calculation, e.g. matrix of genetic divergence between taxa
 * <p>
 *  * Supports only bi-allelic SNPs, 3rd+ allele will be ignored. Allele depth is ignored.
 *
 * @author feilu
 */
public class GenotypeRows implements GenotypeTable {
    /**
     * The taxa in the genotype table
     */
    String[] taxa = null;
    /**
     * SNP, alleles and their presentation in individuals
     */
    SiteGenotype[] geno = null;

    /**
     * Construct an object
     */
    public GenotypeRows() {
        
    }

    /**
     * Construct an object by reading a file
     * @param infileS
     * @param format
     */
    public GenotypeRows(String infileS, GenoIOFormat format) {
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
    }

    /**
     * Construct an object
     * @param geno
     * @param taxa
     */
    public GenotypeRows (SiteGenotype[] geno, String[] taxa) {
        this.taxa = taxa;
        this.geno = geno;
    }
    
    @Override
    public int getTaxaNumber() {
        return taxa.length;
    }

    @Override
    public String getTaxonName(int taxonIndex) {
        return taxa[taxonIndex];
    }

    @Override
    public String[] getTaxaNames () {
//        String nTaxa[]  = new String[this.getTaxaNumber()];
//        System.arraycopy(taxa, 0, nTaxa, 0, nTaxa.length);
//        return nTaxa;
        return this.taxa;
    }

    @Override
    public int getTaxonIndex(String taxon) {
        return Arrays.binarySearch(taxa, taxon);
    }
    
    @Override
    public int getSiteNumber () {
        return geno.length;
    }
    
    @Override
    public short getChromosome(int siteIndex) {
        return geno[siteIndex].getChromosome();
    }

    @Override
    public int getPosition(int siteIndex) {
        return geno[siteIndex].getPosition();
    }

    @Override
    public void sortBySite() {
        Arrays.sort(geno);
    }

    @Override
    public void sortByTaxa() {
        int[] indices = PArrayUtils.getIndicesByAscendingValue(this.taxa);
        Arrays.sort(this.taxa);
        List<SiteGenotype> genoList = Arrays.asList(this.geno);
        genoList.parallelStream().forEach(g -> {
            g.sortByTaxa(indices);
        });
        geno = genoList.toArray(new SiteGenotype[genoList.size()]);
    }

    @Override
    public byte getGenotypeByte (int siteIndex, int taxonIndex) {
        return geno[siteIndex].getGenotypeByte(taxonIndex);
    }

    @Override
    public boolean isHeterozygous(int siteIndex, int taxonIndex) {
        return geno[siteIndex].isHeterozygous(taxonIndex);
    }

    @Override
    public boolean isHomozygous(int siteIndex, int taxonIndex) {
        return geno[siteIndex].isHomozygous(taxonIndex);
    }

    @Override
    public boolean isMissing(int siteIndex, int taxonIndex) {
        return geno[siteIndex].isMissing(taxonIndex);
    }

    @Override
    public boolean isPhase1Alternative(int siteIndex, int taxonIndex) {
        return geno[siteIndex].isPhase1Alternative(taxonIndex);
    }

    @Override
    public boolean isPhase2Alternative(int siteIndex, int taxonIndex) {
        return geno[siteIndex].isPhase2Alternative(taxonIndex);
    }

    @Override
    public boolean isPhase1Reference(int siteIndex, int taxonIndex) {
        return geno[siteIndex].isPhase1Reference(taxonIndex);
    }

    @Override
    public boolean isPhase2Reference(int siteIndex, int taxonIndex) {
        return geno[siteIndex].isPhase2Reference(taxonIndex);
    }

    @Override
    public boolean isAlternativeAlleleTypeOf(AlleleType at, int siteIndex) {
        return geno[siteIndex].isAlternativeAlleleTypeOf(at);
    }

    @Override
    public boolean isReferenceAlleleTypeOf(AlleleType at, int siteIndex) {
        return geno[siteIndex].isReferenceAlleleTypeOf(at);
    }

    @Override
    public int getSiteIndex(short chromosome, int position) {
        ChrPos query = new ChrPos (chromosome, position);
        int index = Arrays.binarySearch(geno, query);
        return index;
    }

    @Override
    public int getMissingNumberBySite(int siteIndex) {
        return geno[siteIndex].getMissingNumber();
    }

    @Override
    public int getMissingNumberByTaxon(int taxonIndex) {
        int cnt = 0;
        for (int i = 0; i < this.getSiteNumber(); i++) {
            if (geno[i].isMissing(taxonIndex)) cnt++;
        }
        return cnt;
    }
    
    @Override
    public int getNonMissingNumberBySite(int siteIndex) {
        return this.getTaxaNumber()-this.getMissingNumberBySite(siteIndex);
    }

    @Override
    public int getNonMissingNumberByTaxon(int taxonIndex) {
        return this.getSiteNumber()-this.getMissingNumberByTaxon(taxonIndex);
    }

    @Override
    public int getHomozygoteNumberBySite(int siteIndex) {
        return geno[siteIndex].getHomozygoteNumber();
    }

    @Override
    public int getHomozygoteNumberByTaxon(int taxonIndex) {
        int cnt = 0;
        for (int i = 0; i < this.getSiteNumber(); i++) {
            if (geno[i].isHomozygous(taxonIndex)) cnt++;
        }
        return cnt;
    }

    @Override
    public int getHeterozygoteNumberBySite(int siteIndex) {
        return geno[siteIndex].getHeterozygoteNumber();
    }

    @Override
    public int getHeterozygoteNumberByTaxon(int taxonIndex) {
        int cnt = 0;
        for (int i = 0; i < this.getSiteNumber(); i++) {
            if (geno[i].isHeterozygous(taxonIndex)) cnt++;
        }
        return cnt;
    }

    @Override
    public int getAlternativeAlleleNumberBySite(int siteIndex) {
        return geno[siteIndex].getAlternativeAlleleNumber();
    }

    @Override
    public int getAlternativeAlleleOccurrenceBySite(int siteIndex) {
        BitSet bs = geno[siteIndex].phase1;
        bs.or(geno[siteIndex].phase2);
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
        return geno[siteIndex].getMinorAlleleByte();
    }

    @Override
    public char getMinorAlleleBase(int siteIndex) {
        return AlleleEncoder.getAlleleBaseFromCoding(this.getMinorAlleleByte(siteIndex));
    }

    @Override
    public float getMinorAlleleFrequency(int siteIndex) {
        return geno[siteIndex].getMinorAlleleFrequency();
    }
    
    @Override
    public byte getMajorAlleleByte(int siteIndex) {
        return geno[siteIndex].getMinorAlleleByte();
    }

    @Override
    public char getMajorAlleleBase(int siteIndex) {
        return AlleleEncoder.getAlleleBaseFromCoding(this.getMinorAlleleByte(siteIndex));
    }

    @Override
    public float getMajorAlleleFrequency(int siteIndex) {
        return geno[siteIndex].getMajorAlleleFrequency();
    }

    @Override
    public byte getReferenceAlleleByte(int siteIndex) {
        return geno[siteIndex].getReferenceAlleleByte();
    }

    @Override
    public char getReferenceAlleleBase(int siteIndex) {
        return AlleleEncoder.getAlleleBaseFromCoding(this.getReferenceAlleleByte(siteIndex));
    }

    @Override
    public float getReferenceAlleleFrequency(int siteIndex) {
        return geno[siteIndex].getReferenceAlleleFrequency();
    }
    
    @Override
    public byte getAlternativeAlleleByte(int siteIndex) {
        return geno[siteIndex].getAlternativeAlleleByte();
    }

    @Override
    public char getAlternativeAlleleBase(int siteIndex) {
        return AlleleEncoder.getAlleleBaseFromCoding(this.getAlternativeAlleleByte(siteIndex));
    }

    @Override
    public float getAlternativeAlleleFrequency(int siteIndex) {
        return geno[siteIndex].getReferenceAlleleFrequency();
    }

    @Override
    public String getUnphasedVCFRecord(int siteIndex) {
        return geno[siteIndex].getUnphasedVCFOutput();
    }

    @Override
    public ByteBuffer getBinaryOutput (int siteIndex, ByteBuffer bb) {
        return geno[siteIndex].getBinaryOutput(bb);
    }

    @Override
    public void setAlternativeAlleleType(AlleleType at, int siteIndex) {
        geno[siteIndex].setAlternativeAlleleType(at);
    }

    @Override
    public void setReferenceAlleleType(AlleleType at, int siteIndex) {
        geno[siteIndex].setReferenceAlleleType(at);
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
        int size = endSiteIndex - startSiteIndex;
        int cnt = 0;
        double dxy = 0;
        for (int i = 0; i < size; i++) {
            float dxySite = this.getDxySite(taxonIndex1, taxonIndex2, startSiteIndex+i);
            if (Float.isNaN(dxySite)) continue;
            dxy+=dxySite;
            cnt++;
        }
        return (float)(dxy/cnt);
    }

    @Override
    public float getIBSDistance(int taxonIndex1, int taxonIndex2, int[] siteIndices) {
        int cnt = 0;
        double dxy = 0;
        for (int i = 0; i < siteIndices.length; i++) {
            float dxySite = this.getDxySite(taxonIndex1,taxonIndex2, siteIndices[i]);
            if (Float.isNaN(dxySite)) continue;
            dxy+=dxySite;
            cnt++;
        }
        return (float)(dxy/cnt);
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

    public float[][] getDxyMatrixFast10K() {
        int size = 10_000;
        if (this.getSiteNumber() > size) {
            int[] siteIndices = new int[size];
            double add = this.getSiteNumber()/size;
            for (int i = 0; i < size; i++) {
                siteIndices[i] = (int)(i*add);
            }
            return this.getIBSDistanceMatrix(siteIndices);
        }
        else {
            return this.getIBSDistanceMatrix();
        }
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

    /**
     * Return the dxy bwteen two taxa at a specific site
     * @param taxonIndex1
     * @param taxonIndex2
     * @param siteIndex
     * @return Float.NaN if no shared non-missing site available
     */
    private float getDxySite (int taxonIndex1, int taxonIndex2, int siteIndex) {
        if (this.isMissing(taxonIndex1, siteIndex)) return Float.NaN;
        if (this.isMissing(taxonIndex2, siteIndex)) return Float.NaN;
        int cnt1 = 0;
        int cnt2 = 0;
        if (this.geno[siteIndex].isPhase1Alternative(taxonIndex1)) cnt1++;
        if (this.geno[siteIndex].isPhase2Alternative(taxonIndex1)) cnt1++;
        if (this.geno[siteIndex].isPhase1Alternative(taxonIndex2)) cnt2++;
        if (this.geno[siteIndex].isPhase2Alternative(taxonIndex2)) cnt2++;
        if (cnt1 > cnt2) return (float)((double)(cnt1-cnt2)/2);
        return (float)((double)(cnt2-cnt1)/2);
    }

    @Override
    public GenotypeTable getSubGenotypeTableBySite(int[] siteIndices) {
        SiteGenotype[] ge = new SiteGenotype[siteIndices.length];
        for (int i = 0; i < siteIndices.length; i++) {
            ge[i] = geno[siteIndices[i]];
        }
        GenotypeRows gb = new GenotypeRows(ge, taxa);
        return gb;
    }

    @Override
    public GenotypeTable getSubGenotypeTableByTaxa(int[] taxaIndices) {
        String[] nTaxa = new String[taxaIndices.length];
        for (int i = 0; i < nTaxa.length; i++) {
            nTaxa[i] = taxa[taxaIndices[i]];
        }
        SiteGenotype[] nGeno = new SiteGenotype[this.getSiteNumber()];
        List<Integer> indexList = new ArrayList<>(this.getSiteNumber());
        for (int i = 0; i < this.getSiteNumber(); i++) {
            indexList.add(i);
        }
        indexList.parallelStream().forEach(index -> {
            nGeno[index] = geno[index].getSubGenotypeByTaxa(taxaIndices);
        });
        return new GenotypeRows(nGeno, nTaxa);
    }

    /**
     * Build an object from a binary genotype file
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
            List<Future<SGBBlockBinary>> resultList = new ArrayList<>();
            byte[][] input = new byte[SGBBlockBinary.blockSize][GenotypeExport.getByteSizeOfSiteInBinary(taxaNumber)];
            int siteCount = 0;
            int startIndex = 0;
            int actBlockSize = 0;
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < siteNumber; i++) {
                dis.read(input[actBlockSize]);
                siteCount++;
                actBlockSize++;
                if (siteCount%SGBBlockBinary.blockSize == 0) {
                    SGBBlockBinary sgb = new SGBBlockBinary(input, startIndex, actBlockSize, taxaNumber);
                    Future<SGBBlockBinary> result = pool.submit(sgb);
                    resultList.add(result);
                    startIndex+=actBlockSize;
                    actBlockSize = 0;
                    input = new byte[SGBBlockBinary.blockSize][GenotypeExport.getByteSizeOfSiteInBinary(taxaNumber)];
                }
                if (siteCount%1000000 == 0) {
                    sb.setLength(0);
                    sb.append("Read in ").append(siteCount).append(" SNPs from ").append(infileS);
                    System.out.println(sb.toString());
                }
            }
            dis.close();
            if (actBlockSize != 0) {
                SGBBlockBinary sgb = new SGBBlockBinary(input, startIndex, actBlockSize, taxaNumber);
                Future<SGBBlockBinary> result = pool.submit(sgb);
                resultList.add(result);
            }
            pool.shutdown();
            pool.awaitTermination(Long.MAX_VALUE, TimeUnit.MICROSECONDS);
            this.geno = new SiteGenotype[siteCount];
            for (int i = 0; i < resultList.size(); i++) {
                SGBBlockBinary block = resultList.get(i).get();
                for (int j = 0; j < block.actBlockSize; j++) {
                    geno[block.getStartIndex()+j] = block.getSiteGenotypes()[j];
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
    }

    /**
     * Class for parallel reading in a binary genotype file
     */
    class SGBBlockBinary implements Callable<SGBBlockBinary> {
        public static final int blockSize = 4096;
        byte[][] lines = null;
        int startIndex = Integer.MIN_VALUE;
        int actBlockSize = Integer.MIN_VALUE;
        short taxaNumber = Short.MIN_VALUE;
        SiteGenotype[] sgbArray = null;

        public SGBBlockBinary (byte[][] lines, int startIndex, int actBlockSize, short taxaNumber) {
            this.lines = lines;
            this.startIndex = startIndex;
            this.actBlockSize = actBlockSize;
            this.taxaNumber = taxaNumber;
        }

        public int getStartIndex () {
            return this.startIndex;
        }

        public SiteGenotype[] getSiteGenotypes () {
            return this.sgbArray;
        }

        @Override
        public SGBBlockBinary call() throws Exception {
            this.sgbArray = new SiteGenotype[this.actBlockSize];
            ByteBuffer bb = ByteBuffer.allocate(lines[0].length);
            for (int i = 0; i < this.actBlockSize; i++) {
                bb.put(lines[i]);
                sgbArray[i] = SiteGenotype.buildFromBinaryLine(bb, taxaNumber);
            }
            lines = null;
            return this;
        }
    }

    /**
     * Build an object from a VCF file
     * @param infileS
     */
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
            List<Future<SGBBlockVCF>> resultList = new ArrayList<>();
            int siteCount = 0;
            int startIndex = 0;
            List<String> lines = new ArrayList ();
            StringBuilder sb = new StringBuilder();
            while ((temp = br.readLine()) != null) {
                lines.add(temp);
                if (lines.size()%SGBBlockVCF.blockSize == 0) {
                    SGBBlockVCF sgb = new SGBBlockVCF(lines, startIndex);
                    Future<SGBBlockVCF> result = pool.submit(sgb);
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
                SGBBlockVCF sgb = new SGBBlockVCF(lines, startIndex);
                Future<SGBBlockVCF> result = pool.submit(sgb);
                resultList.add(result);
            }
            pool.shutdown();
            pool.awaitTermination(Long.MAX_VALUE, TimeUnit.MICROSECONDS);
            this.geno = new SiteGenotype[siteCount];
            for (int i = 0; i < resultList.size(); i++) {
                SGBBlockVCF block = resultList.get(i).get();
                for (int j = 0; j < block.actBlockSize; j++) {
                    geno[block.getStartIndex()+j] = block.getSiteGenotypes()[j];
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
    }

    /**
     * Class for parallel reading in VCF
     */
    class SGBBlockVCF implements Callable<SGBBlockVCF> {
        public static final int blockSize = 4096;
        List<String> lines = null;
        int startIndex = Integer.MIN_VALUE;
        int actBlockSize = Integer.MIN_VALUE;
        SiteGenotype[] sgbArray = null;

        public SGBBlockVCF (List<String> lines, int startIndex) {
            this.lines = lines;
            this.startIndex = startIndex;
            this.actBlockSize = lines.size();
        }

        public int getStartIndex () {
            return this.startIndex;
        }

        public SiteGenotype[] getSiteGenotypes () {
            return this.sgbArray;
        }

        @Override
        public SGBBlockVCF call() throws Exception {
            this.sgbArray = new SiteGenotype[this.actBlockSize];
            for (int i = 0; i < this.actBlockSize; i++) {
                sgbArray[i] = SiteGenotype.buildFromVCFLine(lines.get(i));
            }
            lines = null;
            return this;
        }
    }
}
