package pgl.infra.dna.genot;

import pgl.infra.dna.allele.AlleleEncoder;
import pgl.infra.dna.snp.BiSNP;
import pgl.infra.utils.Dyad;

import java.nio.ByteBuffer;
import java.util.BitSet;
import java.util.concurrent.Callable;

/**
 * Multithreading reading facilities of {@link GenotypeGrid}
 * <p>
 * Support binary format, see {@link GenoIOFormat}
 *
 * @author feilu
 */
public class GenoSiteBlockBinary implements Callable<GenoSiteBlockBinary> {
    public static final int blockSize = 4096;
    byte[][] lines = null;
    int startIndex = Integer.MIN_VALUE;
    int actBlockSize = Integer.MIN_VALUE;
    BitSet[][] genoSiteBlock = null;
    BiSNP[] snpBlock = null;

    public GenoSiteBlockBinary (byte[][] lines, int startIndex, int actBlockSize) {
        this.lines = lines;
        this.startIndex = startIndex;
        this.actBlockSize = actBlockSize;
    }

    public BiSNP[] getSNPBlock () {
        return this.snpBlock;
    }

    public BitSet[][] getGenoSiteBlock () {
        return this.genoSiteBlock;
    }

    public int getStartIndex () {
        return this.startIndex;
    }

    @Override
    public GenoSiteBlockBinary call() throws Exception {
        this.genoSiteBlock = new BitSet[this.actBlockSize][3];
        this.snpBlock = new BiSNP[this.actBlockSize];
        ByteBuffer bb = ByteBuffer.allocate(lines[0].length);
        for (int i = 0; i < this.actBlockSize; i++) {
            bb.put(lines[i]);
            Dyad<BiSNP, BitSet[]> d = buildFromBinaryLine(bb);
            snpBlock[i] = d.getFirstElement();
            genoSiteBlock[i] = d.getSecondElement();
        }
        lines = null;
        return this;
    }

    private Dyad<BiSNP, BitSet[]> buildFromBinaryLine (ByteBuffer bb) {
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
        BitSet[] genoSite = new BitSet[3];
        for (int i = 0; i < genoSite.length; i++) {
            byte[] ba = new byte[size];
            bb.get(ba);
            genoSite[i] = BitSet.valueOf(ba);
        }
        BiSNP snp = new BiSNP(chr, pos, refBase, altBase, null);
        snp.setReferenceAlleleFeature(refFeature);
        snp.setAlternativeAlleleFeature(altFeature);
        bb.clear();
        Dyad<BiSNP, BitSet[]> d = new Dyad<>(snp, genoSite);
        return d;
    }
}
