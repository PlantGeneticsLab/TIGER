package pgl.infra.popg;

import gnu.trove.list.array.TByteArrayList;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.dna.allele.AlleleEncoder;
import pgl.infra.dna.allele.AlleleType;
import pgl.infra.dna.genot.GenotypeGrid;
import pgl.infra.dna.genot.GenotypeOperation;

import java.util.Arrays;

public class ChromosomeFd {

    //The input genotype table containing genotypes of P1, P2, and P3
    GenotypeGrid oriGT = null;

    //The genotype table of P1, only SNPs with ancestral allele were picked up from the original genotype table
    GenotypeGrid p1 = null;

    //The genotype table of P2, only SNPs with ancestral allele were picked up from the original genotype table
    GenotypeGrid p2 = null;

    //The genotype table of P3, only SNPs with ancestral allele were picked up from the original genotype table
    GenotypeGrid p3 = null;

    double regionD = Double.NaN;

    double regionFd = Double.NaN;

    int windowSize = Integer.MIN_VALUE;

    int stepSize = Integer.MIN_VALUE;

    int snpWindowSize = Integer.MIN_VALUE;

    int snpStepSize = Integer.MIN_VALUE;

    int minAncestralSNPNumberInWindow = 5;

    public ChromosomeFd(GenotypeGrid oriGT, String[] p1Taxa, String[] p2Taxa, String[] p3Taxa, int[] ancestralPos, char[] oriAncestralAlleles) {
        this.oriGT = oriGT;
        this.initialize(p1Taxa, p2Taxa, p3Taxa, ancestralPos, oriAncestralAlleles);
    }


    private void initialize (String[] p1Taxa, String[] p2Taxa, String[] p3Taxa, int[] ancestralPos, char[] oriAncestralAlleles) {
        this.oriGT = oriGT;
        TIntArrayList posIndexList = new TIntArrayList();
        TByteArrayList ancestralAlleleByteList = new TByteArrayList();
        int siteNumber = oriGT.getSiteNumber();
        int currentPos = Integer.MIN_VALUE;
        int index = Integer.MIN_VALUE;
        byte currentAlleleByte = Byte.MIN_VALUE;
        for (int i = 0; i < siteNumber; i++) {
            currentPos = oriGT.getPosition(i);
            index = Arrays.binarySearch(ancestralPos, currentPos);
            if (index < 0) continue;
            currentAlleleByte = AlleleEncoder.getAlleleCodingFromBase(oriAncestralAlleles[index]);
            if (currentAlleleByte == oriGT.getReferenceAlleleByte(i) || currentAlleleByte == oriGT.getAlternativeAlleleByte(i)) {
                posIndexList.add(i);
                ancestralAlleleByteList.add(currentAlleleByte);
            }
        }
        GenotypeGrid gt = GenotypeOperation.getSubsetGenotypeBySite(oriGT, posIndexList.toArray());
        byte[] ancestralAlleleBytes = ancestralAlleleByteList.toArray();
        p1 = this.initializeGTPop(gt, p1Taxa, ancestralAlleleBytes);
        p2 = this.initializeGTPop(gt, p2Taxa, ancestralAlleleBytes);
        p3 = this.initializeGTPop(gt, p3Taxa, ancestralAlleleBytes);
    }

    private GenotypeGrid initializeGTPop (GenotypeGrid gt, String[] subTaxa, byte[] ancestralAlleleBytes) {
        Arrays.sort(subTaxa);
        int[] taxaIndex = new int[subTaxa.length];
        for (int i = 0; i < subTaxa.length; i++) {
            taxaIndex[i] = this.oriGT.getTaxonIndex(subTaxa[i]);
        }
        GenotypeGrid gg = GenotypeOperation.getSubsetGenotypeByTaxon(gt, taxaIndex);
        for (int i = 0; i < gg.getSiteNumber(); i++) {
            if (gg.getReferenceAlleleByte(i) == ancestralAlleleBytes[i]) {
                gg.setReferenceAlleleType(AlleleType.Ancestral, i);
            }
            else {
                gg.setAlternativeAlleleType(AlleleType.Ancestral, i);
            }
        }
        return gg;
    }

}
