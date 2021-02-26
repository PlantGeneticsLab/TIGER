package pgl.app.fastCall2;

import com.koloboke.collect.map.hash.HashByteByteMap;
import com.koloboke.collect.map.hash.HashByteByteMaps;
import pgl.AppUtils;
import pgl.infra.dna.allele.AlleleEncoder;
import pgl.infra.utils.*;

import java.util.*;

public class FastCall2 {
    //Current step ID of the pipeline
    int step = Integer.MIN_VALUE;
    //genome block for individual ing
    static int binSize = 50000;
    //A, C, G, T, -, +
    static final byte[] pileupAlleleAscIIs = {65, 67, 71, 84, 45, 43};

    static final HashByteByteMap pileupAscIIToAlleleByteMap =
            HashByteByteMaps.getDefaultFactory().withDefaultValue((byte)-1).newImmutableMap(pileupAlleleAscIIs, AlleleEncoder.alleleBytes);

    public FastCall2 (String parameterFileS) {
        this.runSteps(parameterFileS);
    }

    private void runSteps(String parameterFileS) {
        Dyad<List<String>, List<String>> d = AppUtils.getParameterList(parameterFileS);
        List<String> pLineList = d.getFirstElement();
        List<String> sLineList = d.getSecondElement();
        this.step = Integer.parseInt(sLineList.get(0).split("\\s+")[1]);
        if (step == 1) {
            new DiscoverVariation(pLineList);
        }
        else if (step == 2) {
            new BuildVariationLibrary(pLineList);
        }
    }

    static Dyad<int[][], int[]> getBins (int regionStart, int regionEnd) {
        int actualChrLength = regionEnd - regionStart;
        //starting from actual genome position
        int[][] binBound = PArrayUtils.getSubsetsIndicesBySubsetSize (actualChrLength, binSize);
        int[] binStarts = new int[binBound.length];
        for (int i = 0; i < binBound.length; i++) {
            binBound[i][0] = binBound[i][0]+regionStart;
            binBound[i][1] = binBound[i][1]+regionStart;
            binStarts[i] = binBound[i][0];
        }
        return new Dyad<>(binBound, binStarts);
    }

    static int getCodedPosAlleleIndelLength(int binStart, int position, byte alleleByte, int indelLength) {
        int v = (position-binStart) << 8;
        v = v + (alleleByte << 5);
        if (indelLength > 32) indelLength = 32;
        return (v + indelLength);
    }

    static int getAllelePosition(int codedPosAlleleIndelLength, int binStart) {
        int v = (codedPosAlleleIndelLength >> 8) + binStart;
        return v;
    }

    static byte getCodedAllele (int codedPosAlleleIndelLength) {
        return (byte)codedPosAlleleIndelLength;
    }

    static byte getAlleleByte (int codedPosAlleleIndelLength) {
        int v = (255 & codedPosAlleleIndelLength) >> 5;
        return (byte)v;
    }

    static char getAlleleBase (int codedPosAlleleIndelLength) {
        return AlleleEncoder.alleleByteToBaseMap.get(getAlleleByte(codedPosAlleleIndelLength));
    }

    static byte getIndelLength (int codedPosAlleleIndelLength) {
        return (byte)(31 & codedPosAlleleIndelLength);
    }

    static byte getAlleleByteFromCodedAllele (byte codedAllele) {
        int v = (codedAllele>>>5)&7;
        return (byte)v;
    }

    static char getAlleleBaseFromCodedAllele (byte codedAllele) {
        return AlleleEncoder.alleleByteToBaseMap.get(getAlleleByteFromCodedAllele(codedAllele));
    }

    static byte getIndelLengthFromCodedAllele (byte codedAllele) {
        return (byte)(31 & codedAllele);
    }

}
