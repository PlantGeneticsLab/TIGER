package pgl.app.fastCall2;

import com.itextpdf.text.pdf.qrcode.ByteArray;
import com.mysql.cj.x.protobuf.MysqlxDatatypes;
import it.unimi.dsi.fastutil.bytes.ByteArrayList;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PArrayUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static pgl.infra.dna.allele.AlleleEncoder.alleleByteToBaseMap;

public class VariationLibrary {
    short chrom = Short.MIN_VALUE;
    int binStart = Integer.MIN_VALUE;
    int maoThresh = Integer.MIN_VALUE;
    int maxAltNum = 2;
    int[] potentialPositions = null;
    IntArrayList[] pAlleleLists = null;
    int[] positions = null;
    byte[][] codedAlleles = null;

    public VariationLibrary (List<IndividualGenotype> ingList, int maoThresh, int maxAltNum, short chrom, int binStart) {
        this.maxAltNum = maxAltNum;
        this.maoThresh = maoThresh;
        this.chrom = chrom;
        this.binStart = binStart;
        this.mergeIngs(ingList);
    }

    public void mergeIngs (List<IndividualGenotype> ingList) {
        IntOpenHashSet positionSet = new IntOpenHashSet();
        for (int i = 0; i < ingList.size(); i++) {
            int positionNumber = ingList.get(i).getPositionNumber();
            for (int j = 0; j < positionNumber; j++) {
                positionSet.add(ingList.get(i).getAllelePosition(j));
            }

        }

        potentialPositions = positionSet.toIntArray();
        Arrays.sort(potentialPositions);
        pAlleleLists = new IntArrayList[potentialPositions.length];
        for (int i = 0; i < potentialPositions.length; i++) {
            pAlleleLists[i] = new IntArrayList();
        }
        for (int i = 0; i < ingList.size(); i++) {
            int positionNum = ingList.get(i).getPositionNumber();
            int currentIndex = Integer.MIN_VALUE;
            for (int j = 0; j < positionNum; j++) {
                currentIndex = Arrays.binarySearch(potentialPositions, ingList.get(i).getAllelePosition(j));
                pAlleleLists[currentIndex].add(ingList.get(i).getCodedAlleleInfo(j));
            }
        }
        List<Integer> indexList = new ArrayList<>();
        int[][] alts = new int[potentialPositions.length][maxAltNum];
        int[][] altCounts = new int[potentialPositions.length][maxAltNum];
        for (int i = 0; i < potentialPositions.length; i++) {
            indexList.add(i);
        }
        indexList.parallelStream().forEach(i -> {
            IntOpenHashSet alleleSet = new IntOpenHashSet(pAlleleLists[i]);
            int[] alleles = alleleSet.toIntArray();
            Arrays.sort(alleles);
            int[] alleleCount = new int[alleles.length];
            int index = Integer.MIN_VALUE;
            for (int j = 0; j < pAlleleLists[i].size(); j++) {
                index = Arrays.binarySearch(alleles, pAlleleLists[i].getInt(j));
                alleleCount[index]++;
            }
            int[] countIndex = PArrayUtils.getIndicesByDescendingValue(alleleCount);
            int minNum = alleles.length;
            if (minNum > maxAltNum) minNum = maxAltNum;
            for (int j = 0; j < minNum; j++) {
                alts[i][j] = alleles[countIndex[j]];
                altCounts[i][j] = alleleCount[countIndex[j]];
            }
        });
        IntArrayList positionList = new IntArrayList();
        List<byte[]> codedAlleleLists = new ArrayList<>();
        ByteArrayList codedAlleleList = new ByteArrayList();
        for (int i = 0; i < alts.length; i++) {
            int varifiedNum = this.maxAltNum;
            for (int j = this.maxAltNum; j > 0; j--) {
                if (altCounts[i][j-1] < this.maoThresh) varifiedNum--;
            }
            if (varifiedNum == 0) continue;
            codedAlleleList.clear();
            positionList.add(FastCall2.getAllelePosition(alts[i][0], binStart));
            for (int j = 0; j < varifiedNum; j++) {
                codedAlleleList.add(FastCall2.getCodedAllele(alts[i][j]));
            }
            codedAlleleLists.add(codedAlleleList.toByteArray());
        }
        this.positions = positionList.toIntArray();
        codedAlleles = new byte[positions.length][];
        for (int i = 0; i < positions.length; i++) {
            this.codedAlleles[i] = codedAlleleLists.get(i);
        }
    }

    public void writeTextFileS (String outfileS) {
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Position\tAlts\tIndelLength");
            bw.newLine();
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < positions.length; i++) {
                sb.setLength(0);
                sb.append(positions[i]).append("\t");
                for (int j = 0; j < this.codedAlleles[i].length; j++) {
                    sb.append(alleleByteToBaseMap.get(FastCall2.getAlleleByteFromCodedAllele(this.codedAlleles[i][j]))).append(",");
                }
                sb.deleteCharAt(sb.length()-1).append("\t");
                for (int j = 0; j < this.codedAlleles[i].length; j++) {
                    sb.append(FastCall2.getIndelLengthFromCodedAllele(this.codedAlleles[i][j])).append(",");
                }
                sb.deleteCharAt(sb.length()-1);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.newLine();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
