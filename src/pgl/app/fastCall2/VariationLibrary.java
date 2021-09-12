package pgl.app.fastCall2;

import htsjdk.samtools.util.IOUtil;
import it.unimi.dsi.fastutil.bytes.ByteArrayList;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.shorts.ShortArrayList;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PArrayUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class VariationLibrary implements Comparable<VariationLibrary> {
    short chrom = Short.MIN_VALUE;
    //set to -1, when bin-based vls are concatenated to a chrom-based vl
    int binStart = Integer.MIN_VALUE;
    int[] positions = null;
    short[][] codedAlleles = null;

    private int[] potentialPositions = null;
    private IntArrayList[] pAlleleLists = null;


    public VariationLibrary (short chrom, int binStart, int[] positions, short[][] codedAlleles) {
        this.chrom = chrom;
        this.binStart = binStart;
        this.positions = positions;
        this.codedAlleles = codedAlleles;
    }

    public static VariationLibrary getInstance(List<VariationLibrary> vList) {
        Collections.sort(vList);
        IntArrayList positionList = new IntArrayList();
        List<short[]> alleleList = new ArrayList<>();
        for (int i = 0; i < vList.size(); i++) {
            for (int j = 0; j < vList.get(i).positions.length; j++) {
                positionList.add(vList.get(i).positions[j]);
                alleleList.add(vList.get(i).codedAlleles[j]);
            }
        }
        int[] positions = positionList.toIntArray();
        short[][] codedAlleles = alleleList.toArray(new short[alleleList.size()][]);
        VariationLibrary vl = new VariationLibrary(vList.get(0).chrom, -1, positions, codedAlleles);
        return vl;
    }

    public VariationLibrary (String infileS) {
        this.readBinaryFileS(infileS);
    }

    public VariationLibrary (List<IndividualGenotype> ingList, int maoThresh, int maxAltNum, short chrom, int binStart) {
        this.chrom = chrom;
        this.binStart = binStart;
        this.mergeIngs(ingList, maoThresh, maxAltNum);
    }

    private void mergeIngs (List<IndividualGenotype> ingList, int maoThresh, int maxAltNum) {
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
        List<short[]> codedAlleleLists = new ArrayList<>();
        ShortArrayList codedAlleleList = new ShortArrayList();
        for (int i = 0; i < alts.length; i++) {
            int varifiedNum = maxAltNum;
            for (int j = maxAltNum; j > 0; j--) {
                if (altCounts[i][j-1] < maoThresh) varifiedNum--;
            }
            if (varifiedNum == 0) continue;
            codedAlleleList.clear();
            positionList.add(FastCall2.getAllelePosition(alts[i][0], binStart));
            for (int j = 0; j < varifiedNum; j++) {
                codedAlleleList.add(FastCall2.getCodedAllele(alts[i][j]));
            }
            codedAlleleLists.add(codedAlleleList.toShortArray());
        }
        this.positions = positionList.toIntArray();
        codedAlleles = new short[positions.length][];
        for (int i = 0; i < positions.length; i++) {
            this.codedAlleles[i] = codedAlleleLists.get(i);
        }
    }

    public String[] getAlts (int positionIndex, List<String> l) {
        l.clear();
        for (int i = 0; i < this.codedAlleles[positionIndex].length; i++) {
            l.add(String.valueOf(FastCall2.getAlleleBaseFromCodedAllele(this.codedAlleles[positionIndex][i])));
        }
        String[] result = l.toArray(new String[l.size()]);
        return result;
    }

    public int getPositionIndex (int pos) {
        return Arrays.binarySearch(positions, pos);
    }

    /**
     * Return the index of the next position in the library, inclusive
     * @param pos
     * @return Integer.MIN_VALUE if the query pos is greater than the end of the position list
     */
    public int getStartIndex (int pos) {
        int index = this.getPositionIndex(pos);
        if (index < 0) {
            index = -index -1;
            if (index == positions.length) {
                return Integer.MIN_VALUE;
            }
        }
        return index;
    }

    /**
     * Return the index of the previous position in the library, exclusive
     * @param pos
     * @return Integer.MIN_VALUE if the query pos is less than the start of the position list
     */
    public int getEndIndex (int pos) {
        int index = this.getPositionIndex(pos);
        if (index < 0) {
            index = -index -2;
            if (index < 0) return Integer.MIN_VALUE;
        }
        return index+1;
    }

    public short getChrom () {
        return this.chrom;
    }

    public int getPosition (int positionIndex) {
        return this.positions[positionIndex];
    }

    public void writeBinaryFileS (String outfileS) {
        try {
            DataOutputStream dos = IOUtils.getBinaryGzipWriter(outfileS);
            dos.writeShort(chrom);
            dos.writeInt(binStart);
            dos.writeInt(positions.length);
            for (int i = 0; i < positions.length; i++) {
                dos.writeInt(positions[i]);
                dos.writeByte((byte)codedAlleles[i].length);
                for (int j = 0; j < codedAlleles[i].length; j++) {
                    dos.writeShort(codedAlleles[i][j]);
                }
            }
            dos.flush();
            dos.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        StringBuilder sb = new StringBuilder();
        sb.append(positions.length).append(" polymorphic sites are written to ").append(outfileS);
        System.out.println(sb.toString());
    }

    private void readBinaryFileS (String infileS) {
        try {
            DataInputStream dis = IOUtils.getBinaryGzipReader(infileS);
            chrom = dis.readShort();
            binStart = dis.readInt();
            int positionNum = dis.readInt();
            positions = new int[positionNum];
            codedAlleles = new short[positionNum][];
            for (int i = 0; i < positionNum; i++) {
                positions[i] = dis.readInt();
                int alleleNum = dis.readByte();
                codedAlleles[i] = new short[alleleNum];
                for (int j = 0; j < alleleNum; j++) {
                    codedAlleles[i][j] = dis.readShort();
                }
            }
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public void writeTextFileS (String outfileS, int[] positionIndices) {
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Position\tAlts\tIndelLength");
            bw.newLine();
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < positionIndices.length; i++) {
                sb.setLength(0);
                sb.append(positions[positionIndices[i]]).append("\t");
                for (int j = 0; j < this.codedAlleles[positionIndices[i]].length; j++) {
                    sb.append((FastCall2.getAlleleBaseFromCodedAllele(this.codedAlleles[positionIndices[i]][j]))).append(",");
                }
                sb.deleteCharAt(sb.length()-1).append("\t");
                for (int j = 0; j < this.codedAlleles[positionIndices[i]].length; j++) {
                    sb.append(FastCall2.getIndelLengthFromCodedAllele(this.codedAlleles[positionIndices[i]][j])).append(",");
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
            System.exit(1);
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
                    sb.append((FastCall2.getAlleleBaseFromCodedAllele(this.codedAlleles[i][j]))).append(",");
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
            System.exit(1);
        }
    }

    @Override
    public int compareTo(VariationLibrary o) {
        if (this.binStart < o.binStart) return -1;
        else if (this.binStart > o.binStart) return 1;
        return 0;
    }
}
