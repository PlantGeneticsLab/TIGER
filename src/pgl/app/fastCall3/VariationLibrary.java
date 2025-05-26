package pgl.app.fastCall3;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PArrayUtils;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.util.*;

public class VariationLibrary implements Comparable<VariationLibrary> {
    short chrom = Short.MIN_VALUE;
    //set to -1, when bin-based vls are concatenated to a chrom-based vl
    int binStart = Integer.MIN_VALUE;
    int[] positions = null;
    AllelePackage[][] allelePacks = null;
    private int[] potentialPositions = null;
    private ArrayList<AllelePackage>[] potentialAllelePackLists = null;


    public VariationLibrary(short chrom, int binStart, int[] positions, AllelePackage[][] allelePacks) {
        this.chrom = chrom;
        this.binStart = binStart;
        this.positions = positions;
        this.allelePacks = allelePacks;
    }

    public static VariationLibrary getInstance(List<VariationLibrary> vList) {
        Collections.sort(vList);
        IntArrayList positionList = new IntArrayList();
        List<AllelePackage[]> alleleList = new ArrayList<>();
        for (int i = 0; i < vList.size(); i++) {
            for (int j = 0; j < vList.get(i).positions.length; j++) {
                positionList.add(vList.get(i).positions[j]);
                alleleList.add(vList.get(i).allelePacks[j]);
            }
        }
        int[] positions = positionList.toIntArray();
        AllelePackage[][] allelePacks = alleleList.toArray(new AllelePackage[alleleList.size()][]);
        VariationLibrary vl = new VariationLibrary(vList.get(0).chrom, -1, positions, allelePacks);
        return vl;
    }

    public VariationLibrary(String infileS) {
        this.readBinaryFileS(infileS);
    }

    public VariationLibrary(List<IndividualGenotype> ingList, int maoThresh, int maxAltNum, short chrom, int binStart) {
        this.chrom = chrom;
        this.binStart = binStart;
        this.mergeIngs(ingList, maoThresh, maxAltNum);
    }

    private void mergeIngs (List<IndividualGenotype> ingList, int maoThresh, int maxAltNum) {
        IntOpenHashSet positionSet = new IntOpenHashSet();
        for (int i = 0; i < ingList.size(); i++) {
            int positionNumber = ingList.get(i).getPositionNumber();
            for (int j = 0; j < positionNumber; j++) {
                positionSet.add(ingList.get(i).getAlleleChromPosition(j));
            }

        }
        potentialPositions = positionSet.toIntArray();
        Arrays.sort(potentialPositions);
        potentialAllelePackLists = new ArrayList[potentialPositions.length];
        for (int i = 0; i < potentialPositions.length; i++) {
            potentialAllelePackLists[i] = new ArrayList<>();
        }
        for (int i = 0; i < ingList.size(); i++) {
            int positionNum = ingList.get(i).getPositionNumber();
            int currentIndex = Integer.MIN_VALUE;
            for (int j = 0; j < positionNum; j++) {
                currentIndex = Arrays.binarySearch(potentialPositions, ingList.get(i).getAlleleChromPosition(j));
                potentialAllelePackLists[currentIndex].add(new AllelePackage(ingList.get(i).getAllelePack(j)));
            }
        }
        List<Integer> indexList = new ArrayList<>();
        AllelePackage[][] alts = new AllelePackage[potentialPositions.length][maxAltNum];
        int[][] altCounts = new int[potentialPositions.length][maxAltNum];
        for (int i = 0; i < potentialPositions.length; i++) {
            indexList.add(i);
        }
        indexList.parallelStream().forEach(i -> {
            Set<AllelePackage> alleleSet = new HashSet<>();
            for (int j = 0; j < potentialAllelePackLists[i].size(); j++) {
                alleleSet.add(potentialAllelePackLists[i].get(j));
            }
            AllelePackage[] alleles = alleleSet.toArray(new AllelePackage[alleleSet.size()]);
            Arrays.sort(alleles);
            int[] alleleCount = new int[alleles.length];
            int index = Integer.MIN_VALUE;
            for (int j = 0; j < potentialAllelePackLists[i].size(); j++) {
                index = Arrays.binarySearch(alleles, potentialAllelePackLists[i].get(j));
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
        List<AllelePackage[]> allelePackLists = new ArrayList<AllelePackage[]>();
        List<AllelePackage> allelePackList = new ArrayList<>();
        for (int i = 0; i < alts.length; i++) {
            int varifiedNum = maxAltNum;
            for (int j = maxAltNum; j > 0; j--) {
                if (altCounts[i][j-1] < maoThresh) varifiedNum--;
            }
            if (varifiedNum == 0) continue;
            allelePackList.clear();
            for (int j = 0; j < varifiedNum; j++) {
                allelePackList.add(alts[i][j]);
            }
            positionList.add(alts[i][0].getAlleleChromPosition(binStart));
            allelePackLists.add(allelePackList.toArray(new AllelePackage[allelePackList.size()]));
        }
        this.positions = positionList.toIntArray();
        this.allelePacks = allelePackLists.toArray(new AllelePackage[allelePackLists.size()][]);
    }

    public AllelePackage[] getAllelePacks (int positionIndex) {
        return this.allelePacks[positionIndex];
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
                dos.writeByte((byte) this.allelePacks[i].length);
                for (int j = 0; j < allelePacks[i].length; j++) {
                    int num = allelePacks[i][j].getAllelePackSize();
                    for (int k = 0; k < num; k++) {
                        dos.writeInt(allelePacks[i][j].getAllelePack()[k]);
                    }
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

    public void writeBinaryFileS (String outfileS, int[] customPositions) {
        IntArrayList indexList = new IntArrayList();
        for (int i = 0; i < positions.length; i++) {
            int index = Arrays.binarySearch(customPositions, positions[i]);
            if (index < 0) continue;
            indexList.add(i);
        }
        int[] indices = indexList.toIntArray();
        try {
            DataOutputStream dos = IOUtils.getBinaryGzipWriter(outfileS);
            dos.writeShort(chrom);
            dos.writeInt(binStart);
            dos.writeInt(indices.length);
            for (int i = 0; i < indices.length; i++) {
                dos.writeInt(positions[indices[i]]);
                dos.writeByte((byte) this.allelePacks[indices[i]].length);
                for (int j = 0; j < allelePacks[indices[i]].length; j++) {
                    int num = allelePacks[indices[i]][j].getAllelePackSize();
                    for (int k = 0; k < num; k++) {
                        dos.writeInt(allelePacks[indices[i]][j].getAllelePack()[k]);
                    }
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
            allelePacks = new AllelePackage[positionNum][];
            for (int i = 0; i < positionNum; i++) {
                positions[i] = dis.readInt();
                int alleleNum = dis.readByte();
                allelePacks[i] = new AllelePackage[alleleNum];
                for (int j = 0; j < alleleNum; j++) {
                    int firstInt = dis.readInt();
                    int packSize = AllelePackage.getAllelePackSizeFromFirstInt(firstInt);
                    int[] allelePack = new int[packSize];
                    allelePack[0] = firstInt;
                    for (int k = 1; k < packSize; k++) {
                        allelePack[k] = dis.readInt();
                    }
                    allelePacks[i][j] = new AllelePackage(allelePack);
                }
            }
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public void writeTextFileS (String outfileS) {
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Position\tAlts\tIndelLength\tIndelSeq");
            bw.newLine();
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < positions.length; i++) {
                sb.setLength(0);
                sb.append(positions[i]).append("\t");
                for (int j = 0; j < this.allelePacks[i].length; j++) {
                    sb.append(allelePacks[i][j].getAlleleBase()).append(",");
                }
                sb.deleteCharAt(sb.length()-1).append("\t");
                int totalLength = 0;
                for (int j = 0; j < this.allelePacks[i].length; j++) {
                    int indelLength = allelePacks[i][j].getIndelLength();
                    sb.append(indelLength).append(",");
                    totalLength+=indelLength;
                }
                sb.deleteCharAt(sb.length()-1);
                if (totalLength > 0) {
                    sb.append("\t");
                    for (int j = 0; j < this.allelePacks[i].length; j++) {
                        sb.append(allelePacks[i][j].getIndelSeq()).append(",");
                    }
                    sb.deleteCharAt(sb.length()-1);
                }
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
