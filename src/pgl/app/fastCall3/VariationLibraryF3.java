package pgl.app.fastCall3;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PArrayUtils;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.util.*;

/**
 * A class representing a variation library for the FastCall3 pipeline, storing genetic variant information
 * including positions and allele packs across genomic regions.
 *
 * <p>This class implements {@link Comparable} to enable sorting of variation libraries by their
 * genomic bin start positions. It provides functionality to read and write variation data in both
 * binary and text formats, and supports merging individual genotype data into a consolidated library.</p>
 *
 * <p>Key features include:
 * <ul>
 *   <li>Storage of variant positions and associated allele information</li>
 *   <li>Support for bin-based organization of genomic regions</li>
 *   <li>Efficient binary I/O operations for large datasets</li>
 *   <li>Methods for querying and filtering variants</li>
 *   <li>Support for merging multiple individual genotypes</li>
 * </ul>
 *
 * <p>The binary file format consists of:
 * <ol>
 *   <li>Chromosome number (short)</li>
 *   <li>Bin start position (int)</li>
 *   <li>Number of positions (int)</li>
 *   <li>For each position:
 *     <ul>
 *       <li>Position (int)</li>
 *       <li>Number of alleles (byte)</li>
 *       <li>Allele codes (int[])</li>
 *     </ul>
 *   </li>
 * </ol>
 *
 * @author Fei Lu
 * @version 3.0
 * @since 1.0
 * @see Comparable
 * @see IndividualGenotypeF3
 * @see ScanGenotypeF3
 */
public class VariationLibraryF3 implements Comparable<VariationLibraryF3> {
    short chrom = Short.MIN_VALUE;
    //set to -1, when bin-based vls are concatenated to a chrom-based vl
    int binStart = Integer.MIN_VALUE;
    int[] positions = null;
    int[][] allelePacks = null;
    private int[] potentialPositions = null;
    private IntArrayList[] potentialAllelePackLists = null;
    


    public VariationLibraryF3(short chrom, int binStart, int[] positions, int[][] allelePacks) {
        this.chrom = chrom;
        this.binStart = binStart;
        this.positions = positions;
        this.allelePacks = allelePacks;
    }

    public static VariationLibraryF3 getInstance(List<VariationLibraryF3> vList) {
        Collections.sort(vList);
        IntArrayList positionList = new IntArrayList();
        List<int[]> alleleList = new ArrayList<>();
        for (int i = 0; i < vList.size(); i++) {
            for (int j = 0; j < vList.get(i).getPositionNumber(); j++) {
                positionList.add(vList.get(i).getPosition(j));
                alleleList.add(vList.get(i).allelePacks[j]);
            }
        }
        int[] positions = positionList.toIntArray();
        int[][] allelePacks = alleleList.toArray(new int[alleleList.size()][]);
        VariationLibraryF3 vl = new VariationLibraryF3(vList.get(0).chrom, -1, positions, allelePacks);
        return vl;
    }

    public VariationLibraryF3(String infileS) {
        this.readBinaryFileS(infileS);
    }

    public VariationLibraryF3(List<IndividualGenotypeF3> ingList, int maoThresh, int maxAltNum, short chrom, int binStart) {
        this.chrom = chrom;
        this.binStart = binStart;
        this.mergeIngs(ingList, maoThresh, maxAltNum);
    }

    private void mergeIngs (List<IndividualGenotypeF3> ingList, int maoThresh, int maxAltNum) {
        IntOpenHashSet positionSet = new IntOpenHashSet();
        for (int i = 0; i < ingList.size(); i++) {
            int positionNumber = ingList.get(i).getPositionNumber();
            for (int j = 0; j < positionNumber; j++) {
                positionSet.add(ingList.get(i).getAlleleChromPosition(j));
            }
        }
        potentialPositions = positionSet.toIntArray();
        Arrays.sort(potentialPositions);
        potentialAllelePackLists = new IntArrayList[potentialPositions.length];
        for (int i = 0; i < potentialPositions.length; i++) {
            potentialAllelePackLists[i] = new IntArrayList();
        }
        for (int i = 0; i < ingList.size(); i++) {
            int positionNum = ingList.get(i).getPositionNumber();
            int currentIndex = Integer.MIN_VALUE;
            for (int j = 0; j < positionNum; j++) {
                currentIndex = Arrays.binarySearch(potentialPositions, ingList.get(i).getAlleleChromPosition(j));
                potentialAllelePackLists[currentIndex].add(ingList.get(i).getAllelePack(j));
            }
        }
        List<Integer> indexList = new ArrayList<>();
        int[][] alts = new int[potentialPositions.length][maxAltNum];
        int[][] altCounts = new int[potentialPositions.length][maxAltNum];
        for (int i = 0; i < potentialPositions.length; i++) {
            indexList.add(i);
        }
        indexList.parallelStream().forEach(i -> {
            IntOpenHashSet alleleSet = new IntOpenHashSet();
            for (int j = 0; j < potentialAllelePackLists[i].size(); j++) {
                alleleSet.add(potentialAllelePackLists[i].getInt(j));
            }
            int[] alleles = alleleSet.toIntArray();
            Arrays.sort(alleles);
            int[] alleleCount = new int[alleles.length];
            int index = Integer.MIN_VALUE;
            for (int j = 0; j < potentialAllelePackLists[i].size(); j++) {
                index = Arrays.binarySearch(alleles, potentialAllelePackLists[i].getInt(j));
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
        List<int[]> allelePackLists = new ArrayList<>();
        IntArrayList allelePackList = new IntArrayList();
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
            positionList.add(AllelePackageF3.getAlleleChromPosition(alts[i][0], binStart));
            allelePackLists.add(allelePackList.toIntArray());
        }
        this.positions = positionList.toIntArray();
        this.allelePacks = new int[positions.length][];
        for (int i = 0; i < allelePacks.length; i++) {
            allelePacks[i] = allelePackLists.get(i);
        }
    }

    public int[] getAllelePacks (int positionIndex) {
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

    /**
     * @return the chromosome number of this variation library
     */
    public short getChrom () {
        return this.chrom;
    }

    /**
     * @return the number of positions in this variation library
     */
    public int getPositionNumber () {
        return this.positions.length;
    }

    /**
     * @param positionIndex the index of the position to retrieve
     * @return the position at the given index
     * @since 1.0
     */
    public int getPosition (int positionIndex) {
        return this.positions[positionIndex];
    }

    /**
     * Writes the variation library in binary format to a file.
     *
     * The format is as follows:
     * 1. Short chrom
     * 2. Int binStart
     * 3. Int positionNumber
     * 4. For each position,
     *      Int position
     *      Byte alleleNumber
     *      For each allele, Int alleleCode
     *
     * @param outfileS the output file name
     */
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
                    dos.writeInt(this.allelePacks[i][j]);
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

    /**
     * Writes the variation library to a binary file in gzipped format, with specific positions.
     * The output file will have the same format as the standard binary file written by `writeBinaryFileS(String)`, but only sites present in the given array of positions will be written.
     * The order of the sites in the output file is the same as the order of the positions in the given array.
     * If a position in the given array is not present in the variation library, it is ignored.
     * @param outfileS the output file name
     * @param customPositions an array of positions to be written
     */
    public void writeBinaryFileS (String outfileS, int[] customPositions) {
        IntArrayList indexList = new IntArrayList();
        for (int i = 0; i < positions.length; i++) {
            int index = Arrays.binarySearch(customPositions, positions[i]);
            if (index < 0) continue;
            indexList.add(i);
        }
        int[] indices = indexList.toIntArray();
        int missingN = customPositions.length-indices.length;
        if (missingN > 0) {
            StringBuilder sb = new StringBuilder();
            sb.append(missingN).append(" positions are missing in the variation library");
            System.out.println(sb.toString());
        }
        try {
            DataOutputStream dos = IOUtils.getBinaryGzipWriter(outfileS);
            dos.writeShort(chrom);
            dos.writeInt(binStart);
            dos.writeInt(indices.length);
            for (int i = 0; i < indices.length; i++) {
                dos.writeInt(positions[indices[i]]);
                dos.writeByte((byte) this.allelePacks[indices[i]].length);
                for (int j = 0; j < allelePacks[indices[i]].length; j++) {
                    dos.writeInt(this.allelePacks[indices[i]][j]);
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
        sb.append(indices.length).append(" polymorphic sites are written to ").append(outfileS);
        System.out.println(sb.toString());
    }

    /**
     * Reads a binary file storing the counts of each allele for each site of a single individual.
     * The file should be in the format of gzip compressed binary file.
     * The content of the file is as follows:
     * 1. chrom (short)
     * 2. binStart (int)
     * 3. positionNum (int)
     * 4. For each position:
     *      position (int)
     *      alleleNum (byte): the number of alleles, negative if missing
     *      alleleCounts (int[]): the counts of each allele
     * @param infileS the path of the input file
     */
    private void readBinaryFileS (String infileS) {
        try {
            DataInputStream dis = IOUtils.getBinaryGzipReader(infileS);
            chrom = dis.readShort();
            binStart = dis.readInt();
            int positionNum = dis.readInt();
            positions = new int[positionNum];
            allelePacks = new int[positionNum][];
            for (int i = 0; i < positionNum; i++) {
                positions[i] = dis.readInt();
                int alleleNum = dis.readByte();
                allelePacks[i] = new int[alleleNum];
                for (int j = 0; j < alleleNum; j++) {
                    allelePacks[i][j] = dis.readInt();
                }
            }
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * Writes the variation library to a text file.
     * The format of the file is as follows:
     * Position\tAlts\tIndelLength
     * The Position column is the position of the variation.
     * The Alts column is the list of alleles at the position.
     * The IndelLength column is the list of indel lengths at the position.
     * @param outfileS the path of the output file
     */
    public void writeTextFileS (String outfileS) {
        try {
            new File(outfileS).getParentFile().mkdirs();
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Position\tAlts\tIndelLength");
            bw.newLine();
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < positions.length; i++) {
                sb.setLength(0);
                sb.append(positions[i]).append("\t");
                for (int j = 0; j < this.allelePacks[i].length; j++) {
                    sb.append(AllelePackageF3.getAlleleBase(allelePacks[i][j])).append(",");
                }
                sb.deleteCharAt(sb.length()-1).append("\t");
                for (int j = 0; j < this.allelePacks[i].length; j++) {
                    int indelLength = AllelePackageF3.getIndelLength(allelePacks[i][j]);
                    sb.append(indelLength).append(",");
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

    /**
     * Compare two VariationLibraryF3 objects according to their binStart.
     * @param o another VariationLibraryF3 object
     * @return -1 if this.binStart < o.binStart, 1 if this.binStart > o.binStart, 0 if this.binStart == o.binStart
     */
    @Override
    public int compareTo(VariationLibraryF3 o) {
        if (this.binStart < o.binStart) return -1;
        else if (this.binStart > o.binStart) return 1;
        return 0;
    }
}
