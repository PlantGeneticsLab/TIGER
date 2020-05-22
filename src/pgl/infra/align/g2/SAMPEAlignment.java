/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.align.g2;

import pgl.infra.utils.IOUtils;

import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

/**
 * Class holding paired-end alignment results of next-gen sequencing reads from BWA-MEM
 * Supplementary alignments are ignored while reading SAM from BWA-MEM
 * Can be sorted by position of r1 read
 * @author feilu
 */
public class SAMPEAlignment {
    PEAlignRecord[] pars = null;
    byte sortType = 0;

    /**
     * Construct an object
     */
    public SAMPEAlignment () {

    }

    /**
     * Construct an object from an unsorted SAM file
     * @param inputFileS
     */
    public SAMPEAlignment (String inputFileS) {
        this.readFromBWAMEM(inputFileS);
    }

    /**
     * Read an input SAM file, which needs to be unsorted
     * @param inputFileS
     */
    public void readFromBWAMEM (String inputFileS) {
        System.out.println("Reading SAM format alignment (BWA-MEM) from: " + inputFileS);
        List<SEAlignRecord> rList = new ArrayList();
        try {
            BufferedReader br;
            if (inputFileS.endsWith(".gz")) {
                br = IOUtils.getTextGzipReader(inputFileS);
            } else {
                br = IOUtils.getTextReader(inputFileS);
            }
            while(br.readLine().startsWith("@PG")==false) {};
            String inputStr = null;
            SEAlignRecord sar = null;
            int cnt = 0;
            while((inputStr = br.readLine())!=null) {
                sar = SAMUtils.getSEAlignRecord(inputStr);
                if (sar.isSupplementaryAlignment()) continue; //did not find other type of second-alignment
                rList.add(sar);
                cnt++;
                if (cnt%500000 == 0) System.out.println("Read in " + String.valueOf(cnt) + " lines");
            }
            br.close();
            pars = new PEAlignRecord[cnt/2];
            for (int i = 0; i < cnt; i+=2) {
                pars[i/2] = new PEAlignRecord(rList.get(i), rList.get(i+1));
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        this.sortByPosition();
        System.out.println("SAMPEAlignment object has "+String.valueOf(pars.length) + " alignment records");
    }

    /**
     * Sort the object by alignment position of r1
     */
    public void sortByPosition () {
        Arrays.sort(this.pars);
        sortType = 0;
    }

    /**
     * Return number of alignment of PE reads
     * @return
     */
    public int getAlignmentNumber () {
        return this.pars.length;
    }

    /**
     * Return an alignment record by index
     * @param index
     * @return
     */
    public PEAlignRecord getAlignmentRecord (int index) {
        return pars[index];
    }

    /**
     * Return all the hits, sorted by string
     * @return
     */
    public String[] getHits () {
        HashSet<String> hitSet = new HashSet();
        for (int i = 0; i < this.getAlignmentNumber(); i++) {
            if (this.getAlignmentRecord(i).getR1AlignmentRecord().isMapped()) {
                hitSet.add(this.getAlignmentRecord(i).getR1AlignmentRecord().getHit());
            }
            if (this.getAlignmentRecord(i).getR2AlignmentRecord().isMapped()) {
                hitSet.add(this.getAlignmentRecord(i).getR2AlignmentRecord().getHit());
            }
        }
        String[] hits = hitSet.toArray(new String[hitSet.size()]);
        Arrays.sort(hits);
        return hits;
    }

    /**
     * Note: Need to be sorted by hit and pos first, otherwise program quit
     * @param hit
     * @return
     */
    public int getAlignmentNumberByHit (String hit) {
        if (this.sortType != 0) {
            System.out.println("Alignment should be sorted by hit and pos first. Program quits");
            System.exit(0);
        }
        return getAlignmentEndIndexByHit(hit) - getAlignmentStartIndexByHit(hit);
    }

    /**
     * Return -1 if the hit doesn't exist, inclusive
     * @param hit
     * @return
     */
    public int getAlignmentStartIndexByHit (String hit) {
        if (this.sortType != 0) {
            System.out.println("Alignment should be sorted by hit and pos first. Program quits");
            System.exit(0);
        }
        int index = this.getAlignmentIndexByHitPos(hit, Integer.MIN_VALUE);
        if (index < 0) {
            index = -index -1;
            if (index < this.getAlignmentNumber() && pars[index].getR1AlignmentRecord().getHit().equals(hit)) return index;
            else return -1;
        }
        else {
            while (index > 0 && pars[index-1].getR1AlignmentRecord().getHit().equals(hit)) {
                index--;
            }
            return index;
        }
    }

    /**
     * Return the index of the last alignment of a hit, exclusive.
     * Return -1 if the hit doesn't exist
     * @param hit
     * @return
     */
    public int getAlignmentEndIndexByHit (String hit) {
        if (this.sortType != 0) {
            System.out.println("Alignment should be sorted by hit and pos first. Program quits");
            System.exit(0);
        }
        int index = this.getAlignmentIndexByHitPos(hit, Integer.MAX_VALUE);
        if (index < 0) {
            index = - index - 2;
            if (pars[index].getR1AlignmentRecord().getHit().equals(hit)) return index+1;
            else return -1;
        }
        else {
            while ((index+1) < this.getAlignmentNumber() && pars[index+1].getR1AlignmentRecord().getHit().equals(hit)) {
                index++;
            }
            return index+1;
        }
    }

    /**
     * Return alignment index
     * Note: Need to be sorted by hit and pos first, otherwise program quits
     * @param hit
     * @param pos
     * @return
     */
    public int getAlignmentIndexByHitPos (String hit, int pos) {
        if (this.sortType != 1) {
            System.out.println("Alignment should be sorted by hit and pos first. Program quits");
            System.exit(0);
        }
        PEAlignRecord query = new PEAlignRecord().setR1AlignmentRecord(new SEAlignRecord().setHit(hit).setStartPos(pos).build()).build();
        return Arrays.binarySearch(pars, query);
    }
}
