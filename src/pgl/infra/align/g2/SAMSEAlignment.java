/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.align.g2;

import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;

import pgl.infra.utils.IOUtils;

/**
 * Class holding single-end alignment results of next-gen sequencing reads from BWA-MEM
 * Supplementary alignments are ignored while reading SAM from BWA-MEM
 * Can be sorted by position
 * @author feilu
 */
public class SAMSEAlignment {
    SEAlignRecord[] sars = null;
    /**0 by position*/
    byte sortType = 0;
    
    public SAMSEAlignment () {}

    public SAMSEAlignment (String inputFileS) {
        this.readFromBWAMEM(inputFileS);
    }
    
    /**
     * 
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
                br = br = IOUtils.getTextReader(inputFileS);
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
            sars = rList.toArray(new SEAlignRecord[rList.size()]);
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        this.sortByPosition();
        System.out.println("SAMSEAlignment object has "+String.valueOf(sars.length) + " alignment records");
    }
    
    /**
     * Return number of alignment
     * @return 
     */
    public int getAlignmentNumber () {
        return sars.length;
    }
    
    /**
     * Return non-redundant names of queries
     * @return 
     */
    public String[] getQuerys () {
        String[] queries = new String[this.getAlignmentNumber()];
        for (int i = 0; i < queries.length; i++) {
            queries[i] = sars[i].query;
        }
        return queries;
    }
    
    /**
     * Return non-redundant and sorted names of hits
     * @return 
     */
    public String[] getHits () {
        HashSet<String> hitSet = new HashSet();
        for (int i = 0; i < this.getAlignmentNumber(); i++) {
            if (sars[i].hit.equals("")) continue;
            hitSet.add(sars[i].hit);
        }
        String[] hits = hitSet.toArray(new String[hitSet.size()]);
        Arrays.sort(hits);
        return hits;
    }
    
    /**
     * Return an alignment record by index
     * @param index
     * @return 
     */
    public SEAlignRecord getAlignmentRecord (int index) {
        return sars[index];
    }
    
    /**
     * Sort alignments by start positions
     */
    public void sortByPosition () {
        Arrays.sort(sars);
        sortType = 0;
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
     * Retuen the index of the first alignment of a hit, inclusive
     * @param hit
     * @return -1 if the hit doesn't exist
     */
    public int getAlignmentStartIndexByHit (String hit) {
        if (this.sortType != 0) {
            System.out.println("Alignment should be sorted by hit and pos first. Program quits");
            System.exit(0);
        }
        int index = this.getAlignmentIndexByHitPos(hit, Integer.MIN_VALUE);
        if (index < 0) {
            index = -index -1;
            if (index < this.getAlignmentNumber() && sars[index].hit.equals(hit)) return index;
            else return -1;
        }
        else {
            while (index > 0 && sars[index-1].hit.equals(hit)) {
                index--;
            }
            return index;
        }
    }
    
    /**
     * Return the index of the last alignment of a hit, exclusive.
     * Return -1 if the hit doesn't exist
     * Note: Need to be sorted by hit first
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
            if (sars[index].hit.equals(hit)) return index+1;
            else return -1;
        }
        else {
            while ((index+1) < this.getAlignmentNumber() && sars[index+1].hit.equals(hit)) {
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
        return Arrays.binarySearch(sars, new SEAlignRecord().setHit(hit).setStartPos(pos).build());
    }
}
