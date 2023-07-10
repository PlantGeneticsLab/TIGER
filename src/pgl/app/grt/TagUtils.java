/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.app.grt;

import com.koloboke.collect.map.hash.HashByteByteMap;
import pgl.infra.dna.BaseEncoder;
import java.util.Comparator;

/**
 *
 * @author feilu
 */
public class TagUtils {
    static String polyA = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    public static TagComparator tagCom = new TagComparator();
    
    static long[] getTagFromReads (String readR1, String readR2, HashByteByteMap ascIIByteMap, int tagLengthInLong) {
        int setReadLength = tagLengthInLong*BaseEncoder.longChunkSize;
        long[] tag = new long[tagLengthInLong*2];
        StringBuilder sb = new StringBuilder(readR1);
        if (sb.length()< setReadLength) {
            sb.append(polyA);
            readR1 = sb.substring(0, setReadLength);
        }
        sb = new StringBuilder(readR2);
        if (sb.length()<setReadLength) {
            sb.append(polyA);
            readR2 = sb.substring(0, setReadLength);
        }
        byte[] bArray = readR1.getBytes();
        for (int i = 0; i < bArray.length; i++) {
            bArray[i] = ascIIByteMap.get(bArray[i]);
        }
        for (int i = 0; i < tagLengthInLong; i++) {
            tag[i] = BaseEncoder.getLongFromSubBaseCodingArray(bArray, i*BaseEncoder.longChunkSize, (i+1)*BaseEncoder.longChunkSize);
        }
        bArray = readR2.getBytes();
        for (int i = 0; i < bArray.length; i++) {
            bArray[i] = ascIIByteMap.get(bArray[i]);
        }
        for (int i = 0; i < tagLengthInLong; i++) {
            tag[i+tagLengthInLong] = BaseEncoder.getLongFromSubBaseCodingArray(bArray, i*BaseEncoder.longChunkSize, (i+1)*BaseEncoder.longChunkSize);
        }
        return tag;
    }
    
    static String[] getReadsFromTag (long[] tag, byte r1Length, byte r2Length) {
        String[] reads = new String[2];
        int tagLengthInLong = tag.length/2;
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < tagLengthInLong; i++) {
            sb.append(BaseEncoder.getSequenceFromLong(tag[i]));
        }
        reads[0] = sb.toString().substring(0, r1Length);
        sb = new StringBuilder();
        for (int i = 0; i < tagLengthInLong; i++) {
            sb.append(BaseEncoder.getSequenceFromLong(tag[i+tagLengthInLong]));
        }
        reads[1] = sb.toString().substring(0, r2Length);
        return reads;
    }
    
    static int getGroupIndexFromTag (long[] tag, int offSet, int groupIdentifierLength) {
        long a = tag[0];
        int move = 0;
        long mask = 3;
        int value = 0;
        int index = 0;
        int dOff = offSet * 2;
        for (int i = 0; i < groupIdentifierLength; i++) {
            move = 64 - dOff - i * 2 - 2;
            value = (int)((a >>> move) & mask);
            if (i == 0) {
                index += value;
            }
            else {
                index += Math.pow(4, i) * value;
            }
        }
        return index;
    }
}

class TagComparator implements Comparator <long[]> {
    @Override
    public int compare(long[] o1, long[] o2) {
        for (int i = 0; i < o1.length; i++) {
            if (o1[i] < o2[i]) {
                return -1;
            }
            if (o1[i] > o2[i]) {
                return 1;
            }
        }
        return 0;
    }        
}