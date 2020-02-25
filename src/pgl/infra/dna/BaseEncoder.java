/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.dna;

import com.koloboke.collect.map.hash.HashByteByteMap;
import com.koloboke.collect.map.hash.HashByteByteMaps;
import gnu.trove.list.array.TByteArrayList;




/**
 * Base encoder for A(00), C(01), G(10), T(11). Bases are compressed into 2 bits and loaded into long, int, and short primitives.
 * The primitives are padded with ployA (0000...) at the end, if sequence length is less than 32, 16, or 8.
 * @author fl262
 */
public class BaseEncoder {
    /**
     * Sequence length in long using 2 bits coding
     */
    public static final int longChunkSize = 32;
    /**
     * Sequence length in int using 2 bits coding
     */
    public static final int intChunkSize = 16;
    /**
     * Sequence length in short using 2 bits coding
     */
    public static final int shortChunkSize = 8;
    /**
     * Nucleotide bases in char
     */
    public static final char[] bases = {'A', 'C', 'G', 'T'};
    /**
     * Nucleotide bases in AscII code
     */
    public static final byte[] baseAscIIs = {65, 67, 71, 84};
    /**
     * Nucleotide bases in byte
     */
    public static final byte[] baseBytes = {0, 1, 2, 3};
    /**
     * Converter from AscII to base byte
     */
    public static final HashByteByteMap ascIIBaseByteMap = getAscIIBaseByteMap();
    
    /**
     * Build a base byte converter to convert AscII byte following the base encoding rules
     * A(00000000), C(00000001), G(00000010), T(0000000011), others(00000100)
     * @return 
     */
    public static HashByteByteMap getAscIIBaseByteMap() {
        int size = 128;
        byte[] key = new byte[size];
        byte[] value = new byte[size];
        for (int i = 0; i < key.length; i++) {
            key[i] = (byte)i;
            value[i] = (byte)4;
        }
        value[65] = baseBytes[0];
        value[97] = baseBytes[0];
        value[67] = baseBytes[1];
        value[99] = baseBytes[1];
        value[71] = baseBytes[2];
        value[103] = baseBytes[2];
        value[84] = baseBytes[3];
        value[116] = baseBytes[3];
        HashByteByteMap ascIIByteMap = HashByteByteMaps.newImmutableMap(key, value);
//        HashByteByteMap ascIIByteMap = HashByteByteMaps.newMutableMap();
//        for (int i = 0; i < 128; i++) {
//            ascIIByteMap.put((byte)i, (byte)4);
//        }
//        ascIIByteMap.put((byte)65, (byte)0);
//        ascIIByteMap.put((byte)97, (byte)0);
//        ascIIByteMap.put((byte)67, (byte)1);
//        ascIIByteMap.put((byte)99, (byte)1);
//        ascIIByteMap.put((byte)71, (byte)2);
//        ascIIByteMap.put((byte)103, (byte)2);
//        ascIIByteMap.put((byte)84, (byte)3);
//        ascIIByteMap.put((byte)116, (byte)3);
        return ascIIByteMap;
    }

    /**
     * Convert to an array of base byte from an array of AscII code
     * @param a an array of bytes of AscII code
     * @return
     */
    public static byte[] convertToBaseByteArray (byte[] a) {
        for (int i = 0; i < a.length; i++) {
            a[i] = convertToBaseByte(a[i]);
        }
        return a;
    }

    /**
     * Convert a base byte from an AscII code
     * @param a a byte of AscII code
     * @return
     */
    public static byte convertToBaseByte (byte a) {
        return ascIIBaseByteMap.get(a);
    }

    /**
     * Return a reverse complementary byte array
     * @param a The byte array must be in 0-3 coding for A, C, G, T, other bases are 4
     * @return 
     */
    public static byte[] getReverseComplementary (byte[] a) {
        byte[] b = new byte[a.length];
        for (int i = 0; i < a.length; i++) {
            int t;
            if (a[i] > 3) {
                t = a[i];
            } 
            else {
                t = (~a[i]<<30)>>>30;
            }
            b[a.length-i-1] = (byte)t;
        }
        return b;
    }

    /**
     * Return a long from a byte array. It is padded with polyA at the end if the length of this array is less than 32
     * If the length of this array is greater than 32, return -1 (polyT)
     * @param b The byte array must be in 0-3 coding for A, C, G, T
     * @return 
     */
    public static long getLongSeqFromByteArray(byte[] b) {
        long v = 0;
        if (b.length > longChunkSize) {
            return -1;
        }
        for (int i = 0; i < b.length; i++) {
            v = (v << 2) + b[i];
        }
        v = (v << (2*(longChunkSize-b.length)));
        return v;
    }
    
    /**
     * Return a long from a subset byte array. It is padded with polyA at the end if the length of this subset array is less than 32
     * If the length of this subset array is greater than 32, return -1 (polyT)
     * @param b The byte array must be in 0-3 coding for A, C, G, T
     * @param startIndex inclusive
     * @param endIndex exclusive
     * @return 
     */
    public static long getLongSeqFromSubByteArray(byte[] b, int startIndex, int endIndex) {
        int length = endIndex - startIndex;
        long v = 0;
        if (length > longChunkSize) {
            return -1;
        }
        for (int i = startIndex; i < endIndex; i++) {
            v = (v << 2) + b[i];
        }
        v = (v << (2*(longChunkSize-length)));
        return v;
    }
    
    /**
     * Return an int array from a long array
     * @param val
     * @return 
     */
    public static int[] getIntSeqsFromLongSeqs (long[] val) {
        int[] v = new int[val.length*2];
        for (int i = 0; i < val.length; i++) {
            int[] cv = getIntSeqsFromLongSeq(val[i]);
            v[i*2] = cv[0];
            v[i*2+1] = cv[1];
        }
        return v;
    }
    
    /**
     * Return an int array from a long
     * @return 
     */
    public static int[] getIntSeqsFromLongSeq (long val) {
        int[] v = new int[2];
        v[0] = (int)(val>>>longChunkSize);
        v[1] = (int)val;
        return v;
    }
    
    /**
     * Return an int from a byte array. It is padded with polyA at the end if the length of this array is less than 16
     * If the length of this array is greater than 16, return -1 (polyT)
     * @param b The byte array must be in 0-3 coding for A, C, G, T
     * @return 
     */
    public static int getIntSeqFromByteArray(byte[] b) {
        int v = 0;
        if (b.length > intChunkSize) {
            return -1;
        }
        for (int i = 0; i < b.length; i++) {
            v = (v << 2) + b[i];
        }
        v = (v << (2*(intChunkSize-b.length)));
        return v;
    }
    
    /**
     * Return an int from a subset byte array. It is padded with polyA at the end if the length of this subset array is less than 16
     * If the length of this subset array is greater than 16, return -1 (polyT)
     * @param b The byte array must be in 0-3 coding for A, C, G, T
     * @param startIndex
     * @param endIndex
     * @return 
     */
    public static int getIntSeqFromSubByteArray(byte[] b, int startIndex, int endIndex) {
        int length = endIndex - startIndex;
        int v = 0;
        if (length > intChunkSize) {
            return -1;
        }
        for (int i = startIndex; i < endIndex; i++) {
            v = (v << 2) + b[i];
        }
        v = (v << (2*(intChunkSize-length)));
        return v;
    }
    
    /**
     * Return an int from a byte array. It is padded with polyA at the end if the length of this array is less than 8
     * If the length of this array is greater than 8, return -1 (polyT)
     * @param b The byte array must be in 0-3 coding for A, C, G, T
     * @return 
     */
    public static short getShortSeqFromByteArray(byte[] b) {
        int v = 0;
        if (b.length > shortChunkSize) {
            return -1;
        }
        for (int i = 0; i < b.length; i++) {
            v = (v << 2) + b[i];
        }
        v = (v << (2*(shortChunkSize-b.length)));
        return (short)v;
    }
    
    /**
     * Return a short from a subset byte array. It is padded with polyA at the end if the length of this subset array is less than 8
     * If the length of this subset array is greater than 8, return -1 (polyT)
     * @param b The byte array must be in 0-3 coding for A, C, G, T
     * @param startIndex
     * @param endIndex
     * @return 
     */
    public static short getShortSeqFromSubByteArray(byte[] b, int startIndex, int endIndex) {
        int length = endIndex - startIndex;
        int v = 0;
        if (length > shortChunkSize) {
            return -1;
        }
        for (int i = startIndex; i < endIndex; i++) {
            v = (v << 2) + b[i];
        }
        v = (v << (2*(shortChunkSize-length)));
        return (short)v;
    }
    
    /**
     * Return a reverse complementary long. It is padded with polyA at the end if the length of this array is less than 32
     * If the length of this array is greater than 32 or less than 0, return -1 (polyT)
     * @param seq
     * @param seqLength The actual length of the seq
     * @return 
     */
    public static long getLongReverseComplement(long seq, int seqLength) {
        if (seqLength<0 || seqLength>longChunkSize) return -1;
        long rev = 0;
        long mask = 3;
        seq = ~seq;
        seq = seq >> ((longChunkSize-seqLength)*2);
        for (int i = 0; i < seqLength; i++) {
            rev = (rev << 2) + (seq & mask);
            seq = seq >> 2;
        }
        rev = rev << ((longChunkSize-seqLength)*2);
        return rev;
    }
    
    /**
     * Return a reverse complementary int. It is padded with polyA at the end if the length of this array is less than 16
     * If the length of this array is greater than 16 or less than 0, return -1 (polyT)
     * @param seq
     * @param seqLength The actual length of the seq
     * @return 
     */
    public static int getIntReverseComplement(int seq, int seqLength) {
        if (seqLength<0 || seqLength>intChunkSize) return -1;
        int rev = 0;
        int mask = 3;
        seq = ~seq;
        seq = seq >> ((intChunkSize-seqLength)*2);
        for (int i = 0; i < seqLength; i++) {
            rev = (rev << 2) + (seq & mask);
            seq = seq >> 2;
        }
        rev = rev << ((intChunkSize-seqLength)*2);
        return rev;
    }
    
    /**
     * Return a reverse complementary int. It is padded with polyA at the end if the length of this array is less than 8
     * If the length of this array is greater than 8 or less than 0, return -1 (polyT)
     * @param seq
     * @param seqLength The actual length of the seq
     * @return 
     */
    public static short getShortReverseComplement(int seq, int seqLength) {
        if (seqLength<0 || seqLength>shortChunkSize) return -1;
        int rev = 0;
        int mask = 3;
        seq = ~seq;
        seq = seq >> ((shortChunkSize-seqLength)*2);
        for (int i = 0; i < seqLength; i++) {
            rev = (rev << 2) + (seq & mask);
            seq = seq >> 2;
        }
        rev = rev << ((shortChunkSize-seqLength)*2);
        return (short)rev;
    }
    
    /**
     * Returns the number of bp differences between two 2-bit encoded longs.
     * Maximum divergence is used to save time when only interested in very similar 
     * sequences.
     * @param seq1 2-bit encoded sequence
     * @param seq2 2-bit encoded sequence
     * @param maxDivergence threshold for counting divergence upto
     * @return count of the divergence (above the maxDivergence, chunkSize is returned)
     */
    public static byte getSeqDifferences(long seq1, long seq2, int maxDivergence) {
        long mask = 3;
        byte cnt = 0;
        long diff = seq1 ^ seq2;
        for (int x = 0; x < longChunkSize && cnt <= maxDivergence; x++) {
            if ((diff & mask) > 0) {
                cnt++;
            }
            diff = diff >> 2;
            // System.out.println("v = " + v);
        }
        if (cnt > maxDivergence) {
            cnt = (byte) longChunkSize;
        }
        // if(x<(chunkSize-1)) cnt=(byte)chunkSize;  //if didn't get to the end of the sequence set to maximum
        return cnt;
    }

    
    /**
     * Returns the number of bp differences between two 2-bit encoded longs.
     * @param seq1 2-bit encoded sequence
     * @param seq2 2-bit encoded sequence
     * @return count of the divergence 
     */
    public static byte getSeqDifferences(long seq1, long seq2) {
        long mask = 3;
        byte cnt = 0;
        long diff = seq1 ^ seq2;
        for (int x = 0; x < longChunkSize; x++) {
            if ((diff & mask) > 0) {
                cnt++;
            }
            diff = diff >> 2;
            // System.out.println("v = " + v);
        }
        return cnt;
    }
    
    public static String getSequenceFromLongs (long[] val) {
        StringBuilder sb = new StringBuilder(val.length*longChunkSize+1);
        for (int i = 0; i < val.length; i++) {
            sb.append(getSequenceFromLong(val[i]));
        }
        return sb.toString();
    }
    
    /**
     * Return DNA sequence from long
     * @param val
     * @return 
     */
    public static String getSequenceFromLong(long val) {
    	StringBuilder seq = new StringBuilder(longChunkSize + 1);
    	long mask = 3L << 62;
    	for (int i = 0; i < longChunkSize; i++) {
    		byte base = (byte) (((val & mask) >>> 62));
    		seq.append(bases[base]);          
    		val = val << 2;
    	}
    	return seq.toString();
    }
    
    /**
     * Return DNA sequence from long
     * @param val
     * @return 
     */
    public static String getSequenceFromInt(int val) {
    	StringBuilder seq = new StringBuilder(intChunkSize + 1);
    	int mask = 3 << 30;
    	for (int i = 0; i < intChunkSize; i++) {
    		byte base = (byte) (((val & mask) >>> 30));
    		seq.append(bases[base]);          
    		val = val << 2;
    	}
    	return seq.toString();
    }
    
    /**
     * Return DNA sequence from long
     * @param val
     * @return 
     */
    public static String getSequenceFromShort(int val) {
    	StringBuilder seq = new StringBuilder(shortChunkSize + 1);
    	int mask = 3 << 14;
    	for (int i = 0; i < shortChunkSize; i++) {
    		byte base = (byte) (((val & mask) >>> 14));
    		seq.append(bases[base]);          
    		val = val << 2;
    	}
    	return seq.toString();
    }
    
    public static byte[] getBaseByteArrayFromLong(long val) {
        byte[] array = new byte[longChunkSize];
        long mask = 3L << 62;
    	for (int i = 0; i < longChunkSize; i++) {
            byte base = (byte) (((val & mask) >>> 62));
            array[i] = baseBytes[base];
            val = val << 2;
    	}
        return array;
    }
    
    public static byte[] getBaseByteArrayFromLongs(long[] val) {
        TByteArrayList byteList = new TByteArrayList();
        for (int i = 0; i < val.length; i++) {
            byteList.addAll(getBaseByteArrayFromLong(val[i]));
        }
        return byteList.toArray();
    }
}
