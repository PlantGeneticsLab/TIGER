/*
 * SAMUtils
 */
package pgl.infra.align.g2;

import pgl.app.grt.AlleleInfo;
import pgl.app.grt.SNPCounts;
import pgl.infra.dna.allele.AlleleEncoder;
import pgl.infra.dna.snp.SNPOld;
import gnu.trove.list.array.TByteArrayList;
import gnu.trove.list.array.TIntArrayList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.Dyad;

/**
 *
 * @author Fei Lu
 */
public class SAMUtils {
    /**
     * CIGAR consume operators are '=', 'D', 'M', 'N', 'X'
     */
    public static final byte[] consumeCigarOPByte = {61, 68, 77, 78, 88};
    
    
    /**
     * Return if the seq has multiple segments
     * @param flag
     * @return 
     */
    public static boolean isHavingMultipleSegments(int flag) {
        if (getBinaryValueAtNthBit(flag, 1) == 1) return true;
        return false;
    }
    
    /**
     * Return if a read pair is properly aligned
     * @param flag
     * @return 
     */
    public static boolean isReadsMappedInProperPair(int flag) {
        if (getBinaryValueAtNthBit(flag, 2) == 1) return true;
        return false;
    }
    
    /**
     * Return if the seq is unmapped
     * @param flag
     * @return 
     */
    public static boolean isUnmapped (int flag) {
        if (getBinaryValueAtNthBit(flag, 3) == 1) return true;
        return false;
    }
    
    /**
     * Return if the mate of the seq is unmapped
     * @param flag
     * @return 
     */
    public static boolean isMateUnmapped (int flag) {
        if (getBinaryValueAtNthBit(flag, 4) == 1) return true;
        return false;
    }
    
    /**
     * Return if the seq is reversely aligned
     * @param flag
     * @return 
     */
    public static boolean isReverseAligned (int flag) {
        if (getBinaryValueAtNthBit(flag, 5) == 1) return true;
        return false;
    }
    
    /**
     * Return if the mate of the seq is reversely aligned
     * @param flag
     * @return 
     */
    public static boolean isMateReverseAligned (int flag) {
        if (getBinaryValueAtNthBit(flag, 6) == 1) return true;
        return false;
    }
    
    /**
     * Return if the seq the first read in a pair
     * @param flag
     * @return 
     */
    public static boolean isFirstSegmentInPair (int flag) {
        if (getBinaryValueAtNthBit(flag, 7) == 1) return true;
        return false;
    }
    
    /**
     * Return if the seq the second read in a pair
     * @param flag
     * @return 
     */
    public static boolean isSecondSegmentInPair (int flag) {
        if (getBinaryValueAtNthBit(flag, 8) == 1) return true;
        return false;
    }
    
    /**
     * Return if the seq is secondary alignment
     * @param flag
     * @return 
     */
    public static boolean isSecondaryAlignment (int flag) {
        if (getBinaryValueAtNthBit(flag, 9) == 1) return true;
        return false;
    }
    
    /**
     * Return if the seq is a PCR duplicate
     * @param flag
     * @return 
     */
    public static boolean isPCRDuplicates (int flag) {
        if (getBinaryValueAtNthBit(flag, 11) == 1) return true;
        return false;
    }
    
    /**
     * Return if the alignment is a supplementary alignment
     * @param flag
     * @return 
     */
    public static boolean isSupplementaryAlignment (int flag) {
        if (getBinaryValueAtNthBit(flag, 12) == 1) return true;
        return false;
    }
    
    private static int getBinaryValueAtNthBit (int flag, int bitPosRightMost) {
        int shift = (bitPosRightMost - 1);
        return (flag>>shift)&1;
    }
    
    /**
     * Return a {@link SEAlignRecord} object from a SAM alignment record
     * @param inputStr
     * @return 
     */
    public static SEAlignRecord getSEAlignRecord (String inputStr) {
        return getSEAlignRecord (PStringUtils.fastSplit(inputStr));
    }
    
    /**
     * Return a {@link SEAlignRecord} object from elements of a SAM alignment record
     * @param l
     * @return 
     */
    private static SEAlignRecord getSEAlignRecord (List<String> l) {
        short flag = Short.parseShort(l.get(1));
        String hit = "";
        int startPos = Integer.MIN_VALUE;
        int endPos = Integer.MIN_VALUE;
        short mapQ = Short.MIN_VALUE;
        short alnMatchNumber = Short.MIN_VALUE;
        short editDistance = Short.MIN_VALUE;
        if (l.get(5).equals("*")) {
            //continue;
        }
        else {
            hit = l.get(2);
            startPos = Integer.parseInt(l.get(3));
            String cigar = l.get(5);
            Dyad<TByteArrayList, TIntArrayList> cigarOpPosIndex = getCigarOPAndPosIndex (cigar);
            endPos = getEndPos(cigar, cigarOpPosIndex, startPos);
            mapQ = Short.parseShort(l.get(4));
            alnMatchNumber  = getAlignMatchNumberInCigar (cigar, cigarOpPosIndex);
            editDistance = Short.valueOf(l.get(11).split(":")[2]);
        }
        SEAlignRecord sar = new SEAlignRecord (l.get(0), hit, startPos, endPos, flag, mapQ, alnMatchNumber, editDistance);
        return sar;
    }
    
    /**
     * Return the end position (exclusive) of reference from an alignment
     * Return Integer.MIN_VALUE if the query is not aligned, in which cigar is *
     * @param cigar
     * @param startPos
     * @return 
     */
    public static int getEndPos (String cigar, int startPos) {
        Dyad<TByteArrayList, TIntArrayList> opPosIndex = getCigarOPAndPosIndex(cigar);
        return getEndPos(cigar, opPosIndex, startPos);
    }
    
    /**
     * Return the end position (exclusive) of reference from an alignment
     * Return Integer.MIN_VALUE if the query is not aligned, in which cigar is *
     * @param cigar
     * @param opPosIndex
     * @param startPos
     * @return 
     */
    private static int getEndPos (String cigar, Dyad<TByteArrayList, TIntArrayList> opPosIndex, int startPos) {
        if (opPosIndex == null) return Integer.MIN_VALUE;
        byte[] op = opPosIndex.getFirstElement().toArray();
        int[] posIndex = opPosIndex.getSecondElement().toArray();
        int endPos = startPos - 1;
        for (int i = 0; i < op.length; i++) {
            int index = Arrays.binarySearch(consumeCigarOPByte, op[i]);
            if (index < 0) continue;
            if (i == 0) {
                endPos += Integer.valueOf(cigar.substring(0, posIndex[i]));
            }
            else {
                endPos += Integer.valueOf(cigar.substring(posIndex[i-1]+1, posIndex[i]));
            }
        }
        return endPos+1;
    }
    
    /**
     * Return total length of alignment match in CIGAR ('M'), including both sequence match and mismatch
     * @param cigar
     * @param cigarOpPosIndex
     * @return 
     */
    private static short getAlignMatchNumberInCigar (String cigar, Dyad<TByteArrayList, TIntArrayList> cigarOpPosIndex) {
        byte[] op = cigarOpPosIndex.getFirstElement().toArray();
        int[] posIndex = cigarOpPosIndex.getSecondElement().toArray();
        int len = 0;
        for (int i = 0; i < op.length; i++) {
            if (op[i] != 77) continue;
            if (i == 0) {
                len += Integer.valueOf(cigar.substring(0, posIndex[i]));
            }
            else {
                len += Integer.valueOf(cigar.substring(posIndex[i-1]+1, posIndex[i]));
            }
        }
        return (short)len;
    }
    
    /**
     * Return operators and their position index of CIGAR in a {@link Dyad} format
     * @param cigar
     * @return 
     */
    private static Dyad<TByteArrayList, TIntArrayList> getCigarOPAndPosIndex (String cigar) {
        if (cigar.startsWith("*")) return null;
        TByteArrayList opList = new TByteArrayList();
        TIntArrayList posList = new TIntArrayList();
        byte[] cigarB = cigar.getBytes();
        for (int i = 0; i < cigar.length(); i++) {
            if (cigarB[i] > 64) {
                opList.add(cigarB[i]);
                posList.add(i);
            }
        }
        return new Dyad<TByteArrayList, TIntArrayList> (opList, posList);
    }
    
    /**
     * Return a list of elements from a SAM alignment record
     * @param inputStr
     * @return 
     */
    public static List<String> getAlignElements (String inputStr) {
        return PStringUtils.fastSplit(inputStr);
    }
    
    /**
     * Return a list of called SNPs from a SAM alignment record
     * @param inputStr
     * @param mapQThresh
     * @return 
     */
    public static List<SNPOld> getVariants (String inputStr, int mapQThresh) {
        List<String> l = PStringUtils.fastSplit(inputStr);
        return getVariants(l, mapQThresh);
    }

    /**
     * Working on it
     * @param l
     * @param mapQThresh
     * @param sc
     * @param end
     * @return
     */
    public static List<AlleleInfo> getAlleles3 (List<String> l, int mapQThresh, SNPCounts sc, int end) {
        List<AlleleInfo> tagAlleleList = new ArrayList();
        String cigar = l.get(5);
        if (cigar.startsWith("*")) return null;
        if (Integer.parseInt(l.get(4)) < mapQThresh) return null;
        short chr = (short)Integer.parseInt(l.get(2));
        int startPos = Integer.parseInt(l.get(3));
        Dyad<TByteArrayList, TIntArrayList> cigarOPPosIndex = getCigarOPAndPosIndex(cigar);
        int endPos = getEndPos (cigar, cigarOPPosIndex, startPos);
        int chrIndex = sc.getChrIndex(chr);
        if (chrIndex < 0) return null;
        int sIndex = sc.getSNPIndex(chrIndex, startPos);
        int eIndex = sc.getSNPIndex(chrIndex, endPos);
        if (sIndex == eIndex) return null;
        if (sIndex < 0) {
            sIndex = -sIndex-1;
            if (sIndex == sc.getSNPNumberOnChromosome(chrIndex)) return null;
        }
        if (eIndex < 0) {
            eIndex = -eIndex-1;
            if (eIndex == 0) return null;
        }
        else {
            eIndex++;
        }
        int[] snpRefPos = new int[eIndex-sIndex];
        for (int i = 0; i < snpRefPos.length; i++) {
            snpRefPos[i] = sc.getPositionOfSNP(chrIndex, i+sIndex);
            AlleleInfo ai = new AlleleInfo (sc.getChromosome(chrIndex), sc.getPositionOfSNP(chrIndex, i+sIndex), (byte)end);
            tagAlleleList.add(ai);
        }
        String seq = l.get(9);
        String md = l.get(12).split(":")[2];
        int seqLength = seq.length();
        boolean ifMinus = SAMUtils.isReverseAligned(Integer.parseInt(l.get(1)));
        TByteArrayList opList = cigarOPPosIndex.getFirstElement();
        TIntArrayList posIndexList = cigarOPPosIndex.getSecondElement();
        int currentRefPos = startPos-1;
        int currentAltPos = 0;
        int frontLength = 0;
        int backLength = 0;
        for (int i = 0; i < opList.size(); i++) {
            int length;
            if (i == 0) {
                length = Integer.parseInt(cigar.substring(0, posIndexList.get(0)));
            }
            else {
                length = Integer.parseInt(cigar.substring(posIndexList.get(i-1)+1, posIndexList.get(i)));
            }
            if (opList.get(i) == 83) {//S
                if (i == 0) {
                    seq = seq.substring(length);
                    frontLength=length;
                }
                else if (i == opList.size()-1) {
                    seq = seq.substring(0, seq.length()-length);
                    backLength = length;
                }
            }
            else if (opList.get(i) == 72 || opList.get(i) == 80) {//H, P
                if (i == 0) {
                    frontLength+=length;
                }
                else if (i == opList.size()-1) {
                    backLength = length;
                }
            }
            else if (opList.get(i) == 73) {//I
                int index = Arrays.binarySearch(snpRefPos, currentRefPos);
                if (index < 0) {
                    
                }
                else {
                    byte seqAllele = 5;
                    byte refAllele = sc.getRefAlleleByteOfSNP(chrIndex, index+sIndex);
                    int alleleIndex = sc.getAltAlleleIndex(chrIndex, index+sIndex, seqAllele);
                    if (alleleIndex < 0) {

                    }
                    else {
                        tagAlleleList.get(index).setRelativePosition((byte)(currentAltPos+frontLength+1));
                        tagAlleleList.get(index).setAllele(seqAllele);
                        byte base = AlleleEncoder.alleleBaseToCodingMap.get(seq.charAt(currentAltPos));
                        tagAlleleList.get(index).setBase(base);
                    }
                }
                currentAltPos+=length;
            }
            else if (opList.get(i) == 68 || opList.get(i) == 78) {//D, N
                int index = Arrays.binarySearch(snpRefPos, currentRefPos+1);
                if (index < 0) {
                    
                }
                else {
                    byte seqAllele = 4;
                    byte refAllele = sc.getRefAlleleByteOfSNP(chrIndex, index+sIndex);
                    int alleleIndex = sc.getAltAlleleIndex(chrIndex, index+sIndex, seqAllele);
                    if (alleleIndex < 0) {

                    }
                    else {
                        tagAlleleList.get(index).setRelativePosition((byte)(currentAltPos+frontLength));
                        tagAlleleList.get(index).setAllele(seqAllele);
                        byte base = AlleleEncoder.alleleBaseToCodingMap.get(seq.charAt(currentAltPos-1));
                        tagAlleleList.get(index).setBase(base);
                    }
                }
                currentRefPos+=length;
            }
            else {//M,77; =, X
                for (int j = 0; j < length; j++) {
                    currentRefPos++;
                    currentAltPos++;
                    int index = Arrays.binarySearch(snpRefPos, currentRefPos);
                    if (index < 0) {
                        index = -index-1;
                        if (index == snpRefPos.length) break;
                        int cha = snpRefPos[index]-currentRefPos-1;
                        if (cha + j> length) {
                            cha = length - j-1;
                            currentRefPos+=cha;
                            currentAltPos+=cha;
                            break;
                        }
                        currentRefPos+=cha;
                        currentAltPos+=cha;
                        j+=cha;
                        continue;
                    }
                    byte seqAllele = AlleleEncoder.alleleBaseToCodingMap.get(seq.charAt(currentAltPos-1));
                    byte refAllele = sc.getRefAlleleByteOfSNP(chrIndex, index+sIndex);
                    int insertionAlleleIndex = sc.getAltAlleleIndex(chrIndex, index+sIndex, (byte)5);
                    if (seqAllele == refAllele) {
                        if (insertionAlleleIndex < 0) {
                            
                        }
                        else {
                            tagAlleleList.get(index).setRelativePosition((byte)(currentAltPos+frontLength));
                            tagAlleleList.get(index).setAllele((byte)5);
                            tagAlleleList.get(index).setBase(seqAllele);
                            continue;
                        }
                    }    
                    tagAlleleList.get(index).setRelativePosition((byte)(currentAltPos+frontLength));
                    tagAlleleList.get(index).setAllele(seqAllele);
                    tagAlleleList.get(index).setBase(seqAllele);
                }
            }
        }
        
        if (ifMinus) {
            for (int i = 0; i < tagAlleleList.size(); i++) {
                byte relaPos = tagAlleleList.get(i).getRelativePosition();
                tagAlleleList.get(i).setRelativePosition((byte)(seqLength-relaPos+1));
                byte base = tagAlleleList.get(i).getBase();
                base = (byte)(~base & 3);
                tagAlleleList.get(i).setBase(base);
            }
        }
        for (int i = 0; i < tagAlleleList.size(); i++) {
            if (tagAlleleList.get(i).getAllele() == Byte.MIN_VALUE) {
                tagAlleleList.remove(i);
                i--;
            }
        }
        return tagAlleleList;
    }
    
//    public static List<AlleleInfo> getAlleles2 (List<String> l, int mapQThresh, SNPCounts sc, int end) {
//        List<AlleleInfo> tagAlleleList = new ArrayList();
//        List<ChrPos> allelePosList = new ArrayList();
//        TByteArrayList alleleList = new TByteArrayList();
//        TByteArrayList alleleRelaPosList = new TByteArrayList();
//        String cigar = l.get(5);
//        if (cigar.startsWith("*")) return null;
//        if (Integer.parseInt(l.get(4)) < mapQThresh) return null;
//        short chr = (short)Integer.parseInt(l.get(2));
//        int startPos = Integer.parseInt(l.get(3));
//        Dyad<TByteArrayList, TIntArrayList> cigarOPPosIndex = getCigarOPAndPosIndex(cigar);
//        int endPos = getEndPos (cigar, cigarOPPosIndex, startPos);
//        int chrIndex = sc.getChrIndex(chr);
//        if (chrIndex < 0) return null;
//        int sIndex = sc.getSNPIndex(chrIndex, startPos);
//        int eIndex = sc.getSNPIndex(chrIndex, endPos);
//        if (sIndex == eIndex) return null;
//        if (sIndex < 0) {
//            sIndex = -sIndex-1;
//            if (sIndex == sc.getSNPNumberOnChromosome(chrIndex)) return null;
//        }
//        if (eIndex < 0) {
//            eIndex = -eIndex-1;
//            if (eIndex == 0) return null;
//        }
//        else {
//            eIndex++;
//        }
////*****copy from getVariants method,  to avoid build cigarOPPosIndex twice, with some changes******        
//        String seq = l.get(9);
//        String md = l.get(12).split(":")[2];
//        int seqLength = seq.length();
//        boolean ifMinus = SAMUtils.isReverseAligned(Integer.parseInt(l.get(1)));
//        TByteArrayList opList = cigarOPPosIndex.getFirstElement();
//        TIntArrayList posIndexList = cigarOPPosIndex.getSecondElement();
//        int frontLength = 0;
//        int backLength = 0;
//        if (opList.get(0) == 83) {//S
//            int length = Integer.parseInt(cigar.substring(0, posIndexList.get(0)));
//            frontLength = length;
//            seq = seq.substring(length);
//            cigar = cigar.substring(posIndexList.get(0)+1);
//            if (opList.get(opList.size()-1) == 83 || opList.get(opList.size()-1) == 72) {
//                cigarOPPosIndex = getCigarOPAndPosIndex(cigar);
//                opList = cigarOPPosIndex.getFirstElement();
//                posIndexList = cigarOPPosIndex.getSecondElement();
//            }
//        }
//        else if (opList.get(0) == 72) {//H
//            int length = Integer.parseInt(cigar.substring(0, posIndexList.get(0)));
//            frontLength = length;
//            cigar = cigar.substring(posIndexList.get(0)+1);
//            if (opList.get(opList.size()-1) == 83 || opList.get(opList.size()-1) == 72) {
//                cigarOPPosIndex = getCigarOPAndPosIndex(cigar);
//                opList = cigarOPPosIndex.getFirstElement();
//                posIndexList = cigarOPPosIndex.getSecondElement();
//            }
//        }
//        if (opList.get(opList.size()-1) == 83) {//S
//            int length = Integer.parseInt(cigar.substring(posIndexList.get(opList.size()-2)+1, posIndexList.get(opList.size()-1)));
//            backLength = length;
//            seq = seq.substring(0, seq.length()-length);
//            cigar = cigar.substring(0, posIndexList.get(opList.size()-2)+1);
//        }
//        else if (opList.get(opList.size()-1) == 72) {//H
//            int length = Integer.parseInt(cigar.substring(posIndexList.get(opList.size()-2)+1, posIndexList.get(opList.size()-1)));
//            backLength = length;
//            cigar = cigar.substring(0, posIndexList.get(opList.size()-2)+1);
//        }
//        List<SNP> snpList = new ArrayList();
//        int currentPos = startPos-1;
//        int currentAltPos = 0;
//        
//        TIntArrayList insertionPosList = new TIntArrayList();
//        TIntArrayList insertionLengthList = new TIntArrayList();
//        
//        if (opList.contains((byte)73)) {//I 18M1I31M1D46M
//            for (int i = 0; i < cigar.length(); i++) {
//                char c = cigar.charAt(i);
//                if (Character.isDigit(c)) {
//                    int j;
//                    for (j = i + 1; j < cigar.length(); j++) {
//                        if (!Character.isDigit(cigar.charAt(j))) {
//                            break;
//                        }
//                    }
//                    int length = Integer.parseInt(cigar.substring(i, j));
//                    char cop = cigar.charAt(j);
//                    if (cop == 'I') {
//                        char alt = 'I';
//                        char ref = seq.charAt(currentAltPos-1);
//                        alleleRelaPosList.add((byte)(currentAltPos+frontLength)); 
//                        snpList.add(new SNPOld(chr, currentPos, ref, alt));
//                        StringBuilder sb = new StringBuilder(seq);
//                        sb.delete(currentAltPos, currentAltPos+length);
//                        seq = sb.toString();
//                        insertionPosList.add(currentAltPos);
//                        insertionLengthList.add(length);
//                    }
//                    else if (cop == 'D' || cop == 'N') {
//                        currentPos += length;
//                    }
//                    else { //M,=,X
//                        currentPos += length;
//                        currentAltPos+=length;
//                    }
//                    i = j-1;
//                }
//                else  {
//                    //do nothing
//                }
//            }
//        }
//        currentPos = startPos-1;
//        currentAltPos = 0;
//        for (int i = 0; i < md.length(); i++) {
//            char c = md.charAt(i);
//            if (Character.isDigit(c)) {
//                int j;
//                for (j = i + 1; j < md.length(); j++) {
//                    if (!Character.isDigit(md.charAt(j))) {
//                        break;
//                    }
//                }
//                int length = Integer.parseInt(md.substring(i, j));
//                currentPos += length;
//                currentAltPos += length;
//                i = j-1;
//            }
//            else if (Character.isLetter(c)){
//                currentPos++;
//                currentAltPos++;
//                char alt = seq.charAt(currentAltPos-1);
//                snpList.add(new SNPOld(chr, currentPos, c, alt));
//                int index = insertionPosList.binarySearch(currentAltPos);
//                if (index > -2) {
//                    alleleRelaPosList.add((byte)(currentAltPos+frontLength));
//                }
//                else {
//                    int length = 0;
//                    index = -index-1;
//                    for (int j = 0; j < insertionLengthList.size(); j++) {
//                        length+=insertionLengthList.get(j);
//                    }
//                    alleleRelaPosList.add((byte)(currentAltPos+frontLength+length));
//                }
//            }
//            else if (c == '^'){
//                int j;
//                for (j = i + 1; j < md.length(); j++) {
//                    if (!Character.isLetter(md.charAt(j))) break;
//                }
//                String deletionS = md.substring(i+1, j);
//                char alt = 'D';
//                char ref = md.charAt(i+1);
//                snpList.add(new SNPOld(chr, currentPos+1, ref, alt));
//                int index = insertionPosList.binarySearch(currentAltPos);
//                if (index > -2) {
//                    alleleRelaPosList.add((byte)(currentAltPos+frontLength));
//                }
//                else {
//                    int length = 0;
//                    index = -index-1;
//                    for (int k = 0; k < insertionLengthList.size(); k++) {
//                        length+=insertionLengthList.get(k);
//                    }
//                    alleleRelaPosList.add((byte)(currentAltPos+frontLength+length));
//                }
//                currentPos += deletionS.length();
//                i = j -1;
//            }
//            else {
//                System.out.println(c);
//            }
//        }
//        //alleleRelaPosList.sort();
//        if (ifMinus) {
//            int length = seqLength+frontLength+backLength;
//            for (int i = 0; i < alleleRelaPosList.size(); i++) {
//                int value = length-alleleRelaPosList.get(i)+1;
//                alleleRelaPosList.set(i, (byte)value);
//            }
//        }
//        
//        class SNPWithRelativePos  extends SNPOld {
//            byte relaPos = -1;
//            public SNPWithRelativePos(short chr, int pos, byte ref, byte alt) {
//                super(chr, pos, ref, alt);
//            }
//            
//            public SNPWithRelativePos(short chr, int pos, byte ref, TByteArrayList alts, byte relaPos) {
//                super(chr, pos, ref, alts);
//                this.relaPos = relaPos;
//            }
//            
//        }
//        
//        List<SNPWithRelativePos> srList = new ArrayList();
//        for (int i = 0; i < snpList.size(); i++) {
//            SNPWithRelativePos sr = new SNPWithRelativePos (snpList.get(i).getChromosome(), snpList.get(i).getPosition(),
//                snpList.get(i).getRefAlleleByte(), snpList.get(i).getAltAlleleList(), alleleRelaPosList.get(i));
//            srList.add(sr);
//        }
//        Collections.sort(srList);
////*****************************************************        
//        for (int i = sIndex; i < eIndex; i++) {
//            ChrPos query = new ChrPos (chr, sc.getPositionOfSNP(chrIndex, i));
//            //allelePosList.add(query);
//            int index = Collections.binarySearch(srList, query);
//            byte allele;
//            byte relaPos = -1;
//            if (index < 0) {
//                allele = sc.getRefAlleleByteOfSNP(chrIndex, i);
//            }
//            else {
//                allele = srList.get(index).getAltAlleleByte(0);
//                relaPos = srList.get(index).relaPos;
//            }
//            AlleleInfo ai = new AlleleInfo(query.getChromosome(), query.getPosition(), allele, (byte)end, relaPos);
//            tagAlleleList.add(ai);
//        }
//        if (tagAlleleList.size() == 0) return null;
//        return tagAlleleList;
//    }
//    
//    public static Dyad<List<ChrPos>, TByteArrayList> getAlleles (List<String> l, int mapQThresh, SNPCounts sc) {
//        List<ChrPos> allelePosList = new ArrayList();
//        TByteArrayList alleleList = new TByteArrayList();
//        String cigar = l.get(5);
//        if (cigar.startsWith("*")) return null;
//        if (Integer.parseInt(l.get(4)) < mapQThresh) return null;
//        short chr = (short)Integer.parseInt(l.get(2));
//        int startPos = Integer.parseInt(l.get(3));
//        Dyad<TByteArrayList, TIntArrayList> cigarOPPosIndex = getCigarOPAndPosIndex(cigar);
//        int endPos = getEndPos (cigar, cigarOPPosIndex, startPos);
//        int chrIndex = sc.getChrIndex(chr);
//        if (chrIndex < 0) return null;
//        int sIndex = sc.getSNPIndex(chrIndex, startPos);
//        int eIndex = sc.getSNPIndex(chrIndex, endPos);
//        if (sIndex == eIndex) return null;
//        if (sIndex < 0) {
//            sIndex = -sIndex-1;
//            if (sIndex == sc.getSNPNumberOnChromosome(chrIndex)) return null;
//        }
//        if (eIndex < 0) {
//            eIndex = -eIndex-1;
//            if (eIndex == 0) return null;
//        }
//        else {
//            eIndex++;
//        }
////*****copy from getVariants method,  to avoid build cigarOPPosIndex twice******        
//        String seq = l.get(9);
//        String md = l.get(12).split(":")[2];
//        TByteArrayList opList = cigarOPPosIndex.getFirstElement();
//        TIntArrayList posIndexList = cigarOPPosIndex.getSecondElement();
//        if (opList.get(0) == 83) {//S
//            int length = Integer.parseInt(cigar.substring(0, posIndexList.get(0)));
//            seq = seq.substring(length);
//            cigar = cigar.substring(posIndexList.get(0)+1);
//            if (opList.get(opList.size()-1) == 83 || opList.get(opList.size()-1) == 72) {
//                cigarOPPosIndex = getCigarOPAndPosIndex(cigar);
//                opList = cigarOPPosIndex.getFirstElement();
//                posIndexList = cigarOPPosIndex.getSecondElement();
//            }
//        }
//        else if (opList.get(0) == 72) {//H
//            cigar = cigar.substring(posIndexList.get(0)+1);
//            if (opList.get(opList.size()-1) == 83 || opList.get(opList.size()-1) == 72) {
//                cigarOPPosIndex = getCigarOPAndPosIndex(cigar);
//                opList = cigarOPPosIndex.getFirstElement();
//                posIndexList = cigarOPPosIndex.getSecondElement();
//            }
//        }
//        if (opList.get(opList.size()-1) == 83) {//S
//            int length = Integer.parseInt(cigar.substring(posIndexList.get(opList.size()-2)+1, posIndexList.get(opList.size()-1)));
//            seq = seq.substring(0, seq.length()-length);
//            cigar = cigar.substring(0, posIndexList.get(opList.size()-2)+1);
//        }
//        else if (opList.get(opList.size()-1) == 72) {//H
//            cigar = cigar.substring(0, posIndexList.get(opList.size()-2)+1);
//        }
//        List<SNP> snpList = new ArrayList();
//        int currentPos = startPos-1;
//        int currentAltPos = 0;
//        if (opList.contains((byte)73)) {//I 18M1I31M1D46M
//            for (int i = 0; i < cigar.length(); i++) {
//                char c = cigar.charAt(i);
//                if (Character.isDigit(c)) {
//                    int j;
//                    for (j = i + 1; j < cigar.length(); j++) {
//                        if (!Character.isDigit(cigar.charAt(j))) {
//                            break;
//                        }
//                    }
//                    int length = Integer.parseInt(cigar.substring(i, j));
//                    char cop = cigar.charAt(j);
//                    if (cop == 'I') {
//                        char alt = 'I';
//                        char ref = seq.charAt(currentAltPos-1);
//                        snpList.add(new SNPOld(chr, currentPos, ref, alt));
//                        StringBuilder sb = new StringBuilder(seq);
//                        sb.delete(currentAltPos, currentAltPos+length);
//                        seq = sb.toString();
//                    }
//                    else if (cop == 'D' || cop == 'N') {
//                        currentPos += length;
//                    }
//                    else { //M,=,X
//                        currentPos += length;
//                        currentAltPos+=length;
//                    }
//                    i = j-1;
//                }
//                else  {
//                    //do nothing
//                }
//            }
//        }
//        currentPos = startPos-1;
//        currentAltPos = 0;   
//        for (int i = 0; i < md.length(); i++) {
//            char c = md.charAt(i);
//            if (Character.isDigit(c)) {
//                int j;
//                for (j = i + 1; j < md.length(); j++) {
//                    if (!Character.isDigit(md.charAt(j))) {
//                        break;
//                    }
//                }
//                int length = Integer.parseInt(md.substring(i, j));
//                currentPos += length;
//                currentAltPos += length;
//                i = j-1;
//            }
//            else if (Character.isLetter(c)){
//                currentPos++;
//                currentAltPos++;
//                char alt = seq.charAt(currentAltPos-1);
//                snpList.add(new SNPOld(chr, currentPos, c, alt));
//            }
//            else if (c == '^'){
//                int j;
//                for (j = i + 1; j < md.length(); j++) {
//                    if (!Character.isLetter(md.charAt(j))) break;
//                }
//                String deletionS = md.substring(i+1, j);
//                char alt = 'D';
//                char ref = md.charAt(i+1);
//                snpList.add(new SNPOld(chr, currentPos+1, ref, alt));
//                currentPos += deletionS.length();
//                i = j -1;
//            }
//            else {
//                System.out.println(c);
//            }
//        }
//        Collections.sort(snpList);
////*****************************************************        
//        for (int i = sIndex; i < eIndex; i++) {
//            ChrPos query = new ChrPos (chr, sc.getPositionOfSNP(chrIndex, i));
//            allelePosList.add(query);
//            int index = Collections.binarySearch(snpList, query);
//            if (index < 0) alleleList.add(sc.getRefAlleleByteOfSNP(chrIndex, i));
//            else {
//                alleleList.add(snpList.get(index).getAltAlleleByte(0));
//            }
//        }
//        if (allelePosList.size() == 0) return null;
//        return new Dyad(allelePosList, alleleList);
//
//    }
    
    /**
     * Return a list of called SNPs from elements of a SAM alignment record
     * Return null if the seq is not aligned
     * @param mapQThresh
     * @return 
     */
    public static List<SNPOld> getVariants (List<String> l, int mapQThresh) {
        String cigar = l.get(5);
        if (cigar.startsWith("*")) return null;
        if (Integer.parseInt(l.get(4)) < mapQThresh) return null;
        short chr = (short)Integer.parseInt(l.get(2));
        int startPos = Integer.parseInt(l.get(3));
        String seq = l.get(9);
        String md = l.get(12).split(":")[2];
        Dyad<TByteArrayList, TIntArrayList> cigarOPPosIndex = getCigarOPAndPosIndex(cigar);
        TByteArrayList opList = cigarOPPosIndex.getFirstElement();
        TIntArrayList posIndexList = cigarOPPosIndex.getSecondElement();
        if (opList.get(0) == 83) {//S
            int length = Integer.parseInt(cigar.substring(0, posIndexList.get(0)));
            seq = seq.substring(length);
            cigar = cigar.substring(posIndexList.get(0)+1);
            if (opList.get(opList.size()-1) == 83 || opList.get(opList.size()-1) == 72) {
                cigarOPPosIndex = getCigarOPAndPosIndex(cigar);
                opList = cigarOPPosIndex.getFirstElement();
                posIndexList = cigarOPPosIndex.getSecondElement();
            }
        }
        else if (opList.get(0) == 72) {//H
            cigar = cigar.substring(posIndexList.get(0)+1);
            if (opList.get(opList.size()-1) == 83 || opList.get(opList.size()-1) == 72) {
                cigarOPPosIndex = getCigarOPAndPosIndex(cigar);
                opList = cigarOPPosIndex.getFirstElement();
                posIndexList = cigarOPPosIndex.getSecondElement();
            }
        }
        if (opList.get(opList.size()-1) == 83) {//S
            int length = Integer.parseInt(cigar.substring(posIndexList.get(opList.size()-2)+1, posIndexList.get(opList.size()-1)));
            seq = seq.substring(0, seq.length()-length);
            cigar = cigar.substring(0, posIndexList.get(opList.size()-2)+1);
        }
        else if (opList.get(opList.size()-1) == 72) {//H
            cigar = cigar.substring(0, posIndexList.get(opList.size()-2)+1);
        }
        List<SNPOld> snpList = new ArrayList();
        int currentPos = startPos-1;
        int currentAltPos = 0;
        if (opList.contains((byte)73)) {//I 18M1I31M1D46M
            for (int i = 0; i < cigar.length(); i++) {
                char c = cigar.charAt(i);
                if (Character.isDigit(c)) {
                    int j;
                    for (j = i + 1; j < cigar.length(); j++) {
                        if (!Character.isDigit(cigar.charAt(j))) {
                            break;
                        }
                    }
                    int length = Integer.parseInt(cigar.substring(i, j));
                    char cop = cigar.charAt(j);
                    if (cop == 'I') {
                        char alt = 'I';
                        char ref = seq.charAt(currentAltPos-1);
                        snpList.add(new SNPOld(chr, currentPos, ref, alt));
                        StringBuilder sb = new StringBuilder(seq);
                        sb.delete(currentAltPos, currentAltPos+length);
                        seq = sb.toString();
                    }
                    else if (cop == 'D' || cop == 'N') {
                        currentPos += length;
                    }
                    else { //M,=,X
                        currentPos += length;
                        currentAltPos+=length;
                    }
                    i = j-1;
                }
                else  {
                    //do nothing
                }
            }
        }
        currentPos = startPos-1;
        currentAltPos = 0;   
        for (int i = 0; i < md.length(); i++) {
            char c = md.charAt(i);
            if (Character.isDigit(c)) {
                int j;
                for (j = i + 1; j < md.length(); j++) {
                    if (!Character.isDigit(md.charAt(j))) {
                        break;
                    }
                }
                int length = Integer.parseInt(md.substring(i, j));
                currentPos += length;
                currentAltPos += length;
                i = j-1;
            }
            else if (Character.isLetter(c)){
                currentPos++;
                currentAltPos++;
                char alt = seq.charAt(currentAltPos-1);
                snpList.add(new SNPOld(chr, currentPos, c, alt));
            }
            else if (c == '^'){
                int j;
                for (j = i + 1; j < md.length(); j++) {
                    if (!Character.isLetter(md.charAt(j))) break;
                }
                String deletionS = md.substring(i+1, j);
                char alt = 'D';
                char ref = md.charAt(i+1);
                snpList.add(new SNPOld(chr, currentPos+1, ref, alt));
                currentPos += deletionS.length();
                i = j -1;
            }
            else {
                System.out.println(c);
            }
        }
        Collections.sort(snpList);
        return snpList;
    }
}
