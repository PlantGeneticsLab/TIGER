/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.dna.allele;

import com.koloboke.collect.map.hash.HashByteCharMap;
import com.koloboke.collect.map.hash.HashByteCharMaps;
import com.koloboke.collect.map.hash.HashCharByteMap;
import com.koloboke.collect.map.hash.HashCharByteMaps;

/**
 * Class encoding alleles and site genotype from A, C, G, T, D, I
 * @author feilu
 */
public class AlleleEncoder {

    /**
     * Alleles in char, D is deletion, I is insertion, . is missing
     */
    public static final char[] alleleBases = {'A', 'C', 'G', 'T', 'D', 'I'};
    /**
     * Alleles in AscII code
     */
    public static final byte[] alleleAscIIs = {65, 67, 71, 84, 68, 73};
    /**
     * Alleles in byte code
     */
    public static final byte[] alleleBytes = {0, 1, 2, 3, 4, 5};

    /**
     * The default value of allele missing
     */
    public static final byte alleleMissingByte = 15;

    /**
     * The default char of allele missing
     */
    public static final char alleleMissingBase = 'N';

    /**
     * The default value of genotype missing
     */
    public static final byte genotypeMissingByte = -1;

    /**
     * Converter from char to allele byte
     */
    public static final HashCharByteMap alleleBaseToByteMap =
            HashCharByteMaps.getDefaultFactory().withDefaultValue(alleleMissingByte).newImmutableMap(alleleBases, alleleBytes);

    /**
     * Converter from allele byte to char
     */
    public static final HashByteCharMap alleleByteToBaseMap =
            HashByteCharMaps.getDefaultFactory().withDefaultValue(alleleMissingBase).newImmutableMap(alleleBytes, alleleBases);

    /**
     * Return an allele byte from char
     * @param c
     * @return -1 if not found
     */
    public static byte getAlleleByteFromBase(char c) {
        return alleleBaseToByteMap.get(c);
    }

    /**
     * Return an allele char from byte
     * @param b
     * @return ! if not found
     */
    public static char getAlleleBaseFromByte(byte b) {
        return alleleByteToBaseMap.get(b);
    }

    /**
     * Return a genotype byte from allele bytes of two homologous chromosomes
     * @param b1
     * @param b2
     * @return
     */
    public static byte getGenotypeByte (byte b1, byte b2) {
        return (byte)((b1<<4)+b2);
    }

    /**
     * Return a genotype byte from allele chars of two homologous chromosomes
     * @param c1
     * @param c2
     * @return
     */
    public static byte getGenotypeByte (char c1, char c2) {
        return getGenotypeByte(getAlleleByteFromBase(c1), getAlleleByteFromBase(c2));
    }

    /**
     * Return the allele byte of the 1st homologous chromosome
     * @param g
     * @return
     */
    public static byte getAlleleByte1FromGenotypeByte (byte g) {
        return (byte)(g>>>4);
    }

    /**
     * Return the allele byte of the 2nd homologous chromosome
     * @param g
     * @return
     */
    public static byte getAlleleByte2FromGenotypeByte (byte g) {
        return (byte)(g&15);
    }

    /**
     * Return the allele char of the 1st homologous chromosome
     * @param g
     * @return
     */
    public static char getAlleleBase1FromGenotypeByte(byte g) {
        return getAlleleBaseFromByte(getAlleleByte1FromGenotypeByte(g));
    }

    /**
     * Return the allele char of the 2nd homologous chromosome
     * @param g
     * @return
     */
    public static char getAlleleBase2FromGenotypeByte(byte g) {
        return getAlleleBaseFromByte(getAlleleByte2FromGenotypeByte(g));
    }

    /**
     * Return a genotype string from its genotype byte, e.g. AA or TG
     * @param g
     * @return
     */
    public static String getGenotypeStringFromGenotypeByte (byte g) {
        StringBuilder sb = new StringBuilder();
        sb.append(getAlleleBase1FromGenotypeByte(g));
        sb.append(getAlleleBase2FromGenotypeByte(g));
        return sb.toString();
    }
}
