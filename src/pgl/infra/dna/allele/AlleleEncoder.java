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
 * Class encoding alleles and site genotype from A, C, G, T, D, I.
 * <p>
 * Binary encoding rules: A = 0001, C = 0010, G = 0011, T = 0100, D = 0101, I = 0110, AlleleMissing = N = 1111.
 * @author feilu
 */
public class AlleleEncoder {

    /**
     * Alleles in char, D is deletion, I is insertion, "." is missing.
     */
    public static final char[] alleleBases = {'A', 'C', 'G', 'T', 'D', 'I'};

    /**
     * Alleles in AscII code.
     */
    public static final byte[] alleleAscIIs = {65, 67, 71, 84, 68, 73};

    /**
     * Alleles in coding values.
     */
    public static final byte[] alleleCodings = {0, 1, 2, 3, 4, 5};

    /**
     * The default coding value of allele missing.
     */
    public static final byte alleleMissingCoding = 15;

    /**
     * The default char of allele missing.
     */
    public static final char alleleMissingBase = 'N';

    /**
     * The default coding value of genotype missing.
     */
    public static final byte genotypeMissingCoding = -1;

    /**
     * A pre-defined builder to accelerate building genotype string.
     */
    private static StringBuilder genoSb = new StringBuilder();

    /**
     * Converter from char to allele coding.
     */
    public static final HashCharByteMap alleleBaseToCodingMap =
            HashCharByteMaps.getDefaultFactory().withDefaultValue(alleleMissingCoding).newImmutableMap(alleleBases, alleleCodings);

    /**
     * Converter from allele coding to char.
     */
    public static final HashByteCharMap alleleCodingToBaseMap =
            HashByteCharMaps.getDefaultFactory().withDefaultValue(alleleMissingBase).newImmutableMap(alleleCodings, alleleBases);

    /**
     * Return an allele coding from char, return -1 if the char is "non-ACGTDI".
     * @param c
     * @return -1 if not found
     */
    public static byte getAlleleCodingFromBase(char c) {
        return alleleBaseToCodingMap.get(c);
    }

    /**
     * Return an allele char from coding, return "N" if the coding is not 0, 1, 2, 3, 4, 5.
     * @param b
     * @return "N" if not found
     */
    public static char getAlleleBaseFromCoding(byte b) {
        return alleleCodingToBaseMap.get(b);
    }

    /**
     * Return a genotype coding from allele codings of two homologous chromosomes.
     * @param b1
     * @param b2
     * @return
     */
    public static byte getGenotypeCoding(byte b1, byte b2) {
        return (byte)((b1<<4)+b2);
    }

    /**
     * Return a genotype coding from allele chars of two homologous chromosomes.
     * @param c1
     * @param c2
     * @return
     */
    public static byte getGenotypeCoding(char c1, char c2) {
        return getGenotypeCoding(getAlleleCodingFromBase(c1), getAlleleCodingFromBase(c2));
    }

    /**
     * Return the allele coding of the 1st homologous chromosome.
     * @param g
     * @return
     */
    public static byte getAlleleCoding1FromGenotypeCoding(byte g) {
        return (byte)(g>>>4);
    }

    /**
     * Return the allele coding of the 2nd homologous chromosome.
     * @param g
     * @return
     */
    public static byte getAlleleCoding2FromGenotypeByte(byte g) {
        return (byte)(g&15);
    }

    /**
     * Return the allele char of the 1st homologous chromosome.
     * @param g
     * @return
     */
    public static char getAlleleBase1FromGenotypeCoding(byte g) {
        return getAlleleBaseFromCoding(getAlleleCoding1FromGenotypeCoding(g));
    }

    /**
     * Return the allele char of the 2nd homologous chromosome.
     * @param g
     * @return
     */
    public static char getAlleleBase2FromGenotypeCoding(byte g) {
        return getAlleleBaseFromCoding(getAlleleCoding2FromGenotypeByte(g));
    }

    /**
     * Return a genotype string from its genotype coding, for example, AA or TG
     * @param g
     * @return
     */
    public static String getGenotypeStringFromGenotypeCoding(byte g) {
        genoSb.setLength(0);
        genoSb.append(getAlleleBase1FromGenotypeCoding(g));
        genoSb.append(getAlleleBase2FromGenotypeCoding(g));
        return genoSb.toString();
    }
}
