/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.app.grt;

import static cern.jet.math.Arithmetic.factorial;
import com.koloboke.collect.map.hash.HashIntDoubleMap;
import com.koloboke.collect.map.hash.HashIntDoubleMaps;
import pgl.app.grt.AlleleDepth;

import java.text.SimpleDateFormat;
import java.util.Date;

/**
 * incomplete
 * @author feilu
 */
public class GRTVCFUtils {
    private static final int maxFactorial = 150;
    private static final HashIntDoubleMap factorialMap = createFactorialMap(maxFactorial);
    public static String getVCFHeader (String[] taxaNames) {
        StringBuilder sb = new StringBuilder();
        sb.append("#CHR\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
        for (int i = 0; i < taxaNames.length; i++) {
            sb.append("\t").append(taxaNames[i]);
        }
        return sb.toString();
    }
    
    public static String getVCFAnnotation () {
        StringBuilder sb = new StringBuilder();
        sb.append("##fileformat=VCFv4.2\n");
        SimpleDateFormat sdf = new SimpleDateFormat("MM/dd/yyyy HH:mm:ss.SSS");
        Date dt = new Date();
        String S = sdf.format(dt);
        sb.append("##fileDate=").append(S.split("\\.")[0]).append("\n");
        sb.append("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"").append("Combined depth across samples").append("\">\n");
        sb.append("##INFO=<ID=AD,Number=2+,Type=Integer,Description=\"").append("Allele depth across samples").append("\">\n");
        sb.append("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"").append("Number of samples with data").append("\">\n");
        sb.append("##INFO=<ID=AP,Number=2+,Type=Integer,Description=\"").append("Number of samples in which an allele is present").append("\">\n");
        sb.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"").append("Genotype").append("\">\n");
        sb.append("##FORMAT=<ID=AD,Number=1,Type=String,Description=\"").append("Allele depth of a sample").append("\">\n");
        sb.append("##FORMAT=<ID=PL,Number=1,Type=String,Description=\"").append("Genotype likelihoods").append("\">\n");
        return sb.toString();
    }
    
    private static HashIntDoubleMap createFactorialMap (int maxFactorial) {
        HashIntDoubleMap theMap = HashIntDoubleMaps.getDefaultFactory().withDefaultValue(-1).newMutableMap();
        for (int i = 0; i < maxFactorial+1; i++) {
            theMap.put(i, factorial(i));
        }
        return theMap;
    }
    
    
    public static int getNumberOfTaxaWithAllele (AlleleDepth[] ad, byte allele) {
        int cnt = 0;
        for (int i = 0; i < ad.length; i++) {
            if (ad[i].binarySearch(allele) < 0) continue;
            cnt++;
        }
        return cnt;
    }
    
    public static int getAlleleTotalDepth (AlleleDepth[] ad, byte allele) {
        int cnt = 0;
        for (int i = 0; i < ad.length; i++) {
            cnt+=ad[i].getDepth(allele);
        }
        return cnt;
    }
    
    public static int getNumberOfTaxaWithAlleles (AlleleDepth[] ad) {
        int cnt = 0;
        for (int i = 0; i < ad.length; i++) {
            if (ad[i].getAlleleNumber() == 0) continue;
            cnt++;
        }
        return cnt;
    }
    
    public static int getTotalDepth (AlleleDepth[] ad) {
        int cnt = 0;
        for (int i = 0; i < ad.length; i++) {
            cnt+=ad[i].getTotalDepth();
        }
        return cnt;
    }
    
    /**
     * 
     * @param alleleReadCount read count of alleles, alleles should be in Ref-Alt order
     * @param sequencingErrorRate
     * @return 
     */
    public static String getGenotype (int[] alleleReadCount, double sequencingErrorRate) {
        int n = alleleReadCount.length*(alleleReadCount.length+1)/2;
        int[] likelihood = new int[n];
        int sum = 0;
        for (int i = 0; i < alleleReadCount.length; i++) sum+=alleleReadCount[i];
        if (sum == 0) return "./.";
        else if (sum > maxFactorial) {
            double portion = (double)maxFactorial/sum;
            for (int i = 0; i < alleleReadCount.length; i++) {
                alleleReadCount[i] = (int)(alleleReadCount[i]*portion);
            }
            sum = maxFactorial;
        }
        double coe = factorialMap.get(sum);
        for (int i = 0; i < alleleReadCount.length; i++) coe = coe/factorialMap.get(alleleReadCount[i]);
        double max = Double.MAX_VALUE;
        int a1 = 0;
        int a2 = 0;
        for (int i = 0; i < alleleReadCount.length; i++) {
            for (int j = i; j < alleleReadCount.length; j++) {
                int index = (j*(j+1)/2)+i;
                double value = Double.MAX_VALUE;
                if (i == j) {
                    value = -Math.log10(coe*Math.pow((1-0.75*sequencingErrorRate), alleleReadCount[i])*Math.pow(sequencingErrorRate/4, (sum-alleleReadCount[i])));
                }
                else {
                    value = -Math.log10(coe*Math.pow((0.5-sequencingErrorRate/4), alleleReadCount[i]+alleleReadCount[j])*Math.pow(sequencingErrorRate/4, (sum-alleleReadCount[i]-alleleReadCount[j])));
                }
                if (value < max) {
                    max = value;
                    a1 = i;
                    a2 = j;
                }
                likelihood[index] = (int)Math.round(value);
            }
        }
        StringBuilder sb = new StringBuilder();
        sb.append(a1).append("/").append(a2).append(":");
        for (int i = 0; i < alleleReadCount.length; i++) sb.append(alleleReadCount[i]).append(",");
        sb.deleteCharAt(sb.length()-1); sb.append(":");
        for (int i = 0; i < likelihood.length; i++) sb.append(likelihood[i]).append(",");
        sb.deleteCharAt(sb.length()-1);
        return sb.toString();
    }
}
