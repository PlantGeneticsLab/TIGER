/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.dna.genot;

import com.koloboke.collect.map.hash.HashIntDoubleMap;
import com.koloboke.collect.map.hash.HashIntDoubleMaps;
import pgl.app.grt.AlleleDepth;
import java.text.SimpleDateFormat;
import java.util.Date;
import static cern.jet.math.Arithmetic.factorial;

/**
 * incomplete
 * @author feilu
 */
public class VCFUtils {
    private static final int maxFactorial = 150;
    private static final HashIntDoubleMap factorialMap = createFactorialMap(maxFactorial);
    
    private static HashIntDoubleMap createFactorialMap (int maxFactorial) {
        HashIntDoubleMap theMap = HashIntDoubleMaps.getDefaultFactory().withDefaultValue(-1).newMutableMap();
        for (int i = 0; i < maxFactorial+1; i++) {
            theMap.put(i, factorial(i));
        }
        return theMap;
    }
    
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

    /**
     * Return genotype of VCF via max likelihood of multinomial sampling distribution
     * See https://doi.org/10.1371/journal.pgen.1000862
     * @param cnt allele depth of Ref and Alt alleles
     * @param combinedErrorRate sequencing error rate and alignment error rate. set to 0.05 in general
     * @return
     */
    private String getGenotypeByInt (short[] cnt, double combinedErrorRate) {
        //in case some allele depth is greater than maxFactorial, to keep the original allele counts
        short[] oriCnt = null;
        int n = cnt.length*(cnt.length+1)/2;
        int[] likelihood = new int[n];
        int sum = 0;
        for (int i = 0; i < cnt.length; i++) sum+=cnt[i];
        if (sum == 0) return "./.";
        else if (sum > this.maxFactorial) {
            oriCnt = new short[cnt.length];
            System.arraycopy(cnt, 0, oriCnt, 0, cnt.length);
            double portion = (double)this.maxFactorial/sum;
            for (int i = 0; i < cnt.length; i++) {
                cnt[i] = (short)(cnt[i]*portion);
            }
            sum = this.maxFactorial;
        }
        double coe = this.factorialMap.get(sum);
        for (int i = 0; i < cnt.length; i++) coe = coe/this.factorialMap.get(cnt[i]);
        double max = Double.MAX_VALUE;
        int a1 = 0;
        int a2 = 0;
        for (int i = 0; i < cnt.length; i++) {
            for (int j = i; j < cnt.length; j++) {
                int index = (j*(j+1)/2)+i;
                double value = Double.MAX_VALUE;
                if (i == j) {
                    value = -Math.log10(coe*Math.pow((1-0.75*combinedErrorRate), cnt[i])*Math.pow(combinedErrorRate /4, (sum-cnt[i])));
                }
                else {
                    value = -Math.log10(coe*Math.pow((0.5-combinedErrorRate /4), cnt[i]+cnt[j])*Math.pow(combinedErrorRate /4, (sum-cnt[i]-cnt[j])));
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
        if (sum > this.maxFactorial) {
            for (int i = 0; i < oriCnt.length; i++) sb.append(oriCnt[i]).append(",");
        }
        else {
            for (int i = 0; i < cnt.length; i++) sb.append(cnt[i]).append(",");
        }
        sb.deleteCharAt(sb.length()-1); sb.append(":");
        for (int i = 0; i < likelihood.length; i++) sb.append(likelihood[i]).append(",");
        sb.deleteCharAt(sb.length()-1);
        return sb.toString();
    }
}
