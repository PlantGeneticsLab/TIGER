/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.app.hapScanner;
import static cern.jet.math.Arithmetic.factorial;
import com.koloboke.collect.map.IntDoubleMap;
import com.koloboke.collect.map.hash.HashIntDoubleMaps;
import pgl.infra.table.RowTable;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.atomic.LongAdder;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PArrayUtils;
import pgl.infra.utils.PStringUtils;

/**
 *
 * @author feilu
 */
public class HapScanner {
    //The taxaRefBam file containing information of taxon and its corresponding bam files. The bam file should have .bai file in the same folder
    String taxaRefBamFileS = null;
    //The posAllele file (with header), the format is Chr\tPos\tRef\tAlt (from VCF format). The positions come from haplotype library.
    String posAlleleFileS = null;
    //The pos files (without header), the format is Chr\tPos. The positions come from haplotype library.
    String posFileS = null;
    //The chromosome which will be genotyped
    int chr = -1;
    //The path of samtools
    String samtoolsPath = null;
    //The directory of output
    String outputDirS = null;
    
    int nThreads = -1;
    
    HashMap<String, List<String>> taxaBamsMap = new HashMap<>();
    
    HashMap<String, String> taxaRefMap = new HashMap<>();
    
    String[] subDirS = {"mpileup", "indiVCF", "VCF"};
    
    IntDoubleMap factorialMap = null;
    int maxFactorial = 150;
    //combined: sequencing error and alignment error
    double combinedErrorRate = 0.05;
    
    public HapScanner (String infileS) {
        this.parseParameters(infileS);
        this.mkDir();
        this.scanIndiVCF();
        this.mkFinalVCF();
    }

    /**
     * Multiple threading for merging VCF, but slow
     * @deprecated
     */
    public void mkFinalVCF2 () {
        Set<String> taxaSet = taxaBamsMap.keySet();
        String[] taxa = taxaSet.toArray(new String[taxaSet.size()]);
        Arrays.sort(taxa);
        List<String> taxaList = Arrays.asList(taxa);
        String outfileS = new File(outputDirS, subDirS[2]).getAbsolutePath();
        outfileS = new File(outfileS, "chr"+PStringUtils.getNDigitNumber(3, chr)+".vcf.gz").getAbsolutePath();
        int blockSize = 10000;
        int actualSize = 100000;
        String[][] indiVCF = new String[blockSize][taxa.length];
        String[] blockVCF = new String[blockSize];
        List<Integer> blockIndexList = new ArrayList<>();
        for (int i = 0; i < blockSize; i++) {
            blockIndexList.add(i);
        }
        try {
            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
            bw.write("##fileformat=VCFv4.1\n");
            bw.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
            bw.write("##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the reference and alternate alleles in the order listed\">\n");
            bw.write("##FORMAT=<ID=GL,Number=.,Type=Integer,Description=\"Genotype likelihoods for 0/0, 0/1, 1/1, or  0/0, 0/1, 0/2, 1/1, 1/2, 2/2 if 2 alt alleles\">\n");
            bw.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n");
            bw.write("##INFO=<ID=NZ,Number=1,Type=Integer,Description=\"Number of taxa with called genotypes\">\n");
            bw.write("##INFO=<ID=AD,Number=.,Type=Integer,Description=\"Total allelelic depths in order listed starting with REF\">\n");
            bw.write("##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Numbers of ALT alleles in order listed\">\n");
            bw.write("##INFO=<ID=GN,Number=.,Type=Integer,Description=\"Number of taxa with genotypes AA,AB,BB or AA,AB,AC,BB,BC,CC if 2 alt alleles\">\n");
            bw.write("##INFO=<ID=HT,Number=1,Type=Integer,Description=\"Number of heterozygotes\">\n");
            bw.write("##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Minor allele frequency\">\n");
            bw.write("##ALT=<ID=DEL,Description=\"Deletion\">\n");
            bw.write("##ALT=<ID=INS,Description=\"Insertion\">\n");
            bw.write("##HapMapVersion=\"3.2.1\"\n");
            StringBuilder sb = new StringBuilder("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT");
            for (int i = 0; i < taxa.length; i++) {
                sb.append("\t").append(taxa[i]);
            }
            bw.write(sb.toString());
            bw.newLine();
            
            HashMap<String, Integer> taxaIndexMap = new HashMap<>();
            BufferedReader[] brs = new BufferedReader[taxa.length];
            String indiVCFFolderS = new File(outputDirS, subDirS[1]).getAbsolutePath();
            for (int i = 0; i < brs.length; i++) {
                String indiVCFFileS = new File (indiVCFFolderS, taxa[i]+".indi.vcf").getAbsolutePath();
                brs[i] = IOUtils.getTextReader(indiVCFFileS);
                taxaIndexMap.put(taxa[i], i);
            }
            BufferedReader br = IOUtils.getTextGzipReader(posAlleleFileS);
            List<String> temList = null;
            String temp = br.readLine();
            int cntBlock = 0;
            int cnt = 0;
            StringBuilder[] sbs = new StringBuilder[blockSize];
            String[] alts = new String[blockSize];
            for (int i = 0; i < sbs.length; i++) {
                sbs[i] = new StringBuilder();
            }
            while ((temp = br.readLine()) != null) {
                sbs[cntBlock].setLength(0);
                temList = PStringUtils.fastSplit(temp);
                sbs[cntBlock].append(temList.get(0)).append("\t").append(temList.get(1)).append("\t").append(temList.get(0)).append("-").append(temList.get(1)).append("\t");
                sbs[cntBlock].append(temList.get(2)).append("\t").append(temList.get(3)).append("\t.\t.\t");
                alts[cntBlock] = temList.get(3);
                cntBlock++;
                if (cntBlock%blockSize == 0) {
                    actualSize = blockSize;
                    this.updateBlockVCF(actualSize, brs, taxaList, blockVCF, taxaIndexMap, indiVCF, alts, blockIndexList);
                    for (int i = 0; i < actualSize; i++) {
                        sbs[i].append(blockVCF[i]);
                        bw.write(sbs[i].toString());
                        bw.newLine();
                    }
                    cntBlock = 0;
                    for (int i = 0; i < sbs.length; i++) {
                        sbs[i].setLength(0);
                    }
                }
                cnt++;
                if (cnt%100000 == 0) System.out.println(String.valueOf(cnt)+" SNPs output to " + outfileS);
            }
            if (cntBlock != 0) {
                actualSize = cntBlock;
                this.updateBlockVCF(actualSize, brs, taxaList, blockVCF, taxaIndexMap, indiVCF, alts, blockIndexList);
                for (int i = 0; i < actualSize; i++) {
                    sbs[i].append(blockVCF[i]);
                    bw.write(sbs[i].toString());
                    bw.newLine();
                }
                for (int i = 0; i < sbs.length; i++) {
                    sbs[i].setLength(0);
                }
            }
            bw.flush();
            bw.close();
            br.close();
            for (int i = 0; i < brs.length; i++) {
                brs[i].close();
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("Final VCF is completed at " + outfileS);
    }
    
    public void updateBlockVCF (int actualSize, BufferedReader[] brs, List<String> taxaList, String[] blockVCF, HashMap<String, Integer> taxaIndexMap, String[][] indiVCF, String[] alts, List<Integer> blockIndexList) {
        taxaList.parallelStream().forEach(taxon -> {
            int taxonIndex = taxaIndexMap.get(taxon);
            try {
                for (int i = 0; i < actualSize; i++) {
                    indiVCF[i][taxonIndex] = brs[taxonIndex].readLine();
                }
            }
            catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
        blockIndexList.parallelStream().forEach(i -> {
            StringBuilder sb = new StringBuilder();
            sb.append(this.getInfo(indiVCF[i], alts[i])).append("\tGT:AD:GL");
            for (int j = 0; j < indiVCF[i].length; j++) {
                sb.append("\t").append(indiVCF[i][j]);
            }
            blockVCF[i] = sb.toString();
        });
    }
    
    public void mkFinalVCF () {
        Set<String> taxaSet = taxaBamsMap.keySet();
        String[] taxa = taxaSet.toArray(new String[taxaSet.size()]);
        Arrays.sort(taxa);
        String outfileS = new File(outputDirS, subDirS[2]).getAbsolutePath();
        outfileS = new File(outfileS, "chr"+PStringUtils.getNDigitNumber(3, chr)+".vcf").getAbsolutePath();
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("##fileformat=VCFv4.1\n");
            bw.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
            bw.write("##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the reference and alternate alleles in the order listed\">\n");
            bw.write("##FORMAT=<ID=GL,Number=.,Type=Integer,Description=\"Genotype likelihoods for 0/0, 0/1, 1/1, or  0/0, 0/1, 0/2, 1/1, 1/2, 2/2 if 2 alt alleles\">\n");
            bw.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n");
            bw.write("##INFO=<ID=NZ,Number=1,Type=Integer,Description=\"Number of taxa with called genotypes\">\n");
            bw.write("##INFO=<ID=AD,Number=.,Type=Integer,Description=\"Total allelelic depths in order listed starting with REF\">\n");
            bw.write("##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Numbers of ALT alleles in order listed\">\n");
            bw.write("##INFO=<ID=GN,Number=.,Type=Integer,Description=\"Number of taxa with genotypes AA,AB,BB or AA,AB,AC,BB,BC,CC if 2 alt alleles\">\n");
            bw.write("##INFO=<ID=HT,Number=1,Type=Integer,Description=\"Number of heterozygotes\">\n");
            bw.write("##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Minor allele frequency\">\n");
            bw.write("##ALT=<ID=DEL,Description=\"Deletion\">\n");
            bw.write("##ALT=<ID=INS,Description=\"Insertion\">\n");
            StringBuilder sb = new StringBuilder("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT");
            for (int i = 0; i < taxa.length; i++) {
                sb.append("\t").append(taxa[i]);
            }
            bw.write(sb.toString());
            bw.newLine();
            BufferedReader[] brs = new BufferedReader[taxa.length];
            String indiVCFFolderS = new File(outputDirS, subDirS[1]).getAbsolutePath();
            for (int i = 0; i < brs.length; i++) {
                String indiVCFFileS = new File (indiVCFFolderS, taxa[i]+".chr"+PStringUtils.getNDigitNumber(3, chr)+".indi.vcf").getAbsolutePath();
                brs[i] = new BufferedReader (new FileReader(indiVCFFileS), 4096);
            }
            BufferedReader br = null;
            if (posAlleleFileS.endsWith(".gz")) {
                br = IOUtils.getTextGzipReader(posAlleleFileS);
            }
            else {
                br = IOUtils.getTextReader(posAlleleFileS);
            }
            
            String temp = br.readLine();
            List<String> temList = null;
            int cnt = 0;
            String[] genoArray = new String[brs.length];
            while ((temp = br.readLine()) != null) {
                sb.setLength(0);
                temList = PStringUtils.fastSplit(temp);
                sb.append(temList.get(0)).append("\t").append(temList.get(1)).append("\t").append(temList.get(0)).append("-").append(temList.get(1)).append("\t");
                sb.append(temList.get(2)).append("\t").append(temList.get(3)).append("\t.\t.\t");
                for (int i = 0; i < brs.length; i++) {
                    genoArray[i]= brs[i].readLine();
                }
                sb.append(this.getInfo(genoArray, temList.get(3))).append("\tGT:AD:GL");
                for (int i = 0; i < genoArray.length; i++) {
                    sb.append("\t").append(genoArray[i]);
                }
                bw.write(sb.toString());
                bw.newLine();
                cnt++;
                if (cnt%1000000 == 0) System.out.println(String.valueOf(cnt)+" SNPs output to " + outfileS);
            }
            bw.flush();
            bw.close();
            br.close();
            for (int i = 0; i < brs.length; i++) {
                brs[i].close();
            }
            for (int i = 0; i < taxa.length; i++) {
                String indiVCFFileS = new File (indiVCFFolderS, taxa[i]+".chr"+PStringUtils.getNDigitNumber(3, chr)+".indi.vcf").getAbsolutePath();
                new File(indiVCFFileS).delete();
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
        System.out.println("Final VCF is completed at " + outfileS);
    }
    
    private String getInfo (String[] genoArray, String altList) {
        int dp = 0;
        int nz = 0;
        int nAlt = PStringUtils.fastSplit(altList, ",").size();
        int[] adCnt = new int[1+nAlt];
        int[] acCnt = new int[1+nAlt];
        int[][] gnCnt = new int[1+nAlt][1+nAlt];
        int ht = 0;
        List<String> tempList = null;
        List<String> temList = null;
        for (int i = 0; i < genoArray.length; i++) {
            if (genoArray[i].startsWith(".")) {
                nz++;
                continue;
            }
            tempList = PStringUtils.fastSplit(genoArray[i], ":");
            temList = PStringUtils.fastSplit(tempList.get(1), ",");
            for (int j = 0; j < temList.size(); j++) {
                int c = Integer.parseInt(temList.get(j));
                dp+=c;
                adCnt[j] += c;
            }
            temList = PStringUtils.fastSplit(tempList.get(0), "/");
            for (int j = 0; j < temList.size(); j++) {
                int c = Integer.parseInt(temList.get(j));
                acCnt[c]++;
            }
            int index1 = Integer.parseInt(temList.get(0));
            int index2 = Integer.parseInt(temList.get(1));
            gnCnt[index1][index2]++;
            if (index1 != index2) ht++;
        }
        nz = genoArray.length - nz;
        int sum = 0;
        for (int i = 0; i < acCnt.length; i++) {
            sum+=acCnt[i];
        }
        float maf = (float)((double)acCnt[0]/sum);
        if (maf>0.5) maf = (float)(1-maf);
        StringBuilder sb = new StringBuilder();
        sb.append("DP=").append(dp).append(";NZ=").append(nz).append(";AD=");
        for (int i = 0; i < adCnt.length; i++) {
            sb.append(adCnt[i]).append(",");
        }
        sb.deleteCharAt(sb.length()-1);
        sb.append(";AC=");
        for (int i = 1; i < acCnt.length; i++) {
            sb.append(acCnt[i]).append(",");
        }
        sb.deleteCharAt(sb.length()-1);
        sb.append(";GN=");
        for (int i = 0; i < gnCnt.length; i++) {
            for (int j = i + 1; j < gnCnt.length; j++) {
                sb.append(gnCnt[i][j]).append(",");
            }
        }
        sb.deleteCharAt(sb.length()-1);
        sb.append(";HT=").append(ht).append(";MAF=").append(maf);
        return sb.toString();
    }
    
    public void scanIndiVCF () {
        this.creatFactorialMap();
        RowTable<String> t = new RowTable<> (posAlleleFileS);
        HashMap<Integer, String> posRefMap = new HashMap<>();
        HashMap<Integer, String[]> posAltMap = new HashMap<>();
        int[] positions = new int[t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            positions[i] = t.getCellAsInteger(i, 1);
            posRefMap.put(positions[i], t.getCell(i, 2));
            String[] tem = t.getCell(i, 3).split(",");
            posAltMap.put(t.getCellAsInteger(i, 1), tem);
        }
        Set<String> taxaSet = taxaBamsMap.keySet();
        ArrayList<String> taxaList = new ArrayList(taxaSet);
        Collections.sort(taxaList);
        
        int[][] indices = PArrayUtils.getSubsetsIndicesBySubsetSize(taxaList.size(), this.nThreads);
        LongAdder counter = new LongAdder();
        for (int u = 0; u < indices.length; u++) {
            List<String> subTaxaList = taxaList.subList(indices[u][0], indices[u][1]);
            subTaxaList.parallelStream().forEach(e -> {
                String pileupFolderS = new File(outputDirS, subDirS[0]).getAbsolutePath();
                String pileupFileS = new File (pileupFolderS, e+".chr"+PStringUtils.getNDigitNumber(3, chr)+".pileup.txt").getAbsolutePath();
                String indiVCFFolderS = new File(outputDirS, subDirS[1]).getAbsolutePath();
                String indiVCFFileS = new File (indiVCFFolderS, e+".chr"+PStringUtils.getNDigitNumber(3, chr)+".indi.vcf").getAbsolutePath();
                List<String> bamPaths = taxaBamsMap.get(e);
                StringBuilder sb = new StringBuilder(samtoolsPath);
                sb.append(" mpileup -A -B -q 20 -Q 20 -f ").append(this.taxaRefMap.get(e));
                for (int i = 0; i < bamPaths.size(); i++) {
                    sb.append(" ").append(bamPaths.get(i));
                }
                sb.append(" -l ").append(posFileS).append(" -r ");
                sb.append(chr).append(" -o ").append(pileupFileS);
                String command = sb.toString();
                //System.out.println(command);
                try {
                    Runtime rt = Runtime.getRuntime();
                    Process p = rt.exec(command);
                    p.waitFor();
                }
                catch (Exception ee) {
                    ee.printStackTrace();
                }
                counter.increment();
                int cnt = counter.intValue();
                if (cnt%10 == 0) System.out.println("Pileuped " + String.valueOf(cnt) + " taxa. Total: " + String.valueOf(taxaBamsMap.size()));
                try {
                    File pf = new File (pileupFileS);
                    if (!pf.exists()) {
                        System.out.println("Pileup file does not exist. "+ pf.getName());
                        System.out.println(command);
                    }
                    BufferedReader br = IOUtils.getTextReader(pileupFileS);
                    BufferedWriter bw = IOUtils.getTextWriter(indiVCFFileS);
                    String current = br.readLine();
                    List<String> currentList = null;
                    int currentPosition = -1;
                    if (current != null) {
                        currentList = PStringUtils.fastSplit(current);
                        currentPosition = Integer.parseInt(currentList.get(1));
                    }
                    for (int i = 0; i < positions.length; i++) {
                        if (current == null) {
                            bw.write("./.");
                            bw.newLine();
                        }
                        else {
                            if (positions[i] == currentPosition) {
                                String ref = posRefMap.get(currentPosition);
                                String[] alts = posAltMap.get(currentPosition);
                                char[] alleleC = new char[alts.length+1];
                                alleleC[0] = ref.charAt(0);
                                for (int j = 0; j < alts.length; j++) {
                                    if (alts[j].startsWith("I") || alts[j].startsWith("<I")) {
                                        alleleC[j+1] = '+';
                                    }
                                    else if (alts[j].startsWith("D") || alts[j].startsWith("<D")) {
                                        alleleC[j+1] = '-';
                                    }
                                    alleleC[j+1] = alts[j].charAt(0);
                                }
                                int[] cnts = new int[alts.length+1];
                                sb = new StringBuilder();
                                for (int j = 0; j < bamPaths.size(); j++) {
                                    sb.append(currentList.get(4+j*3));
                                }
                                StringBuilder ssb = new StringBuilder();
                                int curIndex = 0;
                                for (int j = 0; j < sb.length(); j++) {
                                    char cChar = sb.charAt(j);
                                    if (cChar == '+') {
                                        ssb.append(sb.subSequence(curIndex, j+1));
                                        curIndex = j+2+Character.getNumericValue(sb.charAt(j+1));
                                    }
                                    else if (cChar == '-') {
                                        ssb.append(sb.subSequence(curIndex, j+1));
                                        curIndex = j+2+Character.getNumericValue(sb.charAt(j+1));
                                    }
                                }
                                ssb.append(sb.subSequence(curIndex, sb.length()));
                                sb = ssb;
                                String s = sb.toString().toUpperCase();
                                for (int j = 0; j < s.length(); j++) {
                                    char cChar = s.charAt(j);
                                    if (cChar == '.' || cChar == ',') cnts[0]++;
                                    for (int k = 1; k < alleleC.length; k++) {
                                        if (cChar == alleleC[k]) cnts[k]++;
                                    }
                                }
                                for (int j = 1; j < alleleC.length; j++) {
                                    if (alleleC[j] == '+') cnts[0] = cnts[0]-cnts[j];
                                    else if (alleleC[j] == '-') cnts[0] = cnts[0]-cnts[j];
                                }
                                String vcf = this.getGenotype(cnts);
                                bw.write(vcf);
                                bw.newLine();

                                current = br.readLine();
                                if (current != null) {
                                    currentList = PStringUtils.fastSplit(current);
                                    currentPosition = Integer.parseInt(currentList.get(1));
                                } 
                            }
                            else if (positions[i] < currentPosition) {
                                bw.write("./.");
                                bw.newLine();    
                            }
                            else {
                                System.out.println("Currentposition is greater than pileup position. It should not happen. Program quits");
                                System.exit(1);
                            }
                        }
                    }
                    bw.flush();
                    bw.close();
                    br.close();
                    new File(pileupFileS).delete();
                }
                catch (Exception ee) {
                    ee.printStackTrace();
                    System.exit(1);
                }
            });
        }        
       
        
    }
    
    private String getGenotype (int[] cnt) {
        int n = cnt.length*(cnt.length+1)/2;
        int[] likelihood = new int[n];
        int sum = 0;
        for (int i = 0; i < cnt.length; i++) sum+=cnt[i];
        if (sum == 0) return "./.";
        else if (sum > this.maxFactorial) {
            double portion = (double)this.maxFactorial/sum;
            for (int i = 0; i < cnt.length; i++) {
                cnt[i] = (int)(cnt[i]*portion);
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
                    value = -Math.log10(coe*Math.pow((1-0.75*this.combinedErrorRate), cnt[i])*Math.pow(this.combinedErrorRate /4, (sum-cnt[i])));
                }
                else {
                    value = -Math.log10(coe*Math.pow((0.5-this.combinedErrorRate /4), cnt[i]+cnt[j])*Math.pow(this.combinedErrorRate /4, (sum-cnt[i]-cnt[j])));
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
        for (int i = 0; i < cnt.length; i++) sb.append(cnt[i]).append(",");
        sb.deleteCharAt(sb.length()-1); sb.append(":");
        for (int i = 0; i < likelihood.length; i++) sb.append(likelihood[i]).append(",");
        sb.deleteCharAt(sb.length()-1);
        return sb.toString();
    }
    
    private void creatFactorialMap () {
        this.factorialMap = HashIntDoubleMaps.getDefaultFactory().newMutableMap();
        for (int i = 0; i < this.maxFactorial+1; i++) {
            this.factorialMap.put(i, factorial(i));
        }
    }
    
    public void mkDir () {
        for (int i = 0; i < subDirS.length; i++) {
            File f = new File(outputDirS, subDirS[i]);
            f.mkdir();
        }
    }
    
    public void parseParameters (String infileS) {
        List<String> pLineList = new ArrayList<>();
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            String temp = null;
            boolean ifOut = false;
            if (!(temp = br.readLine()).equals("HapScanner")) ifOut = true;
            if (!(temp = br.readLine()).equals("Author: Aoyue Bi, Xuebo Zhao, Fei Lu")) ifOut = true;
            if (!(temp = br.readLine()).equals("Email: biaoyue17@genetics.ac.cn; xuebozhao@genetics.ac.cn; flu@genetics.ac.cn")) ifOut = true;
            if (!(temp = br.readLine()).equals("Homepage: http://plantgeneticslab.weebly.com/")) ifOut = true;
            if (ifOut) {
                System.out.println("Thanks for using HapScanner.");
                System.out.println("Please keep the authorship in the parameter file. Program stops.");
                System.exit(0);
            }
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) continue;
                if (temp.isEmpty()) continue;
                pLineList.add(temp);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        taxaRefBamFileS = pLineList.get(0);
        posAlleleFileS = pLineList.get(1);
        posFileS = pLineList.get(2);
        chr = Integer.valueOf(pLineList.get(3));
        samtoolsPath = pLineList.get(4);
        this.nThreads = Integer.parseInt(pLineList.get(5));
        outputDirS = pLineList.get(6);
        this.combinedErrorRate = Double.parseDouble(pLineList.get(7));
        new File(outputDirS).mkdir();
        RowTable<String> t = new RowTable<>(taxaRefBamFileS);
        Set<String> taxaSet = new HashSet<>();
        for (int i = 0; i < t.getRowNumber(); i++) {
            taxaSet.add(t.getCell(i, 0));
        }
        String[] taxaArray = taxaSet.toArray(new String[taxaSet.size()]);
        Arrays.sort(taxaArray);
        for (int i = 0; i < taxaArray.length; i++) {
            List<String> bamList = new ArrayList<>();
            taxaBamsMap.put(taxaArray[i], bamList);
        }
        for (int i = 0; i < t.getRowNumber(); i++) {
            String key = t.getCell(i, 0);
            List<String> bamList = taxaBamsMap.get(key);
            bamList.add(t.getCell(i, 2));
            taxaBamsMap.replace(key, bamList);
            taxaRefMap.put(key, t.getCell(i, 1));
        }
    }
    
    public static void main (String[] args) {
        new HapScanner(args[0]);
    }
}
