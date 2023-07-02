/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.app.grt;

import pgl.infra.dna.BaseEncoder;
import pgl.infra.dna.allele.AlleleEncoder;
import pgl.infra.pos.ChrPos;
import gnu.trove.list.array.TByteArrayList;
import gnu.trove.list.array.TIntArrayList;
import java.io.BufferedWriter;
import java.io.File;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.List;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PArrayUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.Dyad;

/**
 *
 * @author feilu
 */
public class GBSVCFBuilder {
    TagAnnotations tas = null;
    SNPCounts sc = null;
    int numThreads = 32;
    int identityThreshold = 3;
    int maxAltNumber = 2;
    double sequencingAlignErrorRate = 0.05;
    
    public GBSVCFBuilder (TagAnnotations tas, SNPCounts sc) {
        this.tas = tas;
        this.sc = sc;
    }
    
    public void setThreads (int numThreads) {
        this.numThreads = numThreads;
    }
    
    public void setTagIdentifyThreshold (int identityThreshold) {
        this.identityThreshold = identityThreshold;
    }
    
//    public void callGenotypeAllInMemory (String tagBySampleDirS, String genotypeDirS) {
//        File genoDir = new File(genotypeDirS, "genotype");
//        genoDir.mkdir();
//        File[] sampleFiles = new File (tagBySampleDirS).listFiles();
//        sampleFiles = IOUtils.listFilesEndsWith(sampleFiles, ".tas");
//        Arrays.sort(sampleFiles);
//        String[] sampleNames = new String[sampleFiles.length];       
//        for (int i = 0; i < sampleNames.length; i++) {
//            sampleNames[i] = sampleFiles[i].getName().replaceAll(".tas$", "");
//        }
//        int[][] indices = PArrayUtils.getSubsetsIndicesBySubsetSize(sampleFiles.length, this.numThreads);
//        List<File> sampleFileList = Arrays.asList(sampleFiles);
//        TagFinder tf = new TagFinder(tas);
//        System.out.println("\nStart calling genotype of each individual sample...\n");
//        SimpleGenoInfo[][][] taxaGenoInfo = new SimpleGenoInfo[sampleFiles.length][sc.getChromosomeNumber()][];
//        for (int i = 0; i < indices.length; i++) {
//            List<Integer> subIndexList = new ArrayList<>();
//            for (int j = indices[i][0]; j < indices[i][1]; j++) {
//                subIndexList.add(j);
//            }
//            subIndexList.stream().forEach(index -> {
//                AlleleDepth[][] adt = this.initializeADTable();
//                TagAnnotations ata = new TagAnnotations(sampleFileList.get(index).getAbsolutePath());
//                for (int j = 0; j < ata.getGroupNumber(); j++) {
//                    for (int k = 0; k < ata.getTagNumber(j); k++) {
//                        long[] tag = ata.getTag(j, k);
//                        int readDepth = ata.getReadNumber(j, k);
//                        byte r1Length = ata.getR1TagLength(j, k);
//                        byte r2Length = ata.getR1TagLength(j, k);
//                        Dyad<int[], int[]> result = tf.getMostSimilarTags(tag, r1Length, r2Length, j, identityThreshold);
//                        if (result == null) continue;
//                        int[] divergence = result.getFirstElement();
//                        int[] tagIndices = result.getSecondElement();
//                        int tagIndex = this.getTagIndex(divergence, tagIndices, j, tag);
//                        if (tagIndex < 0) continue;
//                        int alleleNumber = tas.getAlleleNumberOfTag(j, tagIndex);
//                        if (alleleNumber == 0) continue;
//                        short chr = tas.getAlleleOfTag(j, tagIndex).get(0).getChromosome();
//                        int chrIndex = sc.getChrIndex(chr);
//                        if (chrIndex < 0) continue;
//                        for (int u = 0; u < alleleNumber; u++) {
//                            AlleleInfo ai = tas.getAlleleOfTag(j, tagIndex).get(u);
//                            int snpIndex = sc.getSNPIndex(chrIndex, new ChrPos(ai.getChromosome(), ai.getPosition()));
//                            adt[chrIndex][snpIndex].addAllele(ai.getAllele());
//                            adt[chrIndex][snpIndex].addDepth(readDepth);
//                        }
//                    }
//                }
//               
//                for (int j = 0; j < adt.length; j++) {
//                    List<SimpleGenoInfo> taxonGenoInfoList = new ArrayList<>();
//                    for (int k = 0; k < adt[j].length; k++) {
//                        adt[j][k].toArray();
//                        if (adt[j][k].getAlleleNumber() < 1) continue;
//                        SimpleGenoInfo sgi = new SimpleGenoInfo ((short)j, k, adt[j][k].getAlleleNumber(), adt[j][k].getAlleles(), adt[j][k].getDepths());
//                        taxonGenoInfoList.add(sgi);
//                    }
//                    SimpleGenoInfo[] taxonGenoInfo = taxonGenoInfoList.toArray(new SimpleGenoInfo[taxonGenoInfoList.size()]);
//                    taxaGenoInfo[index][j] = taxonGenoInfo;
//                }
//            });
//            System.out.println("Completed individual genotyping calling of " + String.valueOf(indices[i][1])+ " samples");
//        }
//        this.writeGenotypeAllInMemory(taxaGenoInfo, sampleNames, genotypeDirS);
//    }
    
//    class SimpleGenoInfo {
//        short chrIndex = Short.MIN_VALUE;
//        int posIndex = Integer.MIN_VALUE;
//        int alleleNumber = 0;
//        byte[] allele = null;
//        int[] depth = null;
//        
//        public SimpleGenoInfo (short chrIndex, int posIndex) {
//            this.chrIndex = chrIndex;
//            this.posIndex = posIndex;
//        }
//        
//        public SimpleGenoInfo(short chrIndex, int posIndex, int alleleNumber, byte[] allele, int[] depth) {
//            this(chrIndex, posIndex);
//            this.alleleNumber = alleleNumber;
//            this.allele = allele;
//            this.depth = depth;
//        }
//        
//        public int getPosIndex () {
//            return posIndex;
//        }
//        
//        public short getChrIndex () {
//            return this.chrIndex;
//        }
//        
//        public byte getAlleleNumber() {
//            return (byte)this.alleleNumber;
//        }
//        
//        public byte getAllele (int index) {
//            return this.allele[index];
//        }
//        
//        public int getDepth (int index) {
//            return this.depth[index];
//        }
//    }
//    
//    private void writeGenotypeAllInMemory (SimpleGenoInfo[][][] taxaGenoInfo, String[] sampleNames, String genotypeDirS) {
//        System.out.println("Start merging individual genotype into VCF by chromosomes");
//        File genoDir = new File(genotypeDirS, "genotype");
//        File[] fs = genoDir.listFiles();
//        for (int i = 0; i < fs.length; i++) fs[i].delete();
//        genoDir.mkdir();
//        int chrNumber = sc.getChromosomeNumber();
//        String[] outfiles = new String[chrNumber];
//        for (int i = 0; i < outfiles.length; i++) {
//            short chr = sc.getChromosome(i);
//            outfiles[i] = new File (genoDir, "chr"+PStringUtils.getNDigitNumber(3, chr)+".vcf").getAbsolutePath();
//        }
//        try {
//            String annotation = VCFUtils.getVCFAnnotation();
//            String header = VCFUtils.getVCFHeader(sampleNames);
//            for (int i = 0; i < sc.getChromosomeNumber(); i++) {
//                BufferedWriter bw = IOUtils.getTextWriter(outfiles[i]);
//                bw.write(annotation);
//                bw.write(header);
//                bw.newLine();
//                int[] currentPosIndices = new int[taxaGenoInfo.length];
//                int[] currentArrayIndices = new int[taxaGenoInfo.length];
//                for (int j = 0; j < taxaGenoInfo.length; j++) {
//                    if (taxaGenoInfo[j][i].length == 0) {
//                        currentPosIndices[j] = Integer.MIN_VALUE;
//                        currentArrayIndices[j] = Integer.MIN_VALUE;
//                    }
//                    else {
//                        currentPosIndices[j] = taxaGenoInfo[j][i][0].getPosIndex();
//                    }
//                }
//                for (int j = 0; j < sc.getSNPNumberOnChromosome(i); j++) {
//                    int[] depth = new int[6];
//                    AlleleDepth[] sampleAD = new AlleleDepth[taxaGenoInfo.length];
//                    for (int k = 0; k < taxaGenoInfo.length; k++) {
//                        sampleAD[k] = new AlleleDepth();
//                        if (taxaGenoInfo[k][i].length == 0) {
//                            sampleAD[k].toArray();
//                            continue;
//                        }
//                        int posIndex = taxaGenoInfo[k][i][currentArrayIndices[k]].getPosIndex();
//                        if (posIndex > j) {
//
//                        }
//                        else if (posIndex < j) {
//                            int taxonArrayIndex = currentArrayIndices[k];
//                            while (posIndex < j) {
//                                taxonArrayIndex++;
//                                if (taxonArrayIndex == taxaGenoInfo[k][i].length) {
//                                    taxonArrayIndex--;
//                                    break;
//                                }
//                                posIndex = taxaGenoInfo[k][i][taxonArrayIndex].getPosIndex();
//                            }
//                            currentPosIndices[k] = posIndex;
//                            currentArrayIndices[k] = taxonArrayIndex;
//                        }
//                        if (posIndex != j) {
//                            sampleAD[k].toArray();
//                            continue;
//                        }
//                        SimpleGenoInfo sgi = taxaGenoInfo[k][i][currentArrayIndices[k]];
//                        short chrIndex = sgi.getChrIndex();
//                        byte alleleNumber = sgi.getAlleleNumber();
//                        for (int u = 0; u < alleleNumber; u++) {
//                            byte cAllele = sgi.getAllele(u);
//                            int cDepth = sgi.getDepth(u);
//                            sampleAD[k].addAllele(cAllele);
//                            sampleAD[k].addDepth(cDepth);
//                            depth[cAllele]+=cDepth;
//                        }
//                        sampleAD[k].toArray();
//                    }
//                    TByteArrayList alleleList = new TByteArrayList();
//                    TIntArrayList depthList = new TIntArrayList();
//                    for (int k = 0; k < depth.length; k++) {
//                        if (depth[k] == 0) continue;
//                        alleleList.add(AlleleEncoder.alleleBytes[k]);
//                        depthList.add(depth[k]);
//                    }
//                    AlleleDepth siteAD = new AlleleDepth(alleleList.toArray(), depthList.toArray());
//                    AlleleDepth altAD = siteAD.getAltAlleleDepth(sc.getRefAlleleByteOfSNP(i, j));
//                    altAD.sortByDepthDesending();
//                    StringBuilder sb = new StringBuilder();
//                    sb.append(sc.getChromosome(i)).append("\t").append(sc.getPositionOfSNP(i, j)).append("\t")
//                            .append(sc.getChromosome(i)).append("-").append(sc.getPositionOfSNP(i, j)).append("\t")
//                            .append(AlleleEncoder.alleleByteCharMap.get(sc.getRefAlleleByteOfSNP(i, j))).append("\t");
//                    int altNum = altAD.getAlleleNumber();
//                    if (altNum > this.maxAltNumber) altNum = this.maxAltNumber;
//                    for (int k = 0; k < altNum; k++) {
//                        sb.append(AlleleEncoder.alleleByteCharMap.get(altAD.getAllele(k))).append(",");
//                    }
//                    sb.deleteCharAt(sb.length()-1).append("\t.\t.\t");
//                    sb.append("DP=").append(VCFUtils.getTotalDepth(sampleAD)).append(";AD=").append(VCFUtils.getAlleleTotalDepth(sampleAD, sc.getRefAlleleByteOfSNP(i, j))).append(",");
//                    for (int k = 0; k < altNum; k++) {
//                        sb.append(VCFUtils.getAlleleTotalDepth(sampleAD, altAD.getAllele(k))).append(",");
//                    }
//                    sb.deleteCharAt(sb.length()-1).append(";NS=").append(VCFUtils.getNumberOfTaxaWithAlleles(sampleAD)).append(";AP=").append(VCFUtils.getNumberOfTaxaWithAllele(sampleAD, sc.getRefAlleleByteOfSNP(i, j))).append(",");
//                    for (int k = 0; k < altNum; k++) {
//                        sb.append(VCFUtils.getNumberOfTaxaWithAllele(sampleAD, altAD.getAllele(k))).append(",");
//                    }
//                    sb.deleteCharAt(sb.length()-1);
//                    sb.append("\tGT:AD:PL");
//                    for (int k = 0; k < sampleAD.length; k++) {
//                        int n = altNum+1;
//                        int[] readCount = new int[n];
//                        readCount[0] = sampleAD[k].getDepth(sc.getRefAlleleByteOfSNP(i, j));
//                        for (int u = 0; u < altNum; u++) {
//                            readCount[u+1] = sampleAD[k].getDepth(altAD.getAllele(u));
//                        }
//                        sb.append("\t").append(VCFUtils.getGenotype(readCount, sequencingAlignErrorRate));
//                    }
//                    bw.write(sb.toString());
//                    bw.newLine();   
//                }
//                bw.flush();
//                bw.close();
//                System.gc();
//            }
//        }
//        catch (Exception e) {
//            e.printStackTrace();
//        }
//        System.out.println("VCF genotype is output to " + genotypeDirS);
//    }
    
    public void callGenotype (String tagBySampleDirS, String genotypeDirS) {
        File tempDir = new File(genotypeDirS, "temp");
        tempDir.mkdir();
        File genoDir = new File(genotypeDirS, "genotype");
        genoDir.mkdir();
        File[] sampleFiles = new File (tagBySampleDirS).listFiles();
        sampleFiles = IOUtils.listFilesEndsWith(sampleFiles, ".tas");
        Arrays.sort(sampleFiles);
        String[] sampleNames = new String[sampleFiles.length];       
        for (int i = 0; i < sampleNames.length; i++) {
            sampleNames[i] = sampleFiles[i].getName().replaceAll(".tas$", "");
        }
        int[][] indices = PArrayUtils.getSubsetsIndicesBySubsetSize(sampleFiles.length, this.numThreads);
        List<File> sampleFileList = Arrays.asList(sampleFiles);
        TagFinder tf = new TagFinder(tas);
        System.out.println("\nStart calling genotype of each individual sample...\n");
        for (int i = 0; i < indices.length; i++) {
            List<File> subFList = sampleFileList.subList(indices[i][0], indices[i][1]);
            subFList.parallelStream().forEach(f -> {
                String tempFileS = new File(tempDir, f.getName().replaceAll(".tas", ".gen")).getAbsolutePath();
                AlleleDepth[][] adt = this.initializeADTable();
                TagAnnotations ata = new TagAnnotations(f.getAbsolutePath());
                for (int j = 0; j < ata.getGroupNumber(); j++) {
                    for (int k = 0; k < ata.getTagNumber(j); k++) {
                        long[] tag = ata.getTag(j, k);
                        int readDepth = ata.getReadNumber(j, k);
                        byte r1Length = ata.getR1TagLength(j, k);
                        byte r2Length = ata.getR1TagLength(j, k);
                        Dyad<int[], int[]> result = tf.getMostSimilarTags(tag, r1Length, r2Length, j, identityThreshold);
                        if (result == null) continue;
                        int[] divergence = result.getFirstElement();
                        int[] tagIndices = result.getSecondElement();
                        int tagIndex = this.getTagIndex(divergence, tagIndices, j, tag);
                        if (tagIndex < 0) continue;
                        int alleleNumber = tas.getAlleleNumberOfTag(j, tagIndex);
                        if (alleleNumber == 0) continue;
                        short chr = tas.getAlleleOfTag(j, tagIndex).get(0).getChromosome();
                        int chrIndex = sc.getChrIndex(chr);
                        if (chrIndex < 0) continue;
                        for (int u = 0; u < alleleNumber; u++) {
                            AlleleInfo ai = tas.getAlleleOfTag(j, tagIndex).get(u);
                            int snpIndex = sc.getSNPIndex(chrIndex, new ChrPos(ai.getChromosome(), ai.getPosition()));
                            adt[chrIndex][snpIndex].addAllele(ai.getAllele());
                            adt[chrIndex][snpIndex].addDepth(readDepth);
                        }
                    }
                }
                for (int j = 0; j < adt.length; j++) {
                    for (int k = 0; k < adt[j].length; k++) {
                        adt[j][k].toArray();
                    }
                }
                this.writeTempGenotype(tempFileS, adt);
            });
        }
        if (sampleNames.length > 512) {
            this.writeGenotypeViaMergedFile(tempDir, sampleNames, genotypeDirS);
        }
        else {
            this.writeGenotype(tempDir, sampleNames, genotypeDirS);
        }
//        
//        File[] tempfs = tempDir.listFiles();
//        for (int i = 0; i < tempfs.length; i++) tempfs[i].delete();
//        tempDir.delete();
    }
    
    private void writeGenotypeViaMergedFile (File tempDir, String[] sampleNames, String genotypeDirS) {
        System.out.println("Merging individual genotypes");
        File mergeF = new File(genotypeDirS, "mergeGeno.bin");
        Arrays.sort(sampleNames);
        long[] fileSizes = new long[sampleNames.length];
        try {
            Dyad<FileChannel, ByteBuffer> wt = IOUtils.getNIOChannelBufferWriter(mergeF.getAbsolutePath(), 65536);
            FileChannel wf = wt.getFirstElement();
            for (int i = 0; i < sampleNames.length; i++) {
                String inputFileS = new File(tempDir, sampleNames[i]+".gen").getAbsolutePath();
                fileSizes[i] = new File(inputFileS).length();
                Dyad<FileChannel, ByteBuffer> rt = IOUtils.getNIOChannelBufferReader(inputFileS, 65536);
                FileChannel rf = rt.getFirstElement();
                ByteBuffer rb = rt.getSecondElement();
                while ((rf.read(rb)) != -1) {
                    rb.flip();
                    wf.write(rb);
                    rb.clear();
                }
                rf.close();
            }
            wf.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("Start merging individual genotype into VCF by chromosomes");
        File genoDir = new File(genotypeDirS, "genotype");
        File[] fs = genoDir.listFiles();
        for (int i = 0; i < fs.length; i++) fs[i].delete();
        genoDir.mkdir();
        int chrNumber = sc.getChromosomeNumber();
        String[] outfiles = new String[chrNumber];
        for (int i = 0; i < outfiles.length; i++) {
            short chr = sc.getChromosome(i);
            outfiles[i] = new File (genoDir, "chr"+PStringUtils.getNDigitNumber(3, chr)+".vcf").getAbsolutePath();
        }
        int ii = 0;
        int jj = 0;
        int kk = 0;
        try {
            FileChannel fc = FileChannel.open(Paths.get(mergeF.getAbsolutePath()), StandardOpenOption.READ); 
            ByteBuffer[] bbs = new ByteBuffer[sampleNames.length];
            long[] pointers = new long[fileSizes.length]; 
            for (int i = 0; i < bbs.length; i++) {
                bbs[i] = ByteBuffer.allocate(65536);
                if (i == 0) continue;
                pointers[i] = pointers[i-1]+fileSizes[i-1];
            }
            for (int i = 0; i < bbs.length; i++) {
                fc.read(bbs[i], pointers[i]);
                bbs[i].flip();
            }
            String annotation = GRTVCFUtils.getVCFAnnotation();
            String header = GRTVCFUtils.getVCFHeader(sampleNames);
            int totalCnt = 0;
            int dbSNPCnt = 0;
            for (int i = 0; i < sc.getChromosomeNumber(); i++) {
                BufferedWriter bw = IOUtils.getTextWriter(outfiles[i]);
                bw.write(annotation);
                bw.write(header);
                bw.newLine();
                int cnt = 0;
                dbSNPCnt+=sc.getSNPNumberOnChromosome(i);
                for (int j = 0; j < sc.getSNPNumberOnChromosome(i); j++) {
                    int[] depth = new int[6];
                    AlleleDepth[] sampleAD = new AlleleDepth[bbs.length];
                    for (int k = 0; k < bbs.length; k++) {
                        ii = i;jj = j;kk = k;
                        short chrIndex = bbs[k].getShort();
                        if (chrIndex == Short.MIN_VALUE) {
                            bbs[k].clear();
                            fc.read(bbs[k]);
                            bbs[k].flip();
                            chrIndex = bbs[k].getShort();
                        }
                        else if (chrIndex == -1) {
                            
                        }
                        int posIndex = bbs[k].getInt();
                        if (chrIndex != i || posIndex != j) {
                            bbs[k].position(bbs[k].position()-6);
                            sampleAD[k] = new AlleleDepth();
                        }
                        else {
                            sampleAD[k] = new AlleleDepth();
                            byte alleleNumber = bbs[k].get();
                            for (int u = 0; u < alleleNumber; u++) {
                                byte cAllele = bbs[k].get();
                                int cDepth = bbs[k].getInt();
                                sampleAD[k].addAllele(cAllele);
                                sampleAD[k].addDepth(cDepth);
                                depth[cAllele]+=cDepth;
                            }
                        }
                        sampleAD[k].toArray();
                    }
                    TByteArrayList alleleList = new TByteArrayList();
                    TIntArrayList depthList = new TIntArrayList();
                    for (int k = 0; k < depth.length; k++) {
                        if (depth[k] == 0) continue;
                        alleleList.add(AlleleEncoder.alleleCodings[k]);
                        depthList.add(depth[k]);
                    }
                    AlleleDepth siteAD = new AlleleDepth(alleleList.toArray(), depthList.toArray());
                    AlleleDepth altAD = siteAD.getAltAlleleDepth(sc.getRefAlleleByteOfSNP(i, j));
                    int altNum = altAD.getAlleleNumber();
                    if (altNum == 0) continue;
                    altAD.sortByDepthDesending();
                    StringBuilder sb = new StringBuilder();
                    sb.append(sc.getChromosome(i)).append("\t").append(sc.getPositionOfSNP(i, j)).append("\t")
                            .append(sc.getChromosome(i)).append("-").append(sc.getPositionOfSNP(i, j)).append("\t")
                            .append(AlleleEncoder.alleleCodingToBaseMap.get(sc.getRefAlleleByteOfSNP(i, j))).append("\t");
                    if (altNum > this.maxAltNumber) altNum = this.maxAltNumber;
                    for (int k = 0; k < altNum; k++) {
                        sb.append(AlleleEncoder.alleleCodingToBaseMap.get(altAD.getAllele(k))).append(",");
                    }
                    sb.deleteCharAt(sb.length()-1).append("\t.\t.\t");
                    sb.append("DP=").append(GRTVCFUtils.getTotalDepth(sampleAD)).append(";AD=").append(GRTVCFUtils.getAlleleTotalDepth(sampleAD, sc.getRefAlleleByteOfSNP(i, j))).append(",");
                    for (int k = 0; k < altNum; k++) {
                        sb.append(GRTVCFUtils.getAlleleTotalDepth(sampleAD, altAD.getAllele(k))).append(",");
                    }
                    sb.deleteCharAt(sb.length()-1).append(";NS=").append(GRTVCFUtils.getNumberOfTaxaWithAlleles(sampleAD)).append(";AP=").append(GRTVCFUtils.getNumberOfTaxaWithAllele(sampleAD, sc.getRefAlleleByteOfSNP(i, j))).append(",");
                    for (int k = 0; k < altNum; k++) {
                        sb.append(GRTVCFUtils.getNumberOfTaxaWithAllele(sampleAD, altAD.getAllele(k))).append(",");
                    }
                    sb.deleteCharAt(sb.length()-1);
                    sb.append("\tGT:AD:PL");
                    for (int k = 0; k < sampleAD.length; k++) {
                        int n = altNum+1;
                        int[] readCount = new int[n];
                        readCount[0] = sampleAD[k].getDepth(sc.getRefAlleleByteOfSNP(i, j));
                        for (int u = 0; u < altNum; u++) {
                            readCount[u+1] = sampleAD[k].getDepth(altAD.getAllele(u));
                        }
                        sb.append("\t").append(GRTVCFUtils.getGenotype(readCount, sequencingAlignErrorRate));
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                    cnt++;
                }
                totalCnt+=cnt;
                bw.flush();
                bw.close();
                System.out.println("Genotyping on chromosome "+String.valueOf(sc.getChromosome(i))+" with "+String.valueOf(cnt)+ " SNPs finished");
                System.gc();
            }
            fc.close();
            DecimalFormat df = new DecimalFormat("##.##%");
            double percent = ((double)totalCnt/ dbSNPCnt);
            String formattedPercent = df.format(percent);
            System.out.println("A total of " + String.valueOf(dbSNPCnt) + " SNPs are in the database");
            System.out.println("A total of " + String.valueOf(totalCnt) + ", or " +formattedPercent+ " SNPs are segregating and genotyped");
        }
        catch (Exception e) {
            System.out.println(ii+"\t"+jj+"\t"+kk);
            e.printStackTrace();
        }
        mergeF.delete();
        System.out.println("VCF genotype is output to " + genotypeDirS);
    }
    
    private void writeGenotype (File tempDir, String[] sampleNames, String genotypeDirS) {
        System.out.println("Start merging individual genotype into VCF by chromosomes");
        File genoDir = new File(genotypeDirS, "genotype");
        File[] fs = genoDir.listFiles();
        for (int i = 0; i < fs.length; i++) fs[i].delete();
        genoDir.mkdir();
        Arrays.sort(sampleNames);
        int chrNumber = sc.getChromosomeNumber();
        String[] outfiles = new String[chrNumber];
        for (int i = 0; i < outfiles.length; i++) {
            short chr = sc.getChromosome(i);
            outfiles[i] = new File (genoDir, "chr"+PStringUtils.getNDigitNumber(3, chr)+".vcf").getAbsolutePath();
        }
        int ii = 0;
        int jj = 0;
        int kk = 0;
        try {
            FileChannel[] fcs = new FileChannel[sampleNames.length]; 
            ByteBuffer[] bbs = new ByteBuffer[fcs.length];
            for (int i = 0; i < fcs.length; i++) {
                String inputFileS = new File(tempDir, sampleNames[i]+".gen").getAbsolutePath();
                Dyad<FileChannel, ByteBuffer> iot = IOUtils.getNIOChannelBufferReader(inputFileS, 65536);
                fcs[i] = iot.getFirstElement();
                bbs[i] = iot.getSecondElement();
            }
            for (int i = 0; i < fcs.length; i++) {
                fcs[i].read(bbs[i]);
                bbs[i].flip();
            }
            String annotation = GRTVCFUtils.getVCFAnnotation();
            String header = GRTVCFUtils.getVCFHeader(sampleNames);
            int totalCnt = 0;
            int dbSNPCnt = 0;
            for (int i = 0; i < sc.getChromosomeNumber(); i++) {
                BufferedWriter bw = IOUtils.getTextWriter(outfiles[i]);
                bw.write(annotation);
                bw.write(header);
                bw.newLine();
                int cnt = 0;
                dbSNPCnt+=sc.getSNPNumberOnChromosome(i);
                for (int j = 0; j < sc.getSNPNumberOnChromosome(i); j++) {
                    int[] depth = new int[6];
                    AlleleDepth[] sampleAD = new AlleleDepth[fcs.length];
                    for (int k = 0; k < fcs.length; k++) {
                        ii = i;jj = j;kk = k;
                        short chrIndex = bbs[k].getShort();
                        if (chrIndex == Short.MIN_VALUE) {
                            bbs[k].clear();
                            fcs[k].read(bbs[k]);
                            bbs[k].flip();
                            chrIndex = bbs[k].getShort();
                        }
                        else if (chrIndex == -1) {
                            
                        }
                        int posIndex = bbs[k].getInt();
                        if (chrIndex != i || posIndex != j) {
                            bbs[k].position(bbs[k].position()-6);
                            sampleAD[k] = new AlleleDepth();
                        }
                        else {
                            sampleAD[k] = new AlleleDepth();
                            byte alleleNumber = bbs[k].get();
                            for (int u = 0; u < alleleNumber; u++) {
                                byte cAllele = bbs[k].get();
                                int cDepth = bbs[k].getInt();
                                sampleAD[k].addAllele(cAllele);
                                sampleAD[k].addDepth(cDepth);
                                depth[cAllele]+=cDepth;
                            }
                        }
                        sampleAD[k].toArray();
                    }
                    TByteArrayList alleleList = new TByteArrayList();
                    TIntArrayList depthList = new TIntArrayList();
                    for (int k = 0; k < depth.length; k++) {
                        if (depth[k] == 0) continue;
                        alleleList.add(AlleleEncoder.alleleCodings[k]);
                        depthList.add(depth[k]);
                    }
                    AlleleDepth siteAD = new AlleleDepth(alleleList.toArray(), depthList.toArray());
                    AlleleDepth altAD = siteAD.getAltAlleleDepth(sc.getRefAlleleByteOfSNP(i, j));
                    int altNum = altAD.getAlleleNumber();
                    if (altNum == 0) continue;
                    altAD.sortByDepthDesending();
                    StringBuilder sb = new StringBuilder();
                    sb.append(sc.getChromosome(i)).append("\t").append(sc.getPositionOfSNP(i, j)).append("\t")
                            .append(sc.getChromosome(i)).append("-").append(sc.getPositionOfSNP(i, j)).append("\t")
                            .append(AlleleEncoder.alleleCodingToBaseMap.get(sc.getRefAlleleByteOfSNP(i, j))).append("\t");
                    if (altNum > this.maxAltNumber) altNum = this.maxAltNumber;
                    for (int k = 0; k < altNum; k++) {
                        sb.append(AlleleEncoder.alleleCodingToBaseMap.get(altAD.getAllele(k))).append(",");
                    }
                    sb.deleteCharAt(sb.length()-1).append("\t.\t.\t");
                    sb.append("DP=").append(GRTVCFUtils.getTotalDepth(sampleAD)).append(";AD=").append(GRTVCFUtils.getAlleleTotalDepth(sampleAD, sc.getRefAlleleByteOfSNP(i, j))).append(",");
                    for (int k = 0; k < altNum; k++) {
                        sb.append(GRTVCFUtils.getAlleleTotalDepth(sampleAD, altAD.getAllele(k))).append(",");
                    }
                    sb.deleteCharAt(sb.length()-1).append(";NS=").append(GRTVCFUtils.getNumberOfTaxaWithAlleles(sampleAD)).append(";AP=").append(GRTVCFUtils.getNumberOfTaxaWithAllele(sampleAD, sc.getRefAlleleByteOfSNP(i, j))).append(",");
                    for (int k = 0; k < altNum; k++) {
                        sb.append(GRTVCFUtils.getNumberOfTaxaWithAllele(sampleAD, altAD.getAllele(k))).append(",");
                    }
                    sb.deleteCharAt(sb.length()-1);
                    sb.append("\tGT:AD:PL");
                    for (int k = 0; k < sampleAD.length; k++) {
                        int n = altNum+1;
                        int[] readCount = new int[n];
                        readCount[0] = sampleAD[k].getDepth(sc.getRefAlleleByteOfSNP(i, j));
                        for (int u = 0; u < altNum; u++) {
                            readCount[u+1] = sampleAD[k].getDepth(altAD.getAllele(u));
                        }
                        sb.append("\t").append(GRTVCFUtils.getGenotype(readCount, sequencingAlignErrorRate));
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                    cnt++;
                }
                totalCnt+=cnt;
                bw.flush();
                bw.close();
                System.out.println("Genotyping on chromosome "+String.valueOf(sc.getChromosome(i))+" with "+String.valueOf(cnt)+ " SNPs finished");
                System.gc();
            }
            for (int i = 0; i < fcs.length; i++) {
                fcs[i].close();
            }
            DecimalFormat df = new DecimalFormat("##.##%");
            double percent = ((double)totalCnt/ dbSNPCnt);
            String formattedPercent = df.format(percent);
            System.out.println("A total of " + String.valueOf(dbSNPCnt) + " SNPs are in the database");
            System.out.println("A total of " + String.valueOf(totalCnt) + ", or " +formattedPercent+ " SNPs are segregating and genotyped");
        }
        catch (Exception e) {
            System.out.println(ii+"\t"+jj+"\t"+kk);
            e.printStackTrace();
        }
        System.out.println("VCF genotype is output to " + genotypeDirS);
    }
    
    private void wirteTextTempGenotype (String tempFileS, AlleleDepth[][] adt) {
        try {
            BufferedWriter bw = IOUtils.getTextWriter(tempFileS);
            bw.write("ChrIndex\tSNPIndex\tAlleleIndex\tAllele\tDepth");
            bw.newLine();
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < adt.length; i++) {
                for (int j = 0; j < adt[i].length; j++) {
                    int alleleNumber  = adt[i][j].getAlleleNumber();
                    for (int k = 0; k < alleleNumber; k++) {
                        sb = new StringBuilder();
                        sb.append(i).append("\t").append(j).append("\t").append(k).append("\t").append(adt[i][j].getAllele(k)).append("\t").append(adt[i][j].getDepth(k));
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                }
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private void writeTempGenotype (String tempFileS, AlleleDepth[][] adt) {
        try {
            Dyad<FileChannel, ByteBuffer> iot = IOUtils.getNIOChannelBufferWriter(tempFileS, 65536);
            FileChannel fc = iot.getFirstElement();
            ByteBuffer bb = iot.getSecondElement();
            for (int i = 0; i < adt.length; i++) {
                for (int j = 0; j < adt[i].length; j++) {
                    int alleleNumber = adt[i][j].getAlleleNumber();
                    if (alleleNumber == 0) continue;
                    bb.putShort((short)i);
                    bb.putInt(j);
                    bb.put((byte)alleleNumber);
                    for (int k = 0; k < alleleNumber; k++) {
                        bb.put(adt[i][j].getAllele(k));
                        bb.putInt(adt[i][j].getDepth(k));
                    }
                    int remain = bb.remaining();
                    if (remain < 50) {
                        bb.putShort(Short.MIN_VALUE);
                        bb.putInt(Integer.MIN_VALUE);
                        bb.flip();
                        bb.limit(bb.capacity());
                        fc.write(bb);
                        bb.clear();
                    }
                }
            }
            if (bb.hasRemaining()) {
                bb.putShort((short)-1);
                bb.putInt(Integer.MIN_VALUE);
                bb.flip();
                bb.limit(bb.capacity());
                fc.write(bb);
                bb.clear();
            }
            fc.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private AlleleDepth[][] initializeADTable () {
        AlleleDepth[][] adt = new AlleleDepth[sc.getChromosomeNumber()][];
        for (int i = 0; i < adt.length; i++) {
            adt[i] = new AlleleDepth[sc.getSNPNumberOnChromosome(i)];
            for (int j = 0; j < adt[i].length; j++) {
                adt[i][j] = new AlleleDepth();
            }
        }
        return adt;
    }
    
    private int getTagIndex (int[] divergence, int[] tagIndices, int groupIndex, long[] queryTag) {
        byte[][] query = new byte[2][];
        for (int i = 0; i < 2; i++) {
            query[i] = BaseEncoder.getBaseCodingArrayFromLongs(Arrays.copyOfRange(queryTag, i*tas.getTagLengthInLong(), (i+1)*tas.getTagLengthInLong()));
            //System.out.println(BaseEncoder.getSequenceFromLongs(Arrays.copyOfRange(queryTag, i*tas.getTagLengthInLong(), (i+1)*tas.getTagLengthInLong())));
        }
        for (int i = 0; i < tagIndices.length; i++) {
            int cnt = 0;
            List<AlleleInfo> ai = tas.getAlleleOfTag(groupIndex, tagIndices[i]);
            byte[][] dbTag = new byte[2][];
            for (int j = 0; j < dbTag.length; j++) {
                dbTag[j] = BaseEncoder.getBaseCodingArrayFromLongs(Arrays.copyOfRange(tas.getTag(groupIndex, tagIndices[i]), j*tas.getTagLengthInLong(), (j+1)*tas.getTagLengthInLong()));
                //System.out.println(BaseEncoder.getSequenceFromLongs(Arrays.copyOfRange(tas.getTag(groupIndex, tagIndices[i]), j*tas.getTagLengthInLong(), (j+1)*tas.getTagLengthInLong())));
            }
            for (int j = 0; j < ai.size(); j++) {
                byte end = ai.get(j).getEnd();
                byte base = ai.get(j).getBase();
                byte pos = ai.get(j).getRelativePosition();
                if (query[end-1][pos-1] == base) cnt++;
            }
            if (cnt == ai.size()) return tagIndices[i];
        }
        return Integer.MIN_VALUE;
    }
}

