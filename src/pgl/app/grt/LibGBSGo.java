/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.app.grt;

import pgl.infra.pos.ChrPos;

import java.io.BufferedReader;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

import pgl.infra.utils.IOUtils;

/**
 *
 * @author feilu
 */
public class LibGBSGo {
    String workingDirS = null;
    String barcodeFileS = null;
    String libraryFastqMapFileS = null;
    String referenceFileS = null;
    String bwaPath = null;
    String cutter1 = null;
    String cutter2 = null;
    String[] subDirS = {"tagsBySample","tagsLibrary","alignment", "rawGenotype", "filteredGenotype", "queryGenotype"};
//    String[] subDirS = {"tagsBySample_large","tagsLibrary_large","alignment_large", "rawGenotype", "filteredGenotype", "queryGenotype"};
    LibraryInfo li = null;
    
    public LibGBSGo (String parameterFileS) {
        this.initializeParameter(parameterFileS);
//        this.mkTagsBySample();
//        this.mergeTagAnnotations();
//        this.alignTags();
//        this.callSNP();
//        this.removeLowCountSNP();
//        this.callAllele();
//        this.buildVCF();
        //this.filterDatabase();
        //this.retrieveGenotype();
    }
    
    public void retrieveGenotype () {
        String tagBySampleDirS = new File (this.workingDirS, this.subDirS[0]).getAbsolutePath();
        String tagLibraryDirS = new File (this.workingDirS, this.subDirS[1]).getAbsolutePath();
        String filteredTagAnnotationFileS = new File(tagLibraryDirS, "tag.filtered.tas").getAbsolutePath();
        String filteredSNPFileS = new File(tagLibraryDirS, "SNP.filtered.bin").getAbsolutePath();
        String genotypeDirS = new File (this.workingDirS, this.subDirS[5]).getAbsolutePath();
        TagAnnotations tas = new TagAnnotations(filteredTagAnnotationFileS);
        SNPCounts sc = new SNPCounts (filteredSNPFileS);
        GBSVCFBuilder builder = new GBSVCFBuilder(tas, sc);
        builder.setTagIdentifyThreshold(3);
        builder.callGenotype(tagBySampleDirS, genotypeDirS);
    }
    
    public void filterDatabase () {
        String tagLibraryDirS = new File (this.workingDirS, this.subDirS[1]).getAbsolutePath();
        String tagAnnotationFileS = new File(tagLibraryDirS, "tag.tas").getAbsolutePath();
        String rawSNPFileS = new File(tagLibraryDirS, "rawSNP.bin").getAbsolutePath();
        String filteredGenotypeDirS = new File (this.workingDirS, this.subDirS[4]).getAbsolutePath();
        String filteredTagAnnotationFileS = new File(tagLibraryDirS, "tag.filtered.tas").getAbsolutePath();
        String filteredSNPFileS = new File(tagLibraryDirS, "SNP.filtered.bin").getAbsolutePath();
        TagAnnotations tas = new TagAnnotations(tagAnnotationFileS);
        SNPCounts sc = new SNPCounts (rawSNPFileS);
        List<ChrPos> validatedSNPPosList = tas.filterTagAnnotationsWithValidatedGenotype(filteredGenotypeDirS);
        tas.writeBinaryFile(filteredTagAnnotationFileS);
        sc.selectSNPs(validatedSNPPosList);
        sc.writeBinaryFile(filteredSNPFileS);
    }
    
    public void buildVCF () {
        String tagBySampleDirS = new File (this.workingDirS, this.subDirS[0]).getAbsolutePath();
        String tagLibraryDirS = new File (this.workingDirS, this.subDirS[1]).getAbsolutePath();
        String tagAnnotationFileS = new File(tagLibraryDirS, "tag.tas").getAbsolutePath();
        String rawSNPFileS = new File(tagLibraryDirS, "rawSNP.bin").getAbsolutePath();
        String genotypeDirS = new File (this.workingDirS, this.subDirS[3]).getAbsolutePath();
        TagAnnotations tas = new TagAnnotations(tagAnnotationFileS);
        //tas.writeTextFileOfGroup(new File(tagLibraryDirS, "tag.tas.txt").getAbsolutePath(), 799);
        SNPCounts sc = new SNPCounts (rawSNPFileS);
        GBSVCFBuilder builder = new GBSVCFBuilder(tas, sc);
        builder.setTagIdentifyThreshold(3);
        builder.setThreads(8);
        builder.callGenotype(tagBySampleDirS, genotypeDirS);
    }
    
    public void callAllele () {
        String tagLibraryDirS = new File (this.workingDirS, this.subDirS[1]).getAbsolutePath();
        String tagAnnotationFileS = new File(tagLibraryDirS, "tag.tas").getAbsolutePath();
        String alignmentDirS = new File (this.workingDirS, this.subDirS[2]).getAbsolutePath();
        String rawSNPFileS = new File(tagLibraryDirS, "rawSNP.bin").getAbsolutePath();
        String samFileS = new File (alignmentDirS, "tag.sam.gz").getAbsolutePath();
        TagAnnotations tas = new TagAnnotations(tagAnnotationFileS);
        tas.removeAllAllele();
        SNPCounts sc = new SNPCounts (rawSNPFileS);
        int mapQThresh = 30;
        int maxMappingIntervalThresh = 1000;
        tas.callAllele(samFileS, sc, mapQThresh, maxMappingIntervalThresh);
        tas.writeBinaryFile(tagAnnotationFileS);
        //tas.writeTextFile(new File(tagLibraryDirS, "tag.tas.txt").getAbsolutePath());
    }
    
    public void removeLowCountSNP () {
        String tagLibraryDirS = new File (this.workingDirS, this.subDirS[1]).getAbsolutePath();
        String rawSNPFileS = new File(tagLibraryDirS, "rawSNP.bin").getAbsolutePath();
        SNPCounts snpSCs = new SNPCounts(rawSNPFileS);
        int minReadCountAlt = 10;
        snpSCs.writeBinaryFile(rawSNPFileS, minReadCountAlt);
    }
    
    public void callSNP () {
        String tagLibraryDirS = new File (this.workingDirS, this.subDirS[1]).getAbsolutePath();
        String tagAnnotationFileS = new File(tagLibraryDirS, "tag.tas").getAbsolutePath();
        String alignmentDirS = new File (this.workingDirS, this.subDirS[2]).getAbsolutePath();
        String rawSNPFileS = new File(tagLibraryDirS, "rawSNP.bin").getAbsolutePath();
        String samFileS = new File (alignmentDirS, "tag.sam.gz").getAbsolutePath();
        
        int mapQThresh = 30;
        int maxMappingIntervalThresh = 1000;
        int maxDivergence = 7;
        TagAnnotations tas = new TagAnnotations(tagAnnotationFileS);
        tas.removeAllSNP();
        tas.callSNP(samFileS, mapQThresh, maxMappingIntervalThresh, maxDivergence);
        tas.writeBinaryFile(tagAnnotationFileS);
//        tas.writeTextFile(new File(tagLibraryDirS, "tag.tas.txt").getAbsolutePath());
        SNPCounts snpSCs = tas.getSNPCounts();
        snpSCs.writeBinaryFile(rawSNPFileS);
    }
    
    public void alignTags () {
        String tagLibraryDirS = new File (this.workingDirS, this.subDirS[1]).getAbsolutePath();
        String mergedTagAnnotationFileS = new File(tagLibraryDirS, "tag.tas").getAbsolutePath();
        String alignmentDirS = new File (this.workingDirS, this.subDirS[2]).getAbsolutePath();
        new TagAligner(referenceFileS, this.bwaPath, mergedTagAnnotationFileS, alignmentDirS);
    }
    
    public void mergeTagAnnotations () {
        String tagBySampleDirS = new File (this.workingDirS, this.subDirS[0]).getAbsolutePath();
        String tagLibraryDirS = new File (this.workingDirS, this.subDirS[1]).getAbsolutePath();
        String mergedTagCountFileS = new File(tagLibraryDirS, "tag.tas").getAbsolutePath();
        new TagMerger(tagBySampleDirS, mergedTagCountFileS, 1);
    }
    
    public void mkTagsBySample () {
        li = new LibraryInfo(barcodeFileS, libraryFastqMapFileS, this.cutter1, this.cutter2);
        String tagBySampleDirS = new File (this.workingDirS, this.subDirS[0]).getAbsolutePath();
        TagParser tp = new TagParser(li);
        tp.parseFastq(tagBySampleDirS);
        tp.compressTagsBySample(tagBySampleDirS);
    }
    
    public void initializeParameter (String parameterFileS) {
        ArrayList<String> paList = new ArrayList();
        try {
            boolean check = false;
            BufferedReader br = IOUtils.getTextReader(parameterFileS);
            if (!br.readLine().equals("Author: Fei Lu")) check = true;
            if (!br.readLine().equals("Email: flu@genetics.ac.cn; dr.lufei@gmail.com")) check = true;
            if (!br.readLine().equals("Homepage: http://plantgeneticslab.weebly.com/")) check = true;
            if (check) {
                System.out.println("Please keep the author information, or the program quits.");
            }
            String temp = null;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("!Parameter")) {
                    paList.add(br.readLine());
                }
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        this.workingDirS = paList.get(0);
        this.barcodeFileS = paList.get(1);
        this.libraryFastqMapFileS = paList.get(2);
        this.referenceFileS = paList.get(3);
        this.bwaPath = paList.get(4);
        this.cutter1 = paList.get(5).toUpperCase();
        this.cutter2 = paList.get(6).toUpperCase();
        File workingDir = new File(this.workingDirS);
        workingDir.mkdir();
        for (int i = 0; i < this.subDirS.length; i++) {
            File f = new File (this.workingDirS, subDirS[i]);
            f.mkdir();
        }
        System.out.println("Pipeline parameters initialized");
    }
    
    public static void main (String[] args) {
       new LibGBSGo (args[0]);
    }
    
}
