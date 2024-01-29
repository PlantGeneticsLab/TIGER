package pgl.app.illuminaQC;

import pgl.infra.align.g2.ShortreadPEAlignment;
import pgl.infra.dna.FastaByte;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import pgl.graph.r.DensityPlot;
import pgl.infra.utils.IOUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.InputStreamReader;
import java.util.*;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

public class IlluminaQCGo {

    String referencePath = null;
    int mitoChrom = -1;
    int chloroChrom = -1;
    String fastQCPath = null;
    String bwaPath = null;
    String rPath = null;
    String fastqPath = null;
    String workingPath = null;
    boolean ifCoverage = false;
    String[] subDir = {"sampledFastq", "fastQC", "alignment", "insertSize", "sampledFasta"};

    String[][] result = null;
    int fastQCColNum = 0;
    double genomeSize = 0;

    public IlluminaQCGo (String parameterFileS) {
        this.initializeParameter(parameterFileS);
        this.mkSubDirectories();
        this.sampleFastq();
        this.fastQC();
        this.alignBWA();
        this.mkSummary();
    }

    private void mkSummary () {
        String inputDirS = new File (workingPath, subDir[1]).getAbsolutePath();
        String imageDirS = new File (workingPath, subDir[3]).getAbsolutePath();
        String outfileS = new File (workingPath, "qcSummary.txt").getAbsolutePath();

        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".zip");

        String[] column = null;
        String[] fileNames = new String[fs.length];
        String[] addColumn = {"Insert size", "Mitochondria proportion", "Chloroplast proportion", "Coverage"};

        for (int i = 0; i < fs.length; i++) {
            try {
                ZipFile zf = new ZipFile(fs[i].getAbsolutePath());
                Enumeration<? extends ZipEntry> entries = zf.entries();
                while (entries.hasMoreElements()) {
                    ZipEntry ze = entries.nextElement();
                    if (!ze.getName().endsWith("summary.txt")) continue;
                    if (i == 0) {
                        BufferedReader br = new BufferedReader(new InputStreamReader(zf.getInputStream(ze), "UTF-8"));
                        String temp = null;
                        ArrayList<String> columnList = new ArrayList();
                        while ((temp = br.readLine()) != null) {
                            String[] tem = temp.split("\t");
                            columnList.add(tem[1]);
                        }
                        fastQCColNum = columnList.size();
                        for (int j = 0; j < addColumn.length; j++) {
                            columnList.add(addColumn[j]);
                        }
                        column = columnList.toArray(new String[columnList.size()]);
                        result = new String[fs.length][column.length];
                        br.close();
                    }
                    BufferedReader br = new BufferedReader(new InputStreamReader(zf.getInputStream(ze), "UTF-8"));
                    String temp = null;
                    int cnt = 0;
                    while ((temp = br.readLine()) != null) {
                        String[] tem = temp.split("\t");
                        result[i][cnt] = tem[0];
                        fileNames[i] = tem[2];
                        cnt++;
                    }
                    br.close();
                    break;
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        ArrayList<String> taxaList = new ArrayList();
        HashMap<String, Integer> taxonIndexMap = new HashMap();
        int cnt = 0;
        for (int i = 0; i < fileNames.length; i+=2) {
            String taxon = fileNames[i].split("_")[0];
            taxaList.add(taxon);
            taxonIndexMap.put(taxon, cnt);
            cnt++;
        }

        if (ifCoverage) {
            FastaByte genome = new FastaByte(referencePath);
            genomeSize = genome.getTotalSeqLength();
            System.out.println("Genome size is " +String.valueOf(genomeSize)+" bp");
        }

        taxaList.parallelStream().forEach(taxon -> {
            String samFileS = new File (new File (workingPath, subDir[2]).getAbsolutePath(), taxon+".pe.sam").getAbsolutePath();
            String fastqFileS = new File (fastqPath, taxon+"_1.fq.gz").getAbsolutePath();
            String imageFileS = new File (imageDirS, taxon+".pdf").getAbsolutePath();
            int maxSampleSize = 10000;
            int taxonIndex = taxonIndexMap.get(taxon);
            int maxLength = 1000;
            ShortreadPEAlignment sa = new ShortreadPEAlignment ();
            sa.readFromBWAMEM(samFileS);
            TIntArrayList sizeList = new TIntArrayList();
            for (int i = 0; i < sa.getAlignmentNumber(); i++) {
                if (sa.getMappingQualityF(i) < 30) continue;
                if (sa.getMappingQualityB(i) < 30) continue;
                if (!sa.getHitF(i).equals(sa.getHitB(i))) continue;
                sizeList.add(sa.getPEFragmentSize(i));
            }
            int[] fragSize = sizeList.toArray();
            TDoubleArrayList sList = new TDoubleArrayList();
            int count = 0;
            for (int i = 0; i < fragSize.length; i++) {
                if (fragSize[i] > maxLength) continue;
                if (fragSize[i] < 0) continue;
                if (count > maxSampleSize) break;
                sList.add(fragSize[i]);
            }
            DensityPlot d = new DensityPlot (sList.toArray());
            d.setXLim(0, 1000);
            d.setRPath(rPath);
            d.setTitle(taxon);
            d.setYLabel("Density");
            d.setXLabel("Insert size of library");
            d.setSmoothN(10000);
            d.saveGraph(imageFileS);
            System.out.println(taxon + " density plot saved at " + imageFileS);
            int[] sizeDis = new int[maxLength+1];
            for (int i = 0; i < fragSize.length; i++) {
                if (fragSize[i] > maxLength) continue;
                if (fragSize[i] < 0) continue;
                sizeDis[fragSize[i]]++;
            }
            int max = -1;
            int insertSize = -1;
            for (int i = 0; i < sizeDis.length; i++) {
                if (sizeDis[i] > max) {
                    max = sizeDis[i];
                    insertSize = i;
                }
            }
            System.out.println(insertSize);
            result[taxonIndex*2][fastQCColNum] = String.valueOf(insertSize);
            result[taxonIndex*2+1][fastQCColNum] = String.valueOf(insertSize);

            String mitoS = String.valueOf(mitoChrom);
            String chloS = String.valueOf(chloroChrom);
            double mitoCnt = 0;
            double chloCnt = 0;
            for (int j = 0; j < sa.getAlignmentNumber(); j++) {
                if (sa.getHitF(j).equals(mitoS)) mitoCnt++;
                if (sa.getHitF(j).equals(chloS)) chloCnt++;
            }
            double mitoProportion = mitoCnt/sa.getAlignmentNumber();
            result[taxonIndex*2][fastQCColNum+1] = String.valueOf(mitoProportion);
            result[taxonIndex*2+1][fastQCColNum+1] = result[taxonIndex*2][fastQCColNum+1];
            double chloroProportion = chloCnt/sa.getAlignmentNumber();
            result[taxonIndex*2][fastQCColNum+2] = String.valueOf(chloroProportion);
            result[taxonIndex*2+1][fastQCColNum+2] = result[taxonIndex*2][fastQCColNum+2];
            System.out.println(taxon + " mitochondria proportion is " + String.valueOf(mitoProportion));
            System.out.println(taxon + " chloroProportion proportion is " + String.valueOf(chloroProportion));

            if (ifCoverage) {
                try {
                    BufferedReader br = IOUtils.getTextGzipReader(fastqFileS);
                    String temp = null;
                    int readNum = 0;
                    int readLength = 0;
                    while ((temp = br.readLine()) != null) {
                        temp = br.readLine();
                        if (readNum == 0) {
                            readLength = temp.length();
                        }
                        readNum++;
                        br.readLine();br.readLine();
                        if (readNum%10000000==0) System.out.println("Current read number in "+ taxon + " is "+ String.valueOf(readNum));
                    }
                    br.close();
                    double coverage = (double)readNum/genomeSize*readLength*2;
                    result[taxonIndex*2][fastQCColNum+3] = String.valueOf(coverage);
                    result[taxonIndex*2+1][fastQCColNum+3] = String.valueOf(coverage);
                    System.out.println("Coverage of " + taxon + " is " + String.valueOf(coverage));
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
            }
            else {
                result[taxonIndex*2][fastQCColNum+3] = "NA";
                result[taxonIndex*2+1][fastQCColNum+3] = "NA";
            }
        });
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            StringBuilder sb = new StringBuilder("FileName");
            for (int i = 0; i < column.length; i++) {
                sb.append("\t").append(column[i]);
            }
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < result.length; i++) {
                sb = new StringBuilder(fileNames[i]);
                for (int j = 0; j < column.length; j++) {
                    sb.append("\t").append(result[i][j]);
                }
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void alignBWA () {
        String inputDirS = new File (workingPath, subDir[0]).getAbsolutePath();
        String outputDirS = new File (workingPath, subDir[2]).getAbsolutePath();
        String outputPerlS = new File (workingPath, "runAlign.pl").getAbsolutePath();
        File[] fs = new File (inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".gz");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            nameSet.add(fs[i].getName().split("_")[0]);
        }
        String[] names = nameSet.toArray(new String[nameSet.size()]);
        Arrays.sort(names);
        //bwa mem -t  nthreads ref.fa read1.fq read2.fq > aln-pe.sam
        int numCores = Runtime.getRuntime().availableProcessors();
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outputPerlS);
            for (int i = 0; i < names.length; i++) {
                StringBuilder sb = new StringBuilder(bwaPath+" mem ");
                sb.append(" -t ").append(numCores).append(" ").append(referencePath).append(" ");
                sb.append(new File(inputDirS, names[i]+"_1.fq.gz").getAbsolutePath()).append(" ");
                sb.append(new File(inputDirS, names[i]+"_2.fq.gz").getAbsolutePath());
                sb.append(" > ").append(new File(outputDirS, names[i]+".pe.sam").getAbsolutePath());
                String cmd = sb.toString();
                String command = "system(\""+cmd+"\");";
                bw.write(command);
                bw.newLine();

            }
            bw.flush();
            bw.close();
            StringBuilder sb = new StringBuilder("perl ");
            sb.append(outputPerlS);
            Runtime run = Runtime.getRuntime();
            Process p = run.exec(sb.toString());
            System.out.println(sb.toString());
            BufferedReader br = new BufferedReader(new InputStreamReader(p.getErrorStream()));
            String temp = null;
            while ((temp = br.readLine()) != null) {
                System.out.println(temp);
            }
            p.waitFor();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        new File (outputPerlS).delete();
    }

    private void fastQC () {
        String inputDirS = new File (workingPath, subDir[0]).getAbsolutePath();
        String outputDirS = new File (workingPath, subDir[1]).getAbsolutePath();
        File[] fs = new File (inputDirS).listFiles();
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            StringBuilder sb = new StringBuilder(fastQCPath);
            sb.append(" ").append(f.getAbsolutePath()).append(" -o ").append(outputDirS);
            String cmd = sb.toString();
            System.out.println(cmd);
            try {
                Runtime run = Runtime.getRuntime();
                Process p = run.exec(cmd);
                BufferedReader br = new BufferedReader(new InputStreamReader(p.getErrorStream()));
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    System.out.println(temp);
                }
                p.waitFor();
                System.out.println(f.getName()+ " Fastqc evalutation is finished at" + outputDirS);
            }
            catch (Exception e) {
                e.printStackTrace();
            }

        });
    }

    private void sampleFastq () {
        String infileDirS = fastqPath;
        String outputDirS = new File (workingPath, subDir[0]).getAbsolutePath();
        String outputFastaDirS = new File (workingPath, subDir[4]).getAbsolutePath();
        int readNum = 100000;
        int startPoint = 100000;
        int fastaNum = 1000;
        File[] fs = new File(infileDirS).listFiles();
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().split("_")[0]);
        }
        nameSet.parallelStream().forEach(name -> {
            String infile1 = new File (infileDirS, name+"_1.fq.gz").getAbsolutePath();
            String infile2 = new File (infileDirS, name+"_2.fq.gz").getAbsolutePath();
            String outfile1 = new File (outputDirS, name+"_1.fq.gz").getAbsolutePath();
            String outfile2 = new File (outputDirS, name+"_2.fq.gz").getAbsolutePath();
            String outfileFasta = new File (outputFastaDirS, name+"_1.fa").getAbsolutePath();
            try {
                BufferedReader br1 = IOUtils.getTextGzipReader(infile1);
                BufferedReader br2 = IOUtils.getTextGzipReader(infile2);
                BufferedWriter bw1 = IOUtils.getTextGzipWriter(outfile1);
                BufferedWriter bw2 = IOUtils.getTextGzipWriter(outfile2);
                BufferedWriter bwf = IOUtils.getTextGzipWriter(outfileFasta);
                String temp = null;
                int cnt = 0;
                while ((temp = br1.readLine()) != null) {
                    cnt++;
                    if (cnt < startPoint) {
                        br1.readLine();br1.readLine();br1.readLine();
                        br2.readLine();br2.readLine();br2.readLine();br2.readLine();
                    }
                    else {
                        bw1.write(temp+"\n");bw1.write(br1.readLine()+"\n");bw1.write(br1.readLine()+"\n");bw1.write(br1.readLine()+"\n");
                        bw2.write(br2.readLine()+"\n");bw2.write(br2.readLine()+"\n");bw2.write(br2.readLine()+"\n");bw2.write(br2.readLine()+"\n");
                        for (int i = 0; i < readNum-1; i++) {
                            bw1.write(br1.readLine()+"\n");
                            temp = br1.readLine();bw1.write(temp+"\n");
                            bw1.write(br1.readLine()+"\n");bw1.write(br1.readLine()+"\n");
                            bw2.write(br2.readLine()+"\n");bw2.write(br2.readLine()+"\n");bw2.write(br2.readLine()+"\n");bw2.write(br2.readLine()+"\n");
                            if (i > fastaNum) continue;
                            bwf.write(">"+String.valueOf(i));
                            bwf.newLine();
                            bwf.write(temp);
                            bwf.newLine();
                        }
                        bw1.flush();bw1.close();
                        bw2.flush();bw2.close();
                        bwf.flush();bwf.close();
                        br1.close();
                        br2.close();
                        break;
                    }
                }
                System.out.println(String.valueOf(readNum) + " reads are sampled from"+ name);
            }
            catch (Exception e) {
                e.printStackTrace();
            }

        });
    }

    private void mkSubDirectories () {
        for (int i = 0; i < subDir.length; i++) {
            new File(workingPath, subDir[i]).mkdir();
        }
    }

    private void initializeParameter (String parameterFileS) {
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
        referencePath = paList.get(0);
        mitoChrom = Integer.valueOf(paList.get(1));
        chloroChrom = Integer.valueOf(paList.get(2));
        fastQCPath = paList.get(3);
        bwaPath = paList.get(4);
        rPath = paList.get(5);
        fastqPath = paList.get(6);
        workingPath = paList.get(7);
        if (paList.get(8).startsWith("T")) ifCoverage = true;
    }

    public static void main (String[] args) {
        new IlluminaQCGo (args[0]);
    }
}
