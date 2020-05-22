/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.anno.gene;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import pgl.infra.utils.IOUtils;


/**
 *
 * @author feilu
 */
public class GFFUtils {
    
    /**
     * Produce a map file for gene classification
     * @param inputGFF3
     * @param geneAnnotationFileS 
     */
    public static void mkMaizeGeneAnnotationFile (String inputGFF3, String geneAnnotationFileS) {
        try {
            BufferedReader br = IOUtils.getTextGzipReader(inputGFF3);
            BufferedWriter bw = IOUtils.getTextWriter(geneAnnotationFileS);
            bw.write("Gene\tBiotype\tChr\tStart\tEnd\tStrand\tDescription");
            bw.newLine();
            String temp = null;
            StringBuilder sb = new StringBuilder();
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("###")) {
                    sb = new StringBuilder();
                    temp = br.readLine();
                    if (temp == null) {
                        break;
                    }
                    String[] tem = temp.split(";");
                    String[] te = tem[0].split("\t");
                    if (te[2].startsWith("chromo")) continue;
                    sb.append(tem[0].split(":")[1]).append("\t");
                    String s = "NA";
                    for (int i = 1; i < tem.length; i++) {
                        if (tem[i].startsWith("biotype")) {
                            s = tem[i].replaceFirst("biotype=", "");
                            break;
                        }
                    }
                    sb.append(s).append("\t");
                    sb.append(te[0]).append("\t").append(te[3]).append("\t").append(te[4]).append("\t").append(te[6]).append("\t");
                    s = "NA";
                    for (int i = 1; i < tem.length; i++) {
                        if (tem[i].startsWith("description")) {
                            s = tem[i].replaceFirst("description=", "");
                            break;
                        }
                    }
                    sb.append(s);
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.out.println();
        }
    }
    
    /**
     * Standardize maize GFF3 (remove annotation of contigs, change Mt to 11, Pt to 12)
     * @param inputGFF3
     * @param outputGFF3 
     */
    public static void modifyMaizeAGPV4GFF3 (String inputGFF3, String outputGFF3) {
        try {
            BufferedReader br = IOUtils.getTextGzipReader(inputGFF3);
            BufferedWriter bw = IOUtils.getTextWriter(outputGFF3);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("B73")) continue;
                else if (temp.startsWith("##sequence-region   B")) {
                    continue;
                }
                else if (temp.startsWith("##sequence-region   Mt")) {
                    temp = temp.replaceFirst("##sequence-region   Mt", "##sequence-region   11(Mt)");
                }
                else if (temp.startsWith("##sequence-region   Pt")) {
                    temp = temp.replaceFirst("##sequence-region   Pt", "##sequence-region   12(Pt)");
                }
                else if (temp.startsWith("Mt")) {
                    temp = temp.replaceFirst("Mt", "11");
                }
                else if (temp.startsWith("Pt")) {
                    temp = temp.replaceFirst("Pt", "12");
                }
                bw.write(temp);
                bw.newLine();
            }
            bw.flush();
            bw.close();
            br.close();
        }
        catch (Exception e) {
            
        }
    }
    
}
