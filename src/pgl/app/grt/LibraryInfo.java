/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.app.grt;

import pgl.infra.table.RowTable;
import gnu.trove.list.array.TIntArrayList;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 *
 * @author feilu
 */
public class LibraryInfo {
    
    String[] libs = null;
    String[][] taxaNames = null;
    String[][] barcodeR1 = null;
    String[][] barcodeR2 = null;
    HashMap<String, Set<String>>[] barcodeR1TaxaMaps = null;
    HashMap<String, Set<String>>[] barcodeR2TaxaMaps = null;
    String[] libFastqsR1 = null;
    String[] libFastqsR2 = null;
    String cutter1 = null;
    String cutter2 = null;
    
    public LibraryInfo (String barcodeFileS, String libFastqMapFileS, String cutter1, String cutter2) {
        this.parseBarcode(barcodeFileS, libFastqMapFileS);
        this.cutter1 = cutter1;
        this.cutter2 = cutter2;
    }
    
    public int getLibNumber () {
        return libs.length;
    }
    
    public String[] getLibArray () {
        String[] na = new String[libs.length];
        System.arraycopy(libs, 0, na, 0, libs.length);
        return na;
    }
    
    public String getCutter1 () {
        return this.cutter1;
    }
    
    public String getCutter2 () {
        return this.cutter2;
    }
    
    public String getLibName (int index) {
        return libs[index];
    }
    
    public String[] getTaxaNames (int index) {
        return taxaNames[index];
    }
    
    public String[] getLibBarcodeR1 (int index) {
        String[] na = new String[barcodeR1[index].length];
        System.arraycopy(barcodeR1[index], 0, na, 0, barcodeR1[index].length);
        return na;
    }
    
    public String[] getLibBarcodeR2 (int index) {
        String[] na = new String[barcodeR2[index].length];
        System.arraycopy(barcodeR2[index], 0, na, 0, barcodeR2[index].length);
        return na;
    }
    
    public HashMap<String, Set<String>> getbarcodeR1TaxaMap (int index) {
        return barcodeR1TaxaMaps[index];
    }
    
    public HashMap<String, Set<String>> getbarcodeR2TaxaMap (int index) {
        return barcodeR2TaxaMaps[index];
    }
    
    public String getFastqFileSR1 (int index) {
        return libFastqsR1[index];
    }
    
    public String getFastqFileSR2 (int index) {
        return libFastqsR2[index];
    }
    
    private void parseBarcode (String barcodeFileS, String libFastqMapFileS) {
        RowTable<String> t = new RowTable<>(barcodeFileS);
        Set<String> s = new HashSet<>();
        for (int i = 0; i < t.getRowNumber(); i++) {
            StringBuilder sb = new StringBuilder();
            sb.append(t.getCell(i, 1)).append("_").append(t.getCell(i, 2)).append("_").append(t.getCell(i, 3));
            s.add(sb.toString());
        }
        libs = s.toArray(new String[s.size()]);
        Arrays.sort(libs);
        taxaNames = new String[libs.length][];
        barcodeR1 = new String[libs.length][];
        barcodeR2 = new String[libs.length][];
        barcodeR1TaxaMaps = new HashMap[libs.length];
        barcodeR2TaxaMaps = new HashMap[libs.length];
        
        libFastqsR1 = new String[libs.length];
        libFastqsR2 = new String[libs.length];
        for (int i = 0; i < libs.length; i++) {
            List<String> nameList = new ArrayList<>();
            List<String> barcodeR1List = new ArrayList<>();
            List<String> barcodeR2List = new ArrayList<>();
            for (int j = 0; j < t.getRowNumber(); j++) {
                StringBuilder sb = new StringBuilder();
                sb.append(t.getCell(j, 1)).append("_").append(t.getCell(j, 2)).append("_").append(t.getCell(j, 3));
                if (!sb.toString().equals(libs[i])) continue;
                sb = new StringBuilder();
                sb.append(t.getCell(j, 0)).append("_").append(t.getCell(j, 1)).append("_").append(t.getCell(j, 2)).append("_").append(t.getCell(j, 3)).append("_").append(t.getCell(j, 4));
                nameList.add(sb.toString());
                barcodeR1List.add(t.getCell(j, 5));
                barcodeR2List.add(t.getCell(j, 6)); 
            }
            taxaNames[i] = nameList.toArray(new String[nameList.size()]);
            barcodeR1[i] = barcodeR1List.toArray(new String[nameList.size()]);
            barcodeR2[i] = barcodeR2List.toArray(new String[nameList.size()]);
            barcodeR1TaxaMaps[i] = new HashMap<>();
            barcodeR2TaxaMaps[i] = new HashMap<>();
            for (int j = 0; j < barcodeR1List.size(); j++) {
                s = new HashSet<>();
                barcodeR1TaxaMaps[i].put(barcodeR1List.get(j), s);
                s = new HashSet<>();
                barcodeR2TaxaMaps[i].put(barcodeR2List.get(j), s);
            }
            for (int j = 0; j < taxaNames[i].length; j++) {
                s = barcodeR1TaxaMaps[i].get(barcodeR1[i][j]);
                s.add(taxaNames[i][j]);
                barcodeR1TaxaMaps[i].put(barcodeR1[i][j], s);
                s = barcodeR2TaxaMaps[i].get(barcodeR2[i][j]);
                s.add(taxaNames[i][j]);
                barcodeR2TaxaMaps[i].put(barcodeR2[i][j], s);
            }
        }
        t = new RowTable<>(libFastqMapFileS);
        List<String> lList = new ArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            StringBuilder sb = new StringBuilder();
            sb.append(t.getCell(i, 0)).append("_").append(t.getCell(i, 1)).append("_").append(t.getCell(i, 2));
            lList.add(sb.toString());
        }
        Collections.sort(lList);
        TIntArrayList availableIndexList = new TIntArrayList();
        for (int i = 0; i < libs.length; i++) {
            int index = Collections.binarySearch(lList, libs[i]);
            if (index < 0) {
                System.out.println(libs[i] + " does not have corresponding fastqs");
            }
            else {
                if (t.getCell(index, 3).equals("NA") || t.getCell(index, 4).equals("NA")) {
                    System.out.println(libs[i] + " does not have corresponding fastqs");
                }
                else {
                    this.libFastqsR1[i] = t.getCell(index, 3);
                    this.libFastqsR2[i] = t.getCell(index, 4);
                    availableIndexList.add(i);
                }
            } 
        }
        List<String> libList = new ArrayList();
        List<String[]> taxaNameList = new ArrayList();
        List<String[]> barcodeR1List = new ArrayList();
        List<String[]> barcodeR2List = new ArrayList();
        List<HashMap<String, Set<String>>> r1MapList = new ArrayList();
        List<HashMap<String, Set<String>>> r2MapList = new ArrayList();
        List<String> libFqR1List = new ArrayList();
        List<String> libFqR2List = new ArrayList();
        int[] aIndex = availableIndexList.toArray();
        for (int i = 0; i < aIndex.length; i++) {
            libList.add(libs[aIndex[i]]);
            taxaNameList.add(taxaNames[aIndex[i]]);
            barcodeR1List.add(barcodeR1[aIndex[i]]);
            barcodeR2List.add(barcodeR1[aIndex[i]]);
            r1MapList.add(barcodeR1TaxaMaps[aIndex[i]]);
            r2MapList.add(barcodeR2TaxaMaps[aIndex[i]]);
            libFqR1List.add(libFastqsR1[aIndex[i]]);
            libFqR2List.add(libFastqsR2[aIndex[i]]);
        }
        
        libs = libList.toArray(new String[libList.size()]);
        taxaNames = taxaNameList.toArray(new String[taxaNameList.size()][]);
        barcodeR1 = barcodeR1List.toArray(new String[barcodeR1List.size()][]);
        barcodeR2 = barcodeR2List.toArray(new String[barcodeR2List.size()][]);
        barcodeR1TaxaMaps = r1MapList.toArray(new HashMap[r1MapList.size()]);
        barcodeR2TaxaMaps = r2MapList.toArray(new HashMap[r2MapList.size()]);
        libFastqsR1 = libFqR1List.toArray(new String[libFqR1List.size()]);
        libFastqsR2 = libFqR2List.toArray(new String[libFqR2List.size()]);
        
        System.out.println(libs.length+" libraries will be paralell processd. They are:");
        int cnt = 0;
        for (int i = 0; i < libs.length; i++) {
            StringBuilder sb = new StringBuilder();
            sb.append(libs[i]).append(" with ").append(this.taxaNames[i].length).append(" samples");
            System.out.println(sb.toString());
            cnt+=this.taxaNames[i].length;
        }
        System.out.println("A total of " + String.valueOf(cnt) + " samples are in the current batch");
    }
    
}
