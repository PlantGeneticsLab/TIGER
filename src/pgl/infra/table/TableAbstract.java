/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.table;

import gnu.trove.list.array.TDoubleArrayList;
import java.io.BufferedWriter;
import java.util.HashMap;
import java.util.List;

import gnu.trove.list.array.TIntArrayList;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;

/**
 * Defining fields of table and implementing part of the methods
 * @author feilu
 * @param <T>
 */
public abstract class TableAbstract<T> implements TableInterface<T> {
    protected List<String> header = null;
    protected List<List<T>> cells = null;
    protected HashMap<String, Integer> hiMap = new HashMap<>();
    protected int sortColumnIndex  = -1;
    
    protected void buildHIMap () {
        for (int i = 0; i < header.size(); i++) {
            hiMap.put(header.get(i), i);
        }
    }
    
    @Override
    public List<String> getHeader() {
        return this.header;
    }
    
    @Override
    public String getColumnName (int columnIndex) {
        return header.get(columnIndex);
    }
    
    @Override
    public double[] getColumnAsDoubleArray (int columnIndex) {
        TDoubleArrayList l = new TDoubleArrayList();
        Double d = null;
        for (int i = 0; i < this.getRowNumber(); i++) {
            d = this.getCellAsDouble(i, columnIndex);
            if (d == null) return null;
            l.add(d);
        }
        return l.toArray();
    }

    public int[] getColumnAsIntArray (int columnIndex) {
        TIntArrayList l = new TIntArrayList();
        Integer in = null;
        for (int i = 0; i < this.getRowNumber(); i++) {
            in = this.getCellAsInteger(i, columnIndex);
            if (in == null) return null;
            l.add(in);
        }
        return l.toArray();

    }

    @Override
    public String getCellAsString (int rowIndex, int columnIndex) {
        return this.getCell(rowIndex, columnIndex).toString();
    }
    
    @Override
    public Double getCellAsDouble (int rowIndex, int columnIndex) {
        T ob = this.getCell(rowIndex, columnIndex);
        if (ob instanceof Number) {
            return ((Number) ob).doubleValue();
        }
        else if (ob instanceof String) {
            return Double.parseDouble((String)ob);
        }
        else {
            return null;
        }
    }
    
    @Override
    public Integer getCellAsInteger (int rowIndex, int columnIndex) {
        T ob = this.getCell(rowIndex, columnIndex);
        if (ob instanceof Number) {
            return ((Number) ob).intValue();
        }
        else if (ob instanceof String) {
            return Integer.parseInt((String)ob);
        }
        else {
            return null;
        }
    }
    
    @Override
    public int getColumnNumber () {
        return this.header.size();
    }
    
    @Override
    public int getColumnIndex(String columnName) {
        return hiMap.get(columnName);
    }
    
    @Override
    public void sortAsText(String columnName) {
        int columnIndex = this.hiMap.get(columnName);
        this.sortColumnIndex = columnIndex;
        this.sortAsText(columnIndex);
    }

    @Override
    public boolean sortAsNumber (String columnName) {
        int columnIndex = this.hiMap.get(columnName);
        this.sortColumnIndex = columnIndex;
        return this.sortAsNumber(columnIndex);
    }
    
    @Override
    public void removeColumn(String columnName) {
        int columnIndex = this.hiMap.get(columnName);
        this.removeColumn(columnIndex);
    }
    
    @Override
    public void writeTextTable (String outfileS, IOFileFormat format) {
        try {
            BufferedWriter bw = null;
            if (format == IOFileFormat.Text) {
                bw = IOUtils.getTextWriter(outfileS);
            }
            else if (format == IOFileFormat.TextGzip) {
                bw = IOUtils.getTextGzipWriter(outfileS);
            }
            else {
                throw new UnsupportedOperationException("Unsupported format for input");
            }
            StringBuilder sb = new StringBuilder(header.get(0));
            for (int i = 1; i < this.getColumnNumber(); i++) {
                sb.append("\t").append(header.get(i));
            }
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < this.getRowNumber(); i++) {
                sb = new StringBuilder(this.getCellAsString(i, 0));
                for (int j = 1; j < this.getColumnNumber(); j++) {
                    sb.append("\t").append(this.getCellAsString(i, j));
                }
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
            System.out.println("Table is written to " + outfileS);
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    @Override
    public void writeTextTable (String outfileS, IOFileFormat format, boolean[] ifOut) {
        try {
            BufferedWriter bw = null;
            if (format == IOFileFormat.Text) {
                bw = IOUtils.getTextWriter(outfileS);
            }
            else if (format == IOFileFormat.TextGzip) {
                bw = IOUtils.getTextGzipWriter(outfileS);
            }
            else {
                throw new UnsupportedOperationException("Unsupported format for input");
            }
            StringBuilder sb = new StringBuilder(header.get(0));
            for (int i = 1; i < this.getColumnNumber(); i++) {
                sb.append("\t").append(header.get(i));
            }
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < this.getRowNumber(); i++) {
                if (!ifOut[i]) continue;
                sb = new StringBuilder(this.getCellAsString(i, 0));
                for (int j = 1; j < this.getColumnNumber(); j++) {
                    sb.append("\t").append(this.getCellAsString(i, j));
                }
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
            System.out.println("Table is written to " + outfileS);
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
}
