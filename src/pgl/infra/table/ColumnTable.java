/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.table;

import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.List;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

/**
 * A row based implementation of {@link TableInterface}.
 * It is faster on column based operation (e.g. deleting, inserting, adding columns)
 * @author feilu
 * @param <T> 
 */
public class ColumnTable<T> extends TableAbstract<T> {
    /**
     * Constructs a {@link ColumnTable} from a file, with default delimiter of "\t"
     * @param infileS should be text file or ".gz" of a text file
     */
    public ColumnTable (String infileS) {
        if (infileS.endsWith(".gz")) {
            this.readTableFile(infileS, IOFileFormat.TextGzip, "\t");
        }
        else {
            this.readTableFile(infileS, IOFileFormat.Text, "\t");
            String a;
        }
    }
    
    /**
     * Constructs a {@link ColumnTable} from a file with specified format, with default delimiter of "\t"
     * @param infileS
     * @param format 
     */
    public ColumnTable (String infileS, IOFileFormat format) {
        this.readTableFile(infileS, format, "\t");
    }
    
    /**
     * Constructs a {@link ColumnTable} using two dimension list of custom objects
     * @param header
     * @param cells 
     */
    public ColumnTable (List<String> header, List<List<T>> cells) {
        this.header = header;
        this.cells = cells;
    }
    
    private void readTableFile (String infileS, IOFileFormat format, String delimiter) {
        try {
            BufferedReader br = null;
            if (format == IOFileFormat.Text) {
                br = IOUtils.getTextReader(infileS);
            }
            else if (format == IOFileFormat.TextGzip) {
                br = IOUtils.getTextGzipReader(infileS);
            }
            else {
                throw new UnsupportedOperationException("Unsupported format for input");
            }
            String temp = br.readLine();
            this.header = PStringUtils.fastSplit(temp, delimiter);
            this.buildHIMap();
            cells = new ArrayList<>();
            for (int i = 0; i < header.size(); i++) {
                cells.add(new ArrayList<T>());
            }
            while ((temp = br.readLine()) != null) {
                List<T> current = (List<T>)PStringUtils.fastSplit(temp, delimiter);
                for (int i = 0; i < current.size(); i++) {
                    cells.get(i).add(current.get(i));
                }
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    @Override
    public T getCell(int rowIndex, int columnIndex) {
        return this.getColumn(columnIndex).get(rowIndex);
    }
    
    @Override
    public int getRowNumber() {
        return cells.get(0).size();
    }

    @Override
    public void sortAsText(int columnIndex) {
        throw new UnsupportedOperationException("Not supported for a column based table");
    }

    @Override
    public boolean sortAsNumber(int columnIndex) {
        throw new UnsupportedOperationException("Not supported for a column based table");
    }

    @Override
    public void insertColumn(String columnName, int columnIndex, List<T> cList) {
        if (cList.size() != this.getRowNumber()) {
            System.out.println("The column to be inserted has different number of rows. Program quits");
            System.exit(1);
        }
        this.header.add(columnIndex, columnName);
        cells.add(columnIndex, cList);
        this.buildHIMap();
    }

    @Override
    public void removeColumn(int columnIndex) {
        this.header.remove(columnIndex);
        for (int i = 0; i < this.getRowNumber(); i++) {
            cells.get(i).remove(columnIndex);
        }
        cells.remove(columnIndex);
        this.buildHIMap();
    }
    
    @Override
    public void insertRow(int rowIndex, List<T> rList) {
        if (rList.size() != this.getColumnNumber()) {
            System.out.println("The row to be inserted has different number of columns. Program quits");
            System.exit(1);
        }
        for (int i = 0; i < this.getColumnNumber(); i++) {
            cells.get(i).add(rowIndex, rList.get(i));
        }
        cells.add(rowIndex, rList);
    }

    @Override
    public void removeRow(int rowIndex) {
        for (int i = 0; i < this.getColumnNumber(); i++) {
            cells.get(i).remove(rowIndex);
        }
    }

    @Override
    public void setCell(int rowIndex, int columnIndex, T c) {
        cells.get(columnIndex).set(rowIndex, c);
    }

    @Override
    public void setRow(int rowIndex, List<T> rList) {
        for (int i = 0; i < this.getColumnNumber(); i++) {
            cells.get(i).set(rowIndex, rList.get(i));
        }
    }

    @Override
    public void setColumn(int columnIndex, List<T> cList) {
        cells.set(columnIndex, cList);
    }

    @Override
    public void addColumn(String columnName, List cList) {
        if (cList.size() != this.getRowNumber()) {
            System.out.println("The column to be added has different number of rows. Program quits");
            System.exit(1);
        }
        this.header.add(columnName);
        this.cells.add(cList);
        this.hiMap.put(columnName, this.getColumnNumber()-1);
    }

    @Override
    public void addRow(List<T> rList) {
        if (rList.size() != this.getColumnNumber()) {
            System.out.println("The row to be added has different number of columns. Program quits");
            System.exit(1);
        }
        for (int i = 0; i < this.getColumnNumber(); i++) {
            cells.get(i).add(rList.get(i));
        }
    }

    @Override
    public List<T> getColumn(int columnIndex) {
        return cells.get(columnIndex);
    }

    @Override
    public List<T> getRow(int rowIndex) {
        List<T> rList = new ArrayList();
        for (int i = 0; i < this.getColumnNumber(); i++) {
            rList.add(cells.get(i).get(rowIndex));
        }
        return rList;
    }
}
