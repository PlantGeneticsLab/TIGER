/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.table;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.List;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

/**
 * A row based implementation of {@link TableInterface}.
 * It is faster on row based operation (e.g. deleting, inserting, adding rows)
 * @author feilu
 * @param <T> 
 */
public class RowTable<T> extends TableAbstract<T> {
    /**
     * Constructs a {@link RowTable} from a file, with default delimiter of "\t"
     * @param infileS should be text file or ".gz" of a text file
     */
    public RowTable (String infileS) {
        if (infileS.endsWith(".gz")) {
            this.readTableFile(infileS, IOFileFormat.TextGzip, "\t");
        }
        else {
            this.readTableFile(infileS, IOFileFormat.Text, "\t");
        }
    }
    
    /**
     * Constructs a {@link RowTable} from a file, with a custom delimiter
     * @param infileS
     * @param delimiter 
     */
    public RowTable (String infileS, String delimiter) {
        if (infileS.endsWith(".gz")) {
            this.readTableFile(infileS, IOFileFormat.TextGzip, delimiter);
        }
        else {
            this.readTableFile(infileS, IOFileFormat.Text, delimiter);
        }
    }
    
    /**
     * Constructs a {@link RowTable} from a file with specified format, with default delimiter of "\t"
     * @param infileS
     * @param format 
     */
    public RowTable (String infileS, IOFileFormat format) {
        this.readTableFile(infileS, format, "\t");
    }
    
    /**
     * Constructs a {@link RowTable} using two dimension list of custom objects
     * @param header
     * @param cells 
     */
    public RowTable (List<String> header, List<List<T>> cells) {
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
            while ((temp = br.readLine()) != null) {
                List<T> current = (List<T>)PStringUtils.fastSplit(temp, delimiter);
                cells.add(current);
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
        return this.getRow(rowIndex).get(columnIndex);
    }
    
    @Override
    public int getRowNumber() {
        return cells.size();
    }

    @Override
    public void sortAsText(int columnIndex) {
        this.sortColumnIndex = columnIndex;
        GenericSorting.quickSort(0, this.getRowNumber(), compTextByColumn, swapper);
        System.out.println("Table is sorted based on the text of column " + header.get(columnIndex));
    }

    @Override
    public boolean sortAsNumber(int columnIndex) {
        for (int i = 0; i < this.getRowNumber(); i++) {
            T ob = this.getCell(i, columnIndex);
            if (!(ob instanceof String)) {
                System.out.println("Column of " + header.get(columnIndex) + " can not all be converted to number");
                return false;
            }
        }
        this.sortColumnIndex = columnIndex;
        GenericSorting.quickSort(0, this.getRowNumber(), compNumberByColumn, swapper);
        System.out.println("Table is sorted based on the number of column " + header.get(columnIndex));
        return true;
    }

    protected Swapper swapper = new Swapper() {
        @Override
        public void swap(int a, int b) {
            List<T> temp = cells.get(a);
            cells.set(a, cells.get(b));
            cells.set(b, temp);
        }
    };
    
    protected IntComparator compTextByColumn = new IntComparator() {
        @Override
        public int compare(int a, int b) {
            return cells.get(a).get(sortColumnIndex).toString().compareTo(cells.get(b).get(sortColumnIndex).toString());
        }
    };
    
    protected IntComparator compNumberByColumn = new IntComparator() {
        @Override
        public int compare(int a, int b) {
            double va = getCellAsDouble(a, sortColumnIndex);
            double vb = getCellAsDouble(b, sortColumnIndex);
            if (va < vb) return -1;
            else if (va == vb) return 0;
            return 1;
        }
    };

    @Override
    public void insertColumn(String columnName, int columnIndex, List<T> cList) {
        if (cList.size() != this.getRowNumber()) {
            System.out.println("The column to be inserted has different number of rows. Program quits");
            System.exit(1);
        }
        this.header.add(columnIndex, columnName);
        for (int i = 0; i < this.getRowNumber(); i++) {
            cells.get(i).add(columnIndex, (T)cList.get(i));
        }
        this.buildHIMap();
    }

    @Override
    public void removeColumn(int columnIndex) {
        this.header.remove(columnIndex);
        for (int i = 0; i < this.getRowNumber(); i++) {
            cells.get(i).remove(columnIndex);
        }
        this.buildHIMap();
    }
    
    @Override
    public void insertRow(int rowIndex, List<T> rList) {
        if (rList.size() != this.getColumnNumber()) {
            System.out.println("The row to be inserted has different number of columns. Program quits");
            System.exit(1);
        }
        cells.add(rowIndex, rList);
    }

    @Override
    public void removeRow(int rowIndex) {
        cells.remove(rowIndex);
    }

    @Override
    public void setCell(int rowIndex, int columnIndex, T c) {
        cells.get(rowIndex).set(columnIndex, c);
    }

    @Override
    public void setRow(int rowIndex, List<T> rList) {
        cells.set(rowIndex, rList);
    }

    @Override
    public void setColumn(int columnIndex, List<T> cList) {
        for (int i = 0; i < this.getRowNumber(); i++) {
            cells.get(i).set(columnIndex, cList.get(i));
        }
    }

    @Override
    public void addColumn(String columnName, List cList) {
        if (cList.size() != this.getRowNumber()) {
            System.out.println("The column to be added has different number of rows. Program quits");
            System.exit(1);
        }
        this.header.add(columnName);
        for (int i = 0; i < this.getRowNumber(); i++) {
            cells.get(i).add((T)cList.get(i));
        }
        this.hiMap.put(columnName, this.getColumnNumber()-1);
    }

    @Override
    public void addRow(List<T> rList) {
        if (rList.size() != this.getColumnNumber()) {
            System.out.println("The row to be added has different number of columns. Program quits");
            System.exit(1);
        }
        this.cells.add(rList);
    }

    @Override
    public List<T> getColumn(int columnIndex) {
        List<T> cList = new ArrayList();
        for (int i = 0; i < this.getRowNumber(); i++) {
            cList.add(cells.get(i).get(columnIndex));
        }
        return cList;
    }

    @Override
    public List<T> getRow(int rowIndex) {
        return cells.get(rowIndex);
    }
}
