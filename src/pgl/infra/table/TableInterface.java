/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.infra.table;

import java.util.List;
import pgl.infra.utils.IOFileFormat;

/**
 * Interface holding basic methods of a table
 * @author feilu
 * @param <T> The type of data in a cell
 */
public interface TableInterface <T> {
    
    /**
     * Return a cell
     * @param rowIndex
     * @param columnIndex
     * @return 
     */
    public T getCell (int rowIndex, int columnIndex);
    
    /**
     * Return a string value of a cell.
     * @param rowIndex
     * @param columnIndex
     * @return 
     */
    public String getCellAsString (int rowIndex, int columnIndex);
    
    /**
     * Return a Double value of a cell, if the cell is convertible to double. Return null if is not.
     * @param rowIndex
     * @param columnIndex
     * @return 
     */
    public Double getCellAsDouble (int rowIndex, int columnIndex);
    
    /**
     * Return an Integer value of a cell, if the cell is convertible to int. Return null if is not. The Integer value may loose its precision during conversion.
     * @param rowIndex
     * @param columnIndex
     * @return 
     */
    public Integer getCellAsInteger (int rowIndex, int columnIndex);
    
    /**
     * Set a cell
     * @param rowIndex
     * @param columnIndex
     * @param c 
     */
    public void setCell (int rowIndex, int columnIndex, T c);
    
    /**
     * Return the header of table
     * @return 
     */
    public List<String> getHeader();
    
    /**
     * Return a column name
     * @param columnIndex
     * @return 
     */
    public String getColumnName (int columnIndex);
    
    /**
     * Return values as a list from a certain column
     * @param columnIndex
     * @return 
     */
    public List<T> getColumn (int columnIndex);
    
    /**
     * Return a double array from a column
     * @param columnIndex
     * @return 
     */
    public double[] getColumnAsDoubleArray (int columnIndex);

    /**
     * Return a int array from a column
     * @param columnIndex
     * @return
     */
    public int[] getColumnAsIntArray (int columnIndex);

    /**
     * Return a values as list from a certain row
     * @param rowIndex
     * @return 
     */
    public List<T> getRow (int rowIndex);
    
    /**
     * Return the index of header name in table 
     * @param columnName
     * @return 
     */
    public int getColumnIndex (String columnName);
    
    /**
     * Return the total number of column
     * @return 
     */
    public int getColumnNumber ();
    
    /**
     * Return the total number of rows
     * @return 
     */
    public int getRowNumber();
    
    /**
     * Add a column to the table at the end
     * @param columnName
     * @param cList 
     */
    public void addColumn (String columnName, List<T> cList);
    
    /**
     * Add a row to the table at the end
     * @param rList 
     */
    public void addRow (List<T> rList);
    
    /**
     * Insert a column in the table
     * @param columnName
     * @param columnIndex
     * @param cList 
     */
    public void insertColumn (String columnName, int columnIndex, List<T> cList);
    
    /**
     * Remove a column from the table
     * @param columnIndex 
     */
    public void removeColumn (int columnIndex);
    
    /**
     * Remove a column from the table
     * @param columnName 
     */
    public void removeColumn (String columnName);
    
    /**
     * Insert a row to the table
     * @param rowIndex
     * @param rList 
     */
    public void insertRow (int rowIndex, List<T> rList);
    
    /**
     * Remove a row from the table
     * @param rowIndex 
     */
    public void removeRow (int rowIndex);
    
    /**
     * Set values of a row
     * @param rowIndex
     * @param rList 
     */
    public void setRow (int rowIndex, List<T> rList);
    
    /**
     * Set values of a column
     * @param columnIndex
     * @param cList 
     */
    public void setColumn (int columnIndex, List<T> cList);
    
    /**
     * Sort the table by text of a column, natural order
     * @param columnIndex 
     */
    public void sortAsText (int columnIndex);
    
    /**
     * Sort the table by the number of a column, numerical order
     * @param columnIndex
     * @return false if the column cannot be sorted as number
     */
    public boolean sortAsNumber (int columnIndex);
    
    /**
     * Sort the table by text of a column, natural order
     * @param columnName 
     */
    public void sortAsText (String columnName);
    
    /**
     * Sort the table by the number of a column, numerical order
     * @param columnName
     * @return 
     */
    public boolean sortAsNumber (String columnName);
    
    /**
     * Write a text table
     * @param outfileS
     * @param format 
     */
    public void writeTextTable (String outfileS, IOFileFormat format);
    
    /**
     * Write a text table with selected rows
     * @param outfileS
     * @param format
     * @param ifOut 
     */
    public void writeTextTable (String outfileS, IOFileFormat format, boolean[] ifOut);
    
}
