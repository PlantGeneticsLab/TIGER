/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.graph.tSaw;

import pgl.infra.table.TableInterface;
import java.io.IOException;

import tech.tablesaw.api.ColumnType;
import tech.tablesaw.api.Table;
import tech.tablesaw.api.TextColumn;
import tech.tablesaw.io.csv.CsvReadOptions;
import tech.tablesaw.io.csv.CsvWriteOptions;

/**
 *
 * @author feilu
 */
public class TablesawUtils {
    
    /**
     * Return a CsvReadOptions.Builder from a TSV file 
     * @param infileS
     * @return 
     */
    public static CsvReadOptions.Builder getTsvReadOptionBuilder (String infileS) {
        CsvReadOptions.Builder builder = CsvReadOptions.builder(infileS)
                                                        .separator('\t')
                                                        .header(true);
        return builder;
    }
    
    /**
     * Read a TSV file and return a Table
     * @param infileS
     * @return 
     */
    public static Table readTsv (String infileS) {
        CsvReadOptions.Builder builder = getTsvReadOptionBuilder (infileS);
        CsvReadOptions options = builder.build();
        try {
            return Table.read().csv(options);
        } catch (IOException e) {
            e.printStackTrace();
            return null;
        }
    }
    
    /**
     * Read a TSV file, set ColumnType, and return a Table
     * @param infileS
     * @param types e.g. ColumnType[] types = {LOCAL_DATE, INTEGER, FLOAT, FLOAT, CATEGORY}; 
     * e.g. ColumnType[] types = {SKIP, INTEGER, FLOAT, FLOAT, SKIP}
     * @return
     * @throws IOException 
     */
    public static Table readTsv (String infileS, ColumnType[] types) throws IOException {
        CsvReadOptions.Builder builder = getTsvReadOptionBuilder (infileS);
        builder.columnTypes(types);
        CsvReadOptions options = builder.build();
        return Table.read().csv(options);
    }
    
    /**
     * Read a CSV file and return a Table
     * @param infileS
     * @return
     * @throws IOException 
     */
    public static Table readCsv (String infileS) throws IOException {
        return Table.read().file(infileS);
    }
    
    /**
     * Return a Table from TableInferface
     * @param ti
     * @return 
     */
    public static Table getTable (TableInterface<String> ti) {
        Table t = Table.create();
        for (int i = 0; i < ti.getColumnNumber(); i++) {
            TextColumn c = TextColumn.create(ti.getHeader().get(i), ti.getColumn(i));
            t.addColumns(c);
        }
        return t;
    }

    /**
     * Write a Table to a TSV file
     * @param t
     * @param outfileS
     */
    public static void writeTsv (Table t, String outfileS) {
        CsvWriteOptions builder;
        try {
            builder = CsvWriteOptions.builder(outfileS).separator('\t').header(true).build();
            t.write().csv(builder);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
