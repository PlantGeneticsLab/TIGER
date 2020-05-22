/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package pgl.graph;

import com.itextpdf.awt.DefaultFontMapper;
import com.itextpdf.awt.PdfGraphics2D;
import com.itextpdf.text.Document;

import com.itextpdf.text.Rectangle;
import com.itextpdf.text.pdf.PdfContentByte;
import com.itextpdf.text.pdf.PdfWriter;
import pgl.infra.table.TableInterface;
import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Ellipse2D;
import java.io.FileOutputStream;

/**
 *
 * @author fl262
 */
public class HeatMap {
    double[][] matrix = null;
    String[] rowName = null;
    String[] columnName = null;
    int cellX = 50;
    int cellY = 25;
    int xMargin = 200;
    int yMargin = 200;
    double lowValue = 0;
    double highValue = 1;
    int[] lowValueColor = {255, 255, 255};
    int[] highValueColor = {255, 0, 0};
    int lengendNumber = 10;
    String legendName = "";
    
    public HeatMap (double[][] matrix, String[] rowName, String[] columnName) {
        this.matrix = matrix;
        this.rowName = rowName;
        this.columnName = columnName;
    }
    
    public HeatMap (TableInterface t) {
        matrix = new double[t.getRowNumber()][t.getColumnNumber()-1];
        rowName = new String[t.getRowNumber()];
        columnName = new String[t.getColumnNumber()-1];
        for (int i = 0; i < t.getRowNumber(); i++) {
            rowName[i] = t.getCellAsString(i, 0);
            for (int j = 0; j < matrix[0].length; j++) {
                t.getCellAsDouble(i, j+1);
                matrix[i][j] = t.getCellAsDouble(i, j+1);
            }
        }
        for (int i = 0; i < columnName.length; i++) {
            columnName[i] = t.getColumnName(i+1);
        }
    }
    
    public void setCellSize (int cellX, int cellY) {
        this.cellX = cellX;
        this.cellY = cellY;
    }
    
    public void setCapValues (double lowValue, double highValue) {
        this.lowValue = lowValue;
        this.highValue = highValue;
    }
    
    public void setLowValueColor (int r, int g, int b) {
        lowValueColor = new int[3];
        lowValueColor[0] = r;
        lowValueColor[1] = g;
        lowValueColor[2] = b;
    }
    
    public void setLowValueColor (int r, int g, int b, int a) {
        lowValueColor = new int[4];
        lowValueColor[0] = r;
        lowValueColor[1] = g;
        lowValueColor[2] = b;
        lowValueColor[3] = a;
    }
    
    public void setHighValueColor (int r, int g, int b, int a) {
        highValueColor = new int[3];
        highValueColor[0] = r;
        highValueColor[1] = g;
        highValueColor[2] = b;
        highValueColor[3] = a;
    }
    
    public void setLegendName (String legendName) {
        this.legendName = legendName;
    }
    
    public void saveGraph (String outfileS) {
        this.draw(outfileS);
    }
    
    private Color getColorByValue (double v) {
        int[] valueColor = new int[lowValueColor.length];
        double ratio = (v-lowValue)/(highValue-lowValue);
        for (int i = 0; i < valueColor.length; i++) {
            valueColor[i] = lowValueColor[i] + (int)(ratio*(highValueColor[i] - lowValueColor[i]));
        }
        if (valueColor.length == 3) {
            return new Color(valueColor[0], valueColor[1], valueColor[2]);
        }
        else {
            return new Color(valueColor[0], valueColor[1], valueColor[2], valueColor[3]);
        }
    }
    
    private void draw (String outfileS) {
        int fontSize = 8;
        int width = 2*xMargin + cellX*this.matrix[0].length;
        int height = 2*yMargin + cellY*this.matrix.length;
        try {
            Document pd = new Document(new Rectangle(width,height));
            PdfWriter pw;
            pw = PdfWriter.getInstance(pd, new FileOutputStream (outfileS));
            pd.open();
            PdfContentByte canvas = pw.getDirectContent();
            DefaultFontMapper mapper = new DefaultFontMapper();
            PdfGraphics2D g2d = new PdfGraphics2D(canvas, width, height, mapper);
            Font f = new Font("Arial", Font.PLAIN , fontSize);
            g2d.setFont(f);
            for (int i = 0; i < matrix.length; i++) {
                g2d.setColor(Color.black);
                g2d.drawString(rowName[i], xMargin-80, yMargin+i*cellY+cellY/2+fontSize/2);
                int y = yMargin+i*cellY;
                for (int j = 0; j < matrix[0].length; j++) {
                    int x = xMargin+j*cellX;
                    Color co = this.getColorByValue(matrix[i][j]);
                    g2d.setColor(co);
                    g2d.fillRect(x, y, cellX, cellY);
                    g2d.setColor(Color.CYAN);
                    //g2d.setColor(Color.white);
                    g2d.drawRect(x, y, cellX, cellY);
                }
                
            }
            double diameter = ((matrix.length * cellY)*0.5)/this.lengendNumber;
            diameter = (int)(xMargin*0.8);
            int lastX = 0;
            int lastY = 0;
            int yStep = (cellY*matrix.length-(int)diameter)/this.lengendNumber;
            f = new Font("Arial", Font.PLAIN , 20);
            g2d.setFont(f);
            for (int i = 0; i < this.lengendNumber+1; i++) {
                int x = xMargin+matrix[0].length*cellX+(int)(xMargin*0.1);
                //int y = yMargin + (int)((diameter*((double)1.4))*i) + cellY;
                int y = yMargin + yStep*i;
                Ellipse2D.Double s = new Ellipse2D.Double(x,y,diameter,diameter);
                double value = ((this.highValue-this.lowValue)/this.lengendNumber) * i;
                g2d.setColor(this.getColorByValue(value));
                g2d.fill(s);
                g2d.setColor(Color.black);
                g2d.draw(s);
                x += cellX/2;
                y += cellY/2;
                
                g2d.drawString(String.format("%.2g%n", value), x, y);
                lastX = x;
                lastY = y;
            }
            g2d.setColor(Color.black);
            lastX = lastX - cellX/2;
            lastY = lastY+(int)(diameter*((double)1.4));
            g2d.drawString(this.legendName, lastX, lastY);
            f = new Font("Arial", Font.PLAIN , fontSize);
            g2d.setFont(f);
            for (int i = 0; i < matrix[0].length; i++) {
                int x =  xMargin+i*cellX+cellX/2;
                int y = yMargin-cellY/2;
                g2d.rotate(-1.57, x, y);
                g2d.drawString(columnName[i], x, y);
                g2d.rotate(1.57, x, y);
            }
            g2d.dispose();
            pd.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
