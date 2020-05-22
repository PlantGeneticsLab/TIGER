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
import java.awt.Color;
import java.io.FileOutputStream;
import java.util.ArrayList;

/**
 *
 * @author Fei Lu
 */
public class GenomeTrendPlot {
    int chromNum;
    int[] chromLength;
    ArrayList<int[][]> positionList = new ArrayList();
    ArrayList<double[][]> valueList = new ArrayList();
    double vMax;
    double vMin;
    double vMaxTemp = Double.NaN;
    double vMinTemp = Double.NaN;
    String title = "";
     
    int width = 1600;
    int height = 0;
    int maxBarLength = 1200;
    int xMargin = 300;
    int yMargin = 100;
    int trackHeight = 120;
    int intervalHeight = 25;
    int barHeight = 10;
    int corMarkNum = 4;
    int dotSize = 8;
    boolean ifWithLine = true;
    boolean ifWithDot = true;
    int xCoorInterval = 5000000;
    int xUnit = 1000000;
    
    int[] tChrLength;
    ArrayList<int[][]> tPosList = new ArrayList();
    ArrayList<int[][]> tValueList = new ArrayList();
    int[][] tXCoor;
    
    
    public GenomeTrendPlot (int[] chromLength, int[][] position, double[][] value, String title) {
        this.chromNum = chromLength.length;
        this.chromLength = chromLength;
        this.positionList.add(position);
        this.valueList.add(value);
        this.title = title;
    }
    
    public void transform () {
        this.findExtremeValue();
        int maxChrLength = Integer.MIN_VALUE;
        for (int i = 0; i < chromLength.length; i++) {
            if (chromLength[i] > maxChrLength) maxChrLength = chromLength[i];
        }
        double ratio = (double)maxBarLength/maxChrLength;
        tChrLength = new int[chromNum];
        double range = this.vMax-this.vMin;
        double ratioV = this.trackHeight/range;
        for (int i = 0; i < valueList.size(); i++) {
            int[][] tPos = new int[chromNum][];
            int[][] tValue = new int[chromNum][];
            for (int j = 0; j < chromNum; j++) {
                tPos[j] = new int[positionList.get(i)[j].length];
                tValue[j] = new int[valueList.get(i)[j].length];
                for (int k = 0; k < tPos[j].length; k++) {
                    tPos[j][k] = (int)(positionList.get(i)[j][k]*ratio);
                    tValue[j][k] = (int)((valueList.get(i)[j][k]-this.vMin)*ratioV);
                }
            }
            this.tPosList.add(tPos);
            this.tValueList.add(tValue);
        }
        this.tXCoor = new int[chromNum][];
        for (int i = 0; i < chromNum; i++) {
            tChrLength[i] = (int)(chromLength[i]*ratio);
            tXCoor[i] = new int[chromLength[i]/this.xCoorInterval+1];
            for (int j = 0; j < tXCoor[i].length; j++) {
                tXCoor[i][j] = (int)(j*this.xCoorInterval*ratio);
            }
        }
        this.height = 2*yMargin+chromNum*trackHeight+(chromNum-1)*intervalHeight;
    }
    
    private void findExtremeValue () {
        double max = Double.MIN_VALUE;
        double min = Double.MAX_VALUE;
        for (int i = 0; i < valueList.size(); i++) {
            for (int j = 0; j < valueList.get(i).length; j++) {
                for (int k = 0; k < valueList.get(i)[j].length; k++) {
                    if (valueList.get(i)[j][k] < min) min = valueList.get(i)[j][k];
                    if (valueList.get(i)[j][k] > max) max = valueList.get(i)[j][k];
                }
            }
        }
        this.vMax = max;
        this.vMin = min;
        if (!Double.isNaN(vMaxTemp)) vMax = vMaxTemp;
        if (!Double.isNaN(vMinTemp)) vMin = vMinTemp;
    }
    
    public void addData (int[][] position, double[][] value) {
        this.positionList.add(position);
        this.valueList.add(value);
    }
    
    public void setTrackHeight (int height)  {
        this.trackHeight = height;
    }   
    
    public void setIntervalHeight (int height)  {
        this.intervalHeight = height;
    }
    
    public void setLowValueLim (double lowValue) {
        this.vMinTemp = lowValue;
    }
    
    public void setHighValueLim (double highValue)  {
        this.vMaxTemp = highValue;
    }
    
    public void setIfWithDot (boolean ifwithDot) {
        this.ifWithDot = ifwithDot;
    }
    
    public void setIfWithLine (boolean ifwithLine) {
        this.ifWithLine = ifwithLine;
    }
    
    public void saveGragh (String outfileS) {
        this.transform();
        try {
            Document pd = new Document(new Rectangle(width,height));
            PdfWriter pw;
            pw = PdfWriter.getInstance(pd, new FileOutputStream (outfileS));
            pd.open();
            PdfContentByte canvas = pw.getDirectContent();
            DefaultFontMapper mapper = new DefaultFontMapper();
            PdfGraphics2D g2d = new PdfGraphics2D(canvas, width, height, mapper);
            int currentY = this.yMargin;
            for (int i = 0; i < chromLength.length; i++) {
                currentY+=this.trackHeight;
                g2d.setColor(Color.cyan);
                g2d.fillRoundRect(xMargin, currentY, tChrLength[i], barHeight, 15, 15);
                
                if (ifWithDot) {
                    for (int j = 0; j < tPosList.size(); j++) {
                        g2d.setColor(GraphicsUtils.rainbowColor[j]);
                        for (int k = 0; k < tPosList.get(j)[i].length; k++) {
                            int y = currentY-tValueList.get(j)[i][k];
                            g2d.fillOval(tPosList.get(j)[i][k]+xMargin-dotSize/2, y-dotSize/2, dotSize, dotSize);
                        }
                    }
                }
                if (ifWithLine) {
                    for (int j = 0; j < tPosList.size(); j++) {
                        g2d.setColor(GraphicsUtils.rainbowColor[j]);
                        for (int k = 0; k < tPosList.get(j)[i].length-1; k++) {
                            int y = currentY-tValueList.get(j)[i][k];
                            int nextY = currentY-tValueList.get(j)[i][k+1];
                            g2d.drawLine(tPosList.get(j)[i][k]+xMargin, y, tPosList.get(j)[i][k+1]+xMargin, nextY);
                        }
                    }
                }
                g2d.setColor(Color.black);
                g2d.drawLine(xMargin, currentY, xMargin, currentY-yMargin);
                int x = xMargin-(int)(0.15*xMargin);
                double part = this.trackHeight/this.corMarkNum;
                double increment = (this.vMax-this.vMin)/this.corMarkNum;
                for (int j = 0; j < this.corMarkNum; j++) {
                    int y = currentY-(int)(j*part);
                    g2d.drawLine(xMargin, y, xMargin-(int)(0.02*xMargin), y);
                    g2d.drawString(String.format("%.2f", vMin+j*increment), x, y);
                }
                x = xMargin-(int)(0.5*xMargin);
                g2d.drawString("Chromosome "+String.valueOf(i+1), x, currentY-yMargin/2);
                
                for (int j = 0; j < this.tXCoor[i].length; j++) {
                    x = xMargin+this.tXCoor[i][j];
                    g2d.drawLine(x, currentY, x, currentY+this.barHeight);
                    String s = String.valueOf(j*this.xCoorInterval/this.xUnit)+"Mb";
                    g2d.drawString(s, x, currentY+2*this.barHeight);
                }
                currentY+=this.intervalHeight;
            }
            g2d.drawString(title, xMargin, yMargin/2);
            g2d.dispose();
            pd.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
