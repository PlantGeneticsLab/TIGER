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
public class GenomeTrendPlotTwoYAxis {
    int chromNum;
    int[] chromLength;
    ArrayList<int[][]> leftPositionList = new ArrayList();
    ArrayList<double[][]> leftValueList = new ArrayList();
    ArrayList<int[][]> rightPositionList = new ArrayList();
    ArrayList<double[][]> rightValueList = new ArrayList();
    double leftVMax;
    double leftVMin;
    double leftVMaxTemp = Double.NaN;
    double leftVMinTemp = Double.NaN;
    double rightVMax;
    double rightVMin;
    double rightVMaxTemp = Double.NaN;
    double rightVMinTemp = Double.NaN;
    String title = "";
     
    int width = 1800;
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
    ArrayList<int[][]> leftTPosList = new ArrayList();
    ArrayList<int[][]> leftTValueList = new ArrayList();
    ArrayList<int[][]> rightTPosList = new ArrayList();
    ArrayList<int[][]> rightTValueList = new ArrayList();
    int[][] tXCoor;
    
    
    public GenomeTrendPlotTwoYAxis (int[] chromLength, int[][] leftPosition, double[][] leftValue, int[][] rightPosition, double[][] rightValue, String title) {
        this.chromNum = chromLength.length;
        this.chromLength = chromLength;
        this.leftPositionList.add(leftPosition);
        this.leftValueList.add(leftValue);
        this.rightPositionList.add(rightPosition);
        this.rightValueList.add(rightValue);
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
        double range = this.leftVMax-this.leftVMin;
        double ratioV = this.trackHeight/range;
        for (int i = 0; i < leftValueList.size(); i++) {
            int[][] tPos = new int[chromNum][];
            int[][] tValue = new int[chromNum][];
            for (int j = 0; j < chromNum; j++) {
                tPos[j] = new int[leftPositionList.get(i)[j].length];
                tValue[j] = new int[leftValueList.get(i)[j].length];
                for (int k = 0; k < tPos[j].length; k++) {
                    tPos[j][k] = (int)(leftPositionList.get(i)[j][k]*ratio);
                    tValue[j][k] = (int)((leftValueList.get(i)[j][k]-this.leftVMin)*ratioV);
                }
            }
            this.leftTPosList.add(tPos);
            this.leftTValueList.add(tValue);
        }
        range = this.rightVMax-this.rightVMin;
        ratioV = this.trackHeight/range;
        for (int i = 0; i < rightValueList.size(); i++) {
            int[][] tPos = new int[chromNum][];
            int[][] tValue = new int[chromNum][];
            for (int j = 0; j < chromNum; j++) {
                tPos[j] = new int[rightPositionList.get(i)[j].length];
                tValue[j] = new int[rightValueList.get(i)[j].length];
                for (int k = 0; k < tPos[j].length; k++) {
                    tPos[j][k] = (int)(rightPositionList.get(i)[j][k]*ratio);
                    tValue[j][k] = (int)((rightValueList.get(i)[j][k]-this.rightVMin)*ratioV);
                }
            }
            this.rightTPosList.add(tPos);
            this.rightTValueList.add(tValue);
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
        for (int i = 0; i < leftValueList.size(); i++) {
            for (int j = 0; j < leftValueList.get(i).length; j++) {
                for (int k = 0; k < leftValueList.get(i)[j].length; k++) {
                    if (leftValueList.get(i)[j][k] < min) min = leftValueList.get(i)[j][k];
                    if (leftValueList.get(i)[j][k] > max) max = leftValueList.get(i)[j][k];
                }
            }
        }
        this.leftVMax = max;
        this.leftVMin = min;
        if (!Double.isNaN(leftVMaxTemp)) leftVMax = leftVMaxTemp;
        if (!Double.isNaN(leftVMinTemp)) leftVMin = leftVMinTemp;
        max = Double.MIN_VALUE;
        min = Double.MAX_VALUE;
        for (int i = 0; i < rightValueList.size(); i++) {
            for (int j = 0; j < rightValueList.get(i).length; j++) {
                for (int k = 0; k < rightValueList.get(i)[j].length; k++) {
                    if (rightValueList.get(i)[j][k] < min) min = rightValueList.get(i)[j][k];
                    if (rightValueList.get(i)[j][k] > max) max = rightValueList.get(i)[j][k];
                }
            }
        }
        this.rightVMax = max;
        this.rightVMin = min;
        if (!Double.isNaN(rightVMaxTemp)) rightVMax = rightVMaxTemp;
        if (!Double.isNaN(rightVMinTemp)) rightVMin = rightVMinTemp;
    }
    
    public void addLeftData (int[][] position, double[][] value) {
        this.leftPositionList.add(position);
        this.leftValueList.add(value);
    }
    
    public void addRightData (int[][] position, double[][] value) {
        this.rightPositionList.add(position);
        this.rightValueList.add(value);
    }
    
    public void setTrackHeight (int height)  {
        this.trackHeight = height;
    }   
    
    public void setIntervalHeight (int height)  {
        this.intervalHeight = height;
    }
    
    public void setLeftLowValueLim (double lowValue) {
        this.leftVMinTemp = lowValue;
    }
    
    public void setLeftHighValueLim (double highValue)  {
        this.leftVMaxTemp = highValue;
    }
    
    public void setRightLowValueLim (double lowValue) {
        this.rightVMinTemp = lowValue;
    }
    
    public void setRightHighValueLim (double highValue)  {
        this.rightVMaxTemp = highValue;
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
                g2d.setColor(Color.LIGHT_GRAY);
                g2d.fillRoundRect(xMargin, currentY, tChrLength[i], barHeight, 15, 15);
                if (ifWithDot) {
                    for (int j = 0; j < leftTPosList.size(); j++) {
                        g2d.setColor(GraphicsUtils.rainbowColor[j*2]);
                        for (int k = 0; k < leftTPosList.get(j)[i].length; k++) {
                            int y = currentY-leftTValueList.get(j)[i][k];
                            g2d.fillOval(leftTPosList.get(j)[i][k]+xMargin-dotSize/2, y-dotSize/2, dotSize, dotSize);
                        }
                    }
                    for (int j = 0; j < rightTPosList.size(); j++) {
                        g2d.setColor(GraphicsUtils.rainbowColor[j*2+1]);
                        for (int k = 0; k < rightTPosList.get(j)[i].length; k++) {
                            int y = currentY-rightTValueList.get(j)[i][k];
                            g2d.fillOval(rightTPosList.get(j)[i][k]+xMargin-dotSize/2, y-dotSize/2, dotSize, dotSize);
                        }
                    }
                }
                if (ifWithLine) {
                    for (int j = 0; j < leftTPosList.size(); j++) {
                        g2d.setColor(GraphicsUtils.rainbowColor[j*2]);
                        for (int k = 0; k < leftTPosList.get(j)[i].length-1; k++) {
                            int y = currentY-leftTValueList.get(j)[i][k];
                            int nextY = currentY-leftTValueList.get(j)[i][k+1];
                            g2d.drawLine(leftTPosList.get(j)[i][k]+xMargin, y, leftTPosList.get(j)[i][k+1]+xMargin, nextY);
                        }
                    }
                    for (int j = 0; j < rightTPosList.size(); j++) {
                        g2d.setColor(GraphicsUtils.rainbowColor[j*2+1]);
                        for (int k = 0; k < rightTPosList.get(j)[i].length-1; k++) {
                            int y = currentY-rightTValueList.get(j)[i][k];
                            int nextY = currentY-rightTValueList.get(j)[i][k+1];
                            g2d.drawLine(rightTPosList.get(j)[i][k]+xMargin, y, rightTPosList.get(j)[i][k+1]+xMargin, nextY);
                        }
                    }
                }
                g2d.setColor(Color.black);
                g2d.drawLine(xMargin, currentY, xMargin, currentY-yMargin);
                
                int x = xMargin-(int)(0.15*xMargin);
                double part = this.trackHeight/this.corMarkNum;
                double increment = (this.leftVMax-this.leftVMin)/this.corMarkNum;
                for (int j = 0; j < this.corMarkNum; j++) {
                    int y = currentY-(int)(j*part);
                    g2d.drawLine(xMargin, y, xMargin-(int)(0.02*xMargin), y);
                    g2d.drawString(String.format("%.2f", leftVMin+j*increment), x, y);
                }
                x = xMargin-(int)(0.5*xMargin);
                g2d.drawString("Chromosome "+String.valueOf(i+1), x, currentY-yMargin/2);
                
                g2d.drawLine(xMargin+this.tChrLength[i], currentY, xMargin+this.tChrLength[i], currentY-yMargin);
                x = xMargin+this.tChrLength[i]+(int)(0.04*xMargin);
                part = this.trackHeight/this.corMarkNum;
                increment = (this.rightVMax-this.rightVMin)/this.corMarkNum;
                for (int j = 0; j < this.corMarkNum; j++) {
                    int y = currentY-(int)(j*part);
                    g2d.drawLine(xMargin+this.tChrLength[i], y, xMargin+this.tChrLength[i]+(int)(0.02*xMargin), y);
                    g2d.drawString(String.format("%.2f", rightVMin+j*increment), x, y);
                }
                
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
