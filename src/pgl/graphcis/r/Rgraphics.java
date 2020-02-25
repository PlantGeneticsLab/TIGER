/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package pgl.graphcis.r;

import pgl.graphcis.GraphicsUtils;

/**
 *
 * @author Fei Lu
 */
public abstract class Rgraphics implements RgraphicsImpl {
    public static String RPath = "/usr/local/bin/Rscript";
    public static String[] defaultColorArray = {"red", "blue", "gold", "cyan", "orange", "chocolate", "green1", "coral"};
    String rPath = RPath;
    String defaultColor = "red";
    String title = "Graph title";
    String xLab = "Variable x";
    String yLab = "Variavle y";
    //String margin = "par(mar=c(5.1, 5.1, 4.1, 3.1))";
    String margin = "par(mar=c(8, 5.1, 4.1, 3.1))";
    int titleFont = 1;
    double titleSize = 1;
    double labSize = 1;
    double axisSize = 1;
    int pch = 1;
    String font = "sans";
    double xLimLow = Double.NaN;
    double xLimHigh = Double.NaN;
    double yLimLow = Double.NaN;
    double yLimHigh = Double.NaN;
    double xLegend = Double.NaN;
    double yLegend = Double.NaN;
    
    @Override
    public void setRPath (String rPath) {
        this.rPath = rPath;
    }
    
    @Override
    public void setColor (String color) {
        this.defaultColor = color;
    }
    
    public void setTitle (String title) {
        this.title = title;
    }
    
    @Override
    public void setXLab (String xLab) {
        this.xLab = xLab;
    }
    
    @Override
    public void setYLab (String yLab) {
        this.yLab = yLab;
    }
    
    @Override
    public void setTitleSize (double size) {
        this.titleSize = size;
    }
    
    @Override
    public void setYLim (double low, double high) {
        this.yLimLow = low;
        this.yLimHigh = high;
    }
    
    @Override
    public void setXLim (double low, double high) {
        this.xLimLow = low;
        this.xLimHigh = high;
    }
    
    @Override
    public void setLegendPosition (double x, double y) {
        this.xLegend = x;
        this.yLegend = y;
    }
    
    /**
     * 1: ordinary, 2: bold, 3: italics
     * @param n 
     */
    @Override
    public void setTitleFont (int n) {
        this.titleFont = n;
    }
    
    @Override
    public void setLableSize (double labSize) {
        this.labSize = labSize;
    }
    
    @Override
    public void setAxisSize (double axisSize) {
        this.axisSize = axisSize;
    }
    
    @Override
    public void setPlottingCharacter (int pch) {
        this.pch = pch;
    }
    
    @Override
    public void setPaperMode () {
        titleFont = 1;
        titleSize = 1;
        labSize = 1;
        axisSize = 1;
    }
  
    @Override
    public void setSlideMode () {
        titleFont = 1;
        titleSize = 1.5;
        labSize = 1.5;
        axisSize = 1.5;
    }
    
    @Override
    public void setLargerSize () {
        titleSize = titleSize*1.1;
        labSize = labSize*1.1;
        axisSize = axisSize*1.1;
    }
    
}
