/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package pgl.graphcis.r;

import java.io.File;
import rcaller.RCaller;
import rcaller.RCode;

/**
 * Plot scatter plot
 * @author fl262
 */
public class ScatterPlot extends Rgraphics {
    double[] xValues;
    double[] yValues;
    boolean ifTrendLine = false;
    public ScatterPlot () {}
    
    public ScatterPlot (double[] x, double[] y) {
        this.xValues = x;
        this.yValues = y;
    }
    
    public ScatterPlot (int[] x, int[] y) {
        this.setXValues(x);
        this.setYValues(y);
    }
    
    public void setXValues (double[] x) {
        this.xValues = x;
    }
    
    public void setXValues (int[] x) {
        xValues = new double[x.length];
        for (int i = 0; i < x.length; i++) {
            xValues[i] = x[i];            
        }
    }
    
    public void setYValues (double[] y) {
        this.yValues = y;
    }
    
    public void setYValues (int[] y) {
        yValues = new double[y.length];
        for (int i = 0; i < y.length; i++) {
            yValues[i] = y[i];
        }
    }
    
    public void addTrendLine () {
        this.ifTrendLine = true;
    }
    
    public void setColor (int r, int g, int b, int a) {
        this.defaultColor = "rgb("+String.valueOf(r)+","+String.valueOf(g)+","+String.valueOf(b)+","+String.valueOf(a)+",maxColorValue=255)";
    }
    
    @Override
    public void showGraph() {
        RCaller caller = new RCaller();
        caller.setRscriptExecutable(this.rPath);
        RCode rCode = new RCode();
        rCode.addDoubleArray("x", xValues);
        rCode.addDoubleArray("y", yValues);
        try {
            File plotFile = rCode.startPlot();
            rCode.addRCode(margin);
            rCode.addRCode(this.getPlotStatement());
            caller.setRCode(rCode);
            caller.runOnly();
            rCode.endPlot();
            rCode.showPlot(plotFile);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    @Override
    public void saveGraph(String outfileS) {
        outfileS = outfileS.replaceAll("\\\\", "/");
        RCaller caller = new RCaller();
        caller.setRscriptExecutable(this.rPath);
        RCode rCode = new RCode();
        rCode.addDoubleArray("x", xValues);
        rCode.addDoubleArray("y", yValues);
        rCode.addRCode("pdf(\""+outfileS+"\")");
        rCode.addRCode(margin);
        rCode.addRCode(this.getPlotStatement());
        rCode.addRCode("dev.off()");
        caller.setRCode(rCode);
        caller.runOnly();
    }

    @Override
    public String getPlotStatement() {
        StringBuilder sb = new StringBuilder ();
        sb.append("plot(x, y, ");
        if (defaultColor.startsWith("rgb")) {
            sb.append("col=").append(defaultColor).append(",");
        }
        else {
            sb.append("col=\"").append(defaultColor).append("\",");
        }
        sb.append("pch=").append(this.pch).append(",");
        sb.append("main=\"").append(title).append("\",");
        sb.append("xlab=\"").append(xLab).append("\",");
        sb.append("ylab=\"").append(yLab).append("\",");
        sb.append("cex.main=").append(titleSize).append(",");
        sb.append("cex.lab=").append(labSize).append(",");
        sb.append("cex.axis=").append(axisSize).append(",");
        sb.append("font.main=").append(titleFont).append(",");
        if (!Double.isNaN(this.xLimLow)) {
            sb.append("xlim = c(").append(this.xLimLow).append(",").append(this.xLimHigh).append("),");
        }
        if (!Double.isNaN(this.yLimLow)) {
            sb.append("ylim = c(").append(this.yLimLow).append(",").append(this.yLimHigh).append("),");
        }
        if (sb.charAt(sb.length()-1) == ',') {
            sb.deleteCharAt(sb.length()-1);
        }
        sb.append(")");
        if (this.ifTrendLine) {
            sb.append("\n").append("abline(lm(y ~ x))");
        }
        String statement = sb.toString();
        //System.out.println(statement);
        return statement;
    }
}
