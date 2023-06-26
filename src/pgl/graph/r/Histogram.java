/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package pgl.graph.r;

import java.io.File;
import rcaller.RCaller;
import rcaller.RCode;

/**
 * Plot histogram
 * @author Fei Lu
 */
public class Histogram extends Rgraphics {
    
    int breakNum = Integer.MIN_VALUE;
    double[] value = null;
    
    public Histogram () {
        this.yLab = "Frequency";
    }
    
    public Histogram (double[] value) {
        this();
        this.setValue(value);
    }
    
    public Histogram (int[] value) {
        this();
        this.setValue(value);
    }
    
    public void setValue (int[] intValue) {
        value = new double[intValue.length];
        for (int i = 0; i < value.length; i++) {
            value[i] = intValue[i];
        }
    }
    
    public void setValue (short[] shortValue) {
        value = new double[shortValue.length];
        for (int i = 0; i < value.length; i++) {
            value[i] = shortValue[i];
        }
    }
    
    public void setValue (byte[] byteValue) {
        value = new double[byteValue.length];
        for (int i = 0; i < value.length; i++) {
            value[i] = byteValue[i];
        }
    }
    
    public void setValue (double[] doubleValue) {
        this.value = doubleValue;
    }
    
    /**
     * Set breaks in histogram
     * @param breakNumber 
     */
    public void setBreakNumber (int breakNumber) {
        this.breakNum = breakNumber;
    }
    
    @Override
    public void showGraph () {
        RCaller caller = new RCaller();
        caller.setRscriptExecutable(this.rPath);
        RCode rCode = new RCode();
        rCode.addDoubleArray("x", value);
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
    
    /**
     * Save graph to pdf file
     * @param outfileS 
     */
    @Override
    public void saveGraph (String outfileS) {
        outfileS = outfileS.replaceAll("\\\\", "/");
        RCaller caller = new RCaller();
        caller.setRscriptExecutable(this.rPath);
        RCode rCode = new RCode();
        rCode.addDoubleArray("x", value);
        rCode.addRCode("pdf(\""+outfileS+"\")");
        rCode.addRCode(margin);
        rCode.addRCode(this.getPlotStatement());
        rCode.addRCode("dev.off()");
        caller.setRCode(rCode);
        caller.runOnly();
    }

    @Override
    public String getPlotStatement () {
        StringBuilder sb = new StringBuilder ();
        sb.append("hist(x,").append("col=\"").append(defaultColor).append("\",");
        if (breakNum != Integer.MIN_VALUE) sb.append("breaks=").append(breakNum).append(",");
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
        String statement = sb.toString();
        //System.out.println(statement);
        return statement;
    }
}
