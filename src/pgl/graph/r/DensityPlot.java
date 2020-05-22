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
 *
 * @author fl262
 */
public class DensityPlot extends Rgraphics {
    double[] value = null;
    int smoothN = 512;
    public DensityPlot (double[] value) {
        this.value = value;
    }
    
    public DensityPlot (int[] intValue) {
        this.setValue(intValue);
    }
    
    public void setSmoothN (int n) {
        smoothN = n;
    }
    
    public void setValue (int[] intValue) {
        value = new double[intValue.length];
        for (int i = 0; i < value.length; i++) {
            value[i] = intValue[i];
        }
    }

    @Override
    public void showGraph() {
        RCaller caller = new RCaller();
        caller.setRscriptExecutable(this.rPath);
        RCode rCode = new RCode();
        rCode.addDoubleArray("x", value);
        try {
            File plotFile = rCode.startPlot();
            rCode.addRCode(margin);
            StringBuilder sb = new StringBuilder();
            sb.append("d <- density(x, n = ").append(smoothN).append(")");
            rCode.addRCode(sb.toString());
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
        rCode.addDoubleArray("x", value);
        StringBuilder sb = new StringBuilder();
        sb.append("d <- density(x, n = ").append(smoothN).append(")");
        rCode.addRCode(sb.toString());
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
        sb.append("plot(d,");
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
        sb.append(")\n");
        sb.append("polygon(d, col=").append("\"").append(this.defaultColor).append("\", border=\"").append(this.defaultColor).append("\")");
        String statement = sb.toString();
        //System.out.println(statement);
        return statement;
    }
}
