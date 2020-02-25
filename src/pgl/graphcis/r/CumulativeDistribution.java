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
 *
 * @author fl262
 */
public class CumulativeDistribution extends Rgraphics {
    double[] value = null;
    boolean ifOnlyLine = false;
    
    public CumulativeDistribution (double[] value) {
        this.value = value;
    }
    
    public CumulativeDistribution (int[] intValue) {
        this.setValue(intValue);
    }
    
    public void setValue (int[] intValue) {
        value = new double[intValue.length];
        for (int i = 0; i < value.length; i++) {
            value[i] = intValue[i];
        }
    }
    
    public void setIfOnlyLine (boolean ifOnlyLine) {
        this.ifOnlyLine = ifOnlyLine;
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
            rCode.addRCode("d <- ecdf(x)");
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
        rCode.addRCode("d <- ecdf(x)");
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
        sb.append("col=\"").append(this.defaultColor).append("\",");
        sb.append("verticals=").append("T").append(",");
        if (this.ifOnlyLine) sb.append("do.points=").append("F").append(",");
        else sb.append("do.points=").append("T").append(",");
            
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
        String statement = sb.toString();
        //System.out.println(statement);
        return statement;
    }
}
