/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package pgl.graph.r;

import java.io.File;
import org.apache.commons.lang3.math.NumberUtils;
import rcaller.RCaller;
import rcaller.RCode;

/**
 *
 * @author Fei Lu
 */
public class LineChart extends Rgraphics {
    double[] x;
    double[][] ys;
    String[] variableNames;
    double[] xRange = new double[2];
    double[] yRange = new double[2];
    String lineType = "b";
    boolean ifSmooth = false;
    
    public LineChart (double[] x, double[] y, String variableName) {
        this.x = x;
        ys = new double[1][];
        ys[0] = y;
        variableNames = new String[1];
        variableNames[0] = variableName;
        this.creatRange();
    }
    
    public LineChart (double[] x, double[][] ys, String[] variableNames) {
        this.x = x;
        this.ys = ys;
        this.variableNames = variableNames;
        this.creatRange();
    }
    
    private void creatRange () {
        xRange[0] = NumberUtils.min(x); xRange[1] = NumberUtils.max(x);
        double[] maxs = new double[ys.length];
        double[] mins = new double[ys.length];
        for (int i = 0; i < maxs.length; i++) {
            maxs[i] = NumberUtils.max(ys[i]);
            mins[i] = NumberUtils.min(ys[i]);
        }
        yRange[0] = NumberUtils.min(mins); yRange[1] = NumberUtils.max(maxs);
    }
    
    public void setIfSmoothLine (boolean value) {
        this.ifSmooth = value;
    }
    
    public void showOnlyDot () {
        this.lineType = "p";
    }
    
    public void showOnlyLine () {
        this.lineType = "l";
    }
    
    public void showDotAndLine () {
        this.lineType = "b";
    }
    @Override
    public void showGraph() {
        RCaller caller = new RCaller();
        caller.setRscriptExecutable(this.rPath);
        RCode rCode = new RCode();
        for (int i = 0; i < ys.length; i++) {
            rCode.addDoubleArray("y"+String.valueOf(i), ys[i]);
        }
        rCode.addDoubleArray("x", x);
        rCode.addDoubleArray("xRange", this.xRange);
        rCode.addDoubleArray("yRange", this.yRange);
        rCode.addStringArray("variableNames", variableNames);
        rCode.addRCode(this.getPlotStatement());
        
        System.out.println(rCode.getCode().toString());
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
    public String getPlotStatement() {
        StringBuilder sb = new StringBuilder ();
        sb.append("plot(xRange, yRange, type = \"n\",").append("col=\"").append(defaultColor).append("\",");
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
        sb.append("colors <- rainbow(length(variableNames))\n");
        sb.append("linetype <- c(1:length(variableNames))\n");
        sb.append("plotchar <- seq(18, 18+length(variableNames), 1)\n");
        for (int i = 0; i < ys.length; i++) {
            if (this.ifSmooth) {
                sb.append("lines(spline(x, y"+String.valueOf(i)).append("), type=\"").append(this.lineType).append("\","+"lwd=1.5, ");
            }
            else {
                sb.append("lines(x, y"+String.valueOf(i)).append(", type=\"").append(this.lineType).append("\","+"lwd=1.5, ");
            }
            
            sb.append("lty=linetype["+String.valueOf(i+1)+"], col=colors["+String.valueOf(i+1)+"], pch=plotchar["+String.valueOf(i+1)+"])\n");
        }
        if (!Double.isNaN(this.xLegend)) {
            sb.append("legend(").append(this.xLegend).append(", ").append(this.yLegend).append(", ").append("variableNames, cex=0.8, col=colors, pch=plotchar, lty=linetype, title=\"Legend\")");
        }
        else sb.append("legend(xRange[1], yRange[2], variableNames, cex=0.8, col=colors, pch=plotchar, lty=linetype, title=\"Legend\")");
        String statement = sb.toString();
        //System.out.println(statement);
        return statement;
    }

    @Override
    public void saveGraph(String outfileS) {
        outfileS = outfileS.replaceAll("\\\\", "/");
        RCaller caller = new RCaller();
        caller.setRscriptExecutable(this.rPath);
        RCode rCode = new RCode();
        rCode.addRCode("pdf(\""+outfileS+"\")");
        for (int i = 0; i < ys.length; i++) {
            rCode.addDoubleArray("y"+String.valueOf(i), ys[i]);
        }
        rCode.addDoubleArray("x", x);
        rCode.addDoubleArray("xRange", this.xRange);
        rCode.addDoubleArray("yRange", this.yRange);
        rCode.addStringArray("variableNames", variableNames);
        rCode.addRCode(this.getPlotStatement());
        rCode.addRCode("dev.off()");
        caller.setRCode(rCode);
        caller.runOnly();
    }
}
